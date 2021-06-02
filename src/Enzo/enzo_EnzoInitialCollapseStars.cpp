// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialShuCollapse.cpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     2021-03-22
/// @brief    [\ref Enzo] Initial conditions for a uniform density sphere of stars
///


#include "cello.hpp"
#include "enzo.hpp"

EnzoInitialCollapseStars::EnzoInitialCollapseStars
(const EnzoConfig * enzo_config) throw()
  :Initial(enzo_config->initial_cycle, enzo_config->initial_time)
{

  const EnzoUnits * enzo_units = enzo::units();
  
  centre_[0] = enzo_config->initial_collapse_stars_centre[0];
  centre_[1] = enzo_config->initial_collapse_stars_centre[1];
  centre_[2] = enzo_config->initial_collapse_stars_centre[2];

  drift_velocity_[0] = enzo_config->initial_collapse_stars_drift_velocity[0];
  drift_velocity_[1] = enzo_config->initial_collapse_stars_drift_velocity[1];
  drift_velocity_[2] = enzo_config->initial_collapse_stars_drift_velocity[2];
    
  truncation_radius_ = enzo_config->initial_collapse_stars_truncation_radius;
  total_mass_ = enzo_config->initial_collapse_stars_total_mass;
  particle_fraction_ = enzo_config->initial_collapse_stars_particle_fraction;
  random_seed_ = enzo_config->initial_collapse_stars_random_seed;
}

void EnzoInitialCollapseStars::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  PUParray(p,centre_,3);
  PUParray(p,drift_velocity_,3);
  p | truncation_radius_;
  p | total_mass_;
  p | particle_fraction_;
}

void EnzoInitialCollapseStars::enforce_block
( Block * block, Hierarchy * hierarchy ) throw()

{
  if (!block->is_leaf()) return;
  const EnzoConfig * enzo_config = enzo::config();
  ASSERT("EnzoInitialCollapseStars",
	 "Block does not exist",
	 block != NULL);

  // Check we have sensible parameters
  // First check if rank = 3
  ASSERT("EnzoInitialCollapseStars::EnzoInitialCollapseStars()",
	 "Cannot run CollapseStars with rank < 3",cello::rank() == 3);

  // Check we are in unigrid mode
  ASSERT("EnzoInitialCollapseStars::EnzoInitialCollapseStars()",
	 "CollapseStars requires unigrid mode (Adapt : max_level = 0). ",
	 enzo_config->mesh_max_level == 0);


  // Check if we have periodic boundary conditions
  int periodic_x,periodic_y,periodic_z;
  hierarchy->get_periodicity(&periodic_x,&periodic_y,&periodic_z);
  ASSERT("EnzoInitialCollapseStars::EnzoInitialCollapseStars()",
	 "CollapseStars must have periodic boundary conditions",
	 periodic_x && periodic_y && periodic_z);

  // Check if the truncation radius is less than half the domain size
  double dxm,dxp,dym,dyp,dzm,dzp;
  hierarchy->lower(&dxm,&dym,&dzm);
  hierarchy->upper(&dxp,&dyp,&dzp);
  double domain_width_x = dxp - dxm;
  double domain_width_y = dyp - dym;
  double domain_width_z = dzp - dzm;
  
  
  ASSERT("EnzoInitialCollapseStars::EnzoInitialCollapseStars()",
	   "Truncation radius must be less than half the domain width in all "
	   "dimensions",
	   (truncation_radius_ < 0.5 * domain_width_x) &&
	   (truncation_radius_ < 0.5 * domain_width_y) &&
	   (truncation_radius_ < 0.5 * domain_width_z));

  // TODO: Check required fields?

  double folded_centre_position[3];
  hierarchy->get_folded_position(centre_, folded_centre_position);

  centre_[0] = folded_centre_position[0];
  centre_[1] = folded_centre_position[1];
  centre_[2] = folded_centre_position[2];

  Field field = block->data()->field();

  // Get Field parameters
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);
  const int mx = nx + 2*gx;
  const int my = ny + 2*gy;
  const int mz = nz + 2*gz;
  const int m = mx*my*mz;

  // Block extents
  double bxm,bym,bzm;
  double bxp,byp,bzp;
  block->data()->lower(&bxm,&bym,&bzm);
  block->data()->upper(&bxp,&byp,&bzp);
  double hx,hy,hz;
  field.cell_width(bxm,bxp,&hx,
		   bym,byp,&hy,
		   bzm,bzp,&hz);

  // Get pointers to fields
  enzo_float * dt = (enzo_float *) field.values ("density_total");
  enzo_float * po = (enzo_float *) field.values ("potential");
  enzo_float * ax = (enzo_float *) field.values ("acceleration_x");
  enzo_float * ay = (enzo_float *) field.values ("acceleration_y");
  enzo_float * az = (enzo_float *) field.values ("acceleration_z");
  enzo_float *  x = (enzo_float *) field.values ("X");
  enzo_float *  b = (enzo_float *) field.values ("B");

  // For some of the fields, we initialise their values to zero
  std::fill_n(dt,m,0.0);
  std::fill_n(po,m,0.0);
  std::fill_n(ax,m,0.0);
  std::fill_n(ay,m,0.0);
  std::fill_n(az,m,0.0);
  std::fill_n(x,m,0.0);
  std::fill_n(b,m,0.0);

  ParticleDescr * particle_descr = cello::particle_descr();
  Particle particle              = block->data()->particle();
  
  // Attribute indices
  const int it   = particle.type_index("star");
  const int ia_m = particle.attribute_index (it, "mass");
  const int ia_x = particle.attribute_index (it, "x");
  const int ia_y = particle.attribute_index (it, "y");
  const int ia_z = particle.attribute_index (it, "z");
  const int ia_vx = particle.attribute_index (it, "vx");
  const int ia_vy = particle.attribute_index (it, "vy");
  const int ia_vz = particle.attribute_index (it, "vz");
  const int ia_copy = particle.attribute_index (it, "is_copy");
  
  // Attribrute stride lengths
  const int dm   = particle.stride(it, ia_m);
  const int dp   = particle.stride(it, ia_x);
  const int dv   = particle.stride(it, ia_vx);
  const int dloc = particle.stride(it, ia_copy);
  
  /// Initialise pointers for particle attribute arrays
  enzo_float * pmass = 0;
  enzo_float * px   = 0;
  enzo_float * py   = 0;
  enzo_float * pz   = 0;
  enzo_float * pvx  = 0;
  enzo_float * pvy  = 0;
  enzo_float * pvz  = 0;
  int64_t * is_copy = 0;
      
      
  // Set random seed
  srand(random_seed_);


  // Compute the mass of each star particle

  const double total_volume = 4./3. * cello::pi * truncation_radius_ *
    truncation_radius_ * truncation_radius_;

  const double cell_volume = hx*hy*hz;

  const double particle_mass = total_mass_ * cell_volume /
    (particle_fraction_ * total_volume);

  // Loop over active cells

  for (int iz = gz; iz < gz + nz; iz++){
    const double z = bzm + (iz + 0.5)*hz;
    for (int iy = gy; iy < gy + ny; iy++){
      const double y = bym + (iy + 0.5)*hy;
      for (int ix = gx; ix < gx + nx; ix++){
	const double x = bxm + (ix + 0.5)*hx;

	// Check if within the truncation radius
	double cell_pos[3] = {x,y,z};
	double npi[3];
	hierarchy->get_nearest_periodic_image(cell_pos,centre_,npi);
	const double r2 =
	  (npi[0] - centre_[0]) * (npi[0] - centre_[0]) +
	  (npi[1] - centre_[1]) * (npi[1] - centre_[1]) +
	  (npi[2] - centre_[2]) * (npi[2] - centre_[2]);

	if (r2 < truncation_radius_ * truncation_radius_){
	  double rnum = (double(rand())) / (double(RAND_MAX));
	  if (rnum < particle_fraction_){
	    // insert particle
	    int ib,ipp  = 0;
	    const int new_particle = particle.insert_particles(it, 1);
	    particle.index(new_particle,&ib,&ipp);
	    
	    // Get pointers to particle attribute arrays
	    pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
	    px    = (enzo_float *) particle.attribute_array(it, ia_x, ib);
	    py    = (enzo_float *) particle.attribute_array(it, ia_y, ib);
	    pz    = (enzo_float *) particle.attribute_array(it, ia_z, ib);
	    pvx   = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
	    pvy   = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
	    pvz   = (enzo_float *) particle.attribute_array(it, ia_vz, ib);
	    is_copy   = (int64_t *) particle.attribute_array(it, ia_copy, ib);
	    
	    // Now assign values to attributes
	    pmass[ipp*dm] = particle_mass;
	    px[ipp*dp] = x;
	    py[ipp*dp] = y;
	    pz[ipp*dp] = z;
	    pvx[ipp*dp] = drift_velocity_[0];
	    pvy[ipp*dp] = drift_velocity_[1];
	    pvz[ipp*dp] = drift_velocity_[2];
	    is_copy[ipp*dloc] = 0;
	    
	  } // if rnum < particle_fraction
	} // if r2 < truncation_radius * truncation_radius
	
      } // ix
    } // iy
  } // iz
  return;
}

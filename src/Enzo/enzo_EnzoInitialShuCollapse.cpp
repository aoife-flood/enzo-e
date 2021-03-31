// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialShuCollapse.cpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     2021-03-22
/// @brief    Implementation of the Shu Collapse problem (Shu 1977)
///


#include "cello.hpp"
#include "enzo.hpp"

EnzoInitialShuCollapse::EnzoInitialShuCollapse
(const EnzoConfig * enzo_config) throw()
  :Initial(enzo_config->initial_cycle, enzo_config->initial_time)
{

  const EnzoUnits * enzo_units = enzo::units();
  
  centre_[0] = enzo_config->initial_shu_collapse_centre[0];
  centre_[1] = enzo_config->initial_shu_collapse_centre[1];
  centre_[2] = enzo_config->initial_shu_collapse_centre[2];

  drift_velocity_[0] = enzo_config->initial_shu_collapse_drift_velocity[0];
  drift_velocity_[1] = enzo_config->initial_shu_collapse_drift_velocity[1];
  drift_velocity_[2] = enzo_config->initial_shu_collapse_drift_velocity[2];
    
  truncation_radius_ = enzo_config->initial_shu_collapse_truncation_radius;
  sound_speed_ = enzo_config->initial_shu_collapse_sound_speed;
  instability_parameter_ = enzo_config->initial_shu_collapse_instability_parameter;
  central_particle_ = enzo_config->initial_shu_collapse_central_particle;
  central_particle_mass_ = enzo_config->initial_shu_collapse_central_particle_mass;
  
  // Get physics attributes
  gamma_ = enzo_config->field_gamma;
  ggm1_ = gamma_ * (gamma_ - 1.0);
  grav_constant_internal_units_ =
    cello::grav_constant * enzo_units->mass() *
    enzo_units->time() * enzo_units->time() /
    (enzo_units->length() * enzo_units->length() *
     enzo_units->length());

}

void EnzoInitialShuCollapse::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  PUParray(p,centre_,3);
  PUParray(p,drift_velocity_,3);
  p | truncation_radius_;
  p | sound_speed_;
  p | instability_parameter_;
  p | central_particle_;
  
}

void EnzoInitialShuCollapse::enforce_block
( Block * block, Hierarchy * hierarchy ) throw()

{
  if (!block->is_leaf()) return;
  const EnzoConfig * enzo_config = enzo::config();
  ASSERT("EnzoInitialShuCollapse",
	 "Block does not exist",
	 block != NULL);

  // Check we have sensible parameters
  // First check if rank = 3
  ASSERT("EnzoInitialShuCollapse::EnzoInitialShuCollapse()",
	 "Cannot run Shu Collapse with rank < 3",cello::rank() == 3);

  // Check we are in unigrid mode
  ASSERT("EnzoInitialShuCollapse::EnzoInitialShuCollapse()",
	 "Shu Collapse requires unigrid mode (Adapt : max_level = 0). ",
	 enzo_config->mesh_max_level == 0);


  // Check if we have periodic boundary conditions
  int px,py,pz;
  hierarchy->get_periodicity(&px,&py,&pz);
  ASSERT("EnzoInitialShuCollapse::EnzoInitialShuCollapse()",
	 "Shu Collapse must have periodic boundary conditions",px && py && pz);

  // Check if the truncation radius is less than half the domain size
  double dxm,dxp,dym,dyp,dzm,dzp;
  hierarchy->lower(&dxm,&dym,&dzm);
  hierarchy->upper(&dxp,&dyp,&dzp);
  double domain_width_x = dxp - dxm;
  double domain_width_y = dyp - dym;
  double domain_width_z = dzp - dzm;
  
  
  ASSERT("EnzoInitialShuCollapse::EnzoInitialShuCollapse()",
	   "Truncation radius must be less than half the domain width in all "
	   "dimensions",
	   (truncation_radius_ < 0.5 * domain_width_x) &&
	   (truncation_radius_ < 0.5 * domain_width_y) &&
	   (truncation_radius_ < 0.5 * domain_width_z));

  // If we have a central particle, need to have star formation turned off
  if (central_particle_){
    const int num_methods = enzo_config->num_method;
    bool star_maker_on = false;
    for (int i = 0; i < num_methods; i++){
      if (enzo_config->method_list[i] == "star_maker"){
	ERROR("EnzoInitialShuCollapse::EnzoInitialShuCollapse()",
	      "If central_particle = true, then cannot have star_maker method");
      }
    }
  }

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
  enzo_float *  d = (enzo_float *) field.values ("density");
  enzo_float * dt = (enzo_float *) field.values ("density_total");
  enzo_float *  p = (enzo_float *) field.values ("pressure");
  enzo_float * po = (enzo_float *) field.values ("potential");
  enzo_float * te = (enzo_float *) field.values ("total_energy");
  enzo_float * ie = (enzo_float *) field.values ("internal_energy");
  enzo_float * ax = (enzo_float *) field.values ("acceleration_x");
  enzo_float * ay = (enzo_float *) field.values ("acceleration_y");
  enzo_float * az = (enzo_float *) field.values ("acceleration_z");
  enzo_float * vx = (enzo_float *) field.values ("velocity_x");
  enzo_float * vy = (enzo_float *) field.values ("velocity_y");
  enzo_float * vz = (enzo_float *) field.values ("velocity_z");
  enzo_float *  x = (enzo_float *) field.values ("X");
  enzo_float *  b = (enzo_float *) field.values ("B");

  // For some of the fields, we initialise their values to zero
  std::fill_n(dt,m,0.0);
  std::fill_n(p,m,0.0);
  std::fill_n(po,m,0.0);
  std::fill_n(ax,m,0.0);
  std::fill_n(ay,m,0.0);
  std::fill_n(az,m,0.0);
  std::fill_n(x,m,0.0);
  std::fill_n(b,m,0.0);

  // Set specific internal energy
  const enzo_float ie_value = sound_speed_ * sound_speed_ / ggm1_;
  std::fill_n(ie,m,ie_value);

  // Set velocity
  std::fill_n(vx,m,drift_velocity_[0]);
  std::fill_n(vy,m,drift_velocity_[1]);
  std::fill_n(vz,m,drift_velocity_[2]);

  // Set total energy
  const enzo_float ke = 0.5 * (drift_velocity_[0] * drift_velocity_[0] +
			       drift_velocity_[1] * drift_velocity_[1] +
			       drift_velocity_[2] * drift_velocity_[2]);

  const enzo_float te_value = ke + ie_value;
  std::fill_n(te,m,te_value);

  // Now to initialise the density field
  const double density_profile_factor =
    instability_parameter_ *
    sound_speed_ * sound_speed_ /
    (4.0 * cello::pi * grav_constant_internal_units_);

  for (int iz = 0; iz < mz ; iz++){
    const double z = bzm + (iz - gz + 0.5)*hz;
    for (int iy = 0; iy < my; iy++){
      const double y = bym + (iy - gy + 0.5)*hy;
      for (int ix = 0; ix < mx; ix++){
	const double x = bxm + (ix - gx + 0.5)*hx;
	double cell_pos[3] = {x,y,z};
	double npi[3];
	hierarchy->get_nearest_periodic_image(cell_pos,centre_,npi);
	const double r2 =
	  (npi[0] - centre_[0]) * (npi[0] - centre_[0]) +
	  (npi[1] - centre_[1]) * (npi[1] - centre_[1]) +
	  (npi[2] - centre_[2]) * (npi[2] - centre_[2]);
	const int i = INDEX(ix,iy,iz,mx,my);
	d[i] =
	  (r2 < truncation_radius_ * truncation_radius_) ?
	  density_profile_factor / r2 :
	  density_profile_factor / (truncation_radius_ * truncation_radius_);	      
      } //ix
    } //iy
  } //iz

  // If central_particle is true and collapse centre is in this block, we
  // add a particle at the collapse centre
  if (central_particle_ &&
      block->check_position_in_block(centre_[0],centre_[1],centre_[2]))
    {
    
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
      const int ia_loc = particle.attribute_index (it, "is_local");

      // Attribrute stride lengths
      const int dm   = particle.stride(it, ia_m);
      const int dp   = particle.stride(it, ia_x);
      const int dv   = particle.stride(it, ia_vx);
      const int dloc = particle.stride(it, ia_loc);
  
      /// Initialise pointers for particle attribute arrays
      enzo_float * pmass = 0;
      enzo_float * px   = 0;
      enzo_float * py   = 0;
      enzo_float * pz   = 0;
      enzo_float * pvx  = 0;
      enzo_float * pvy  = 0;
      enzo_float * pvz  = 0;
      int64_t * is_local = 0;
     
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
      is_local   = (int64_t *) particle.attribute_array(it, ia_loc, ib);

      // Now assign values to attributes
      pmass[ipp*dm] = central_particle_mass_;
      px[ipp*dp] = centre_[0];
      py[ipp*dp] = centre_[1];
      pz[ipp*dp] = centre_[2];
      pvx[ipp*dp] = drift_velocity_[0];
      pvy[ipp*dp] = drift_velocity_[1];
      pvz[ipp*dp] = drift_velocity_[2];
      is_local[ipp*dloc] = 1;

    } // Is there are central particle to place in this block?
  return;
}

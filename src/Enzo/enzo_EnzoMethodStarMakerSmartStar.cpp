/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodStarMakerSmartStar.cpp
/// @author     John Regan (john.regan@mu.ie)
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date
/// @brief      Functionality to allow for SMS, PopIII, PopII and Black Hole Formation
///
///     Derived star maker class that actually makes stars. This is partially
///     adapted after the star_maker_ssn method from Enzo and from the active
///     particle type SmartStars from Enzo

#include "cello.hpp"
#include "enzo.hpp"
#include <time.h>

#define DEBUG_SMARTSTARS
#define SCALAR 1

//-------------------------------------------------------------------

EnzoMethodStarMakerSmartStar::EnzoMethodStarMakerSmartStar
()
  : EnzoMethodStarMaker() 
{
  const EnzoConfig * enzo_config = enzo::config();
  check_potential_minimum_ = enzo_config->method_smart_stars_check_potential_minimum;
}

//-------------------------------------------------------------------

void EnzoMethodStarMakerSmartStar::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  EnzoMethodStarMaker::pup(p); // call parent class pup

  return;
}

//------------------------------------------------------------------

void EnzoMethodStarMakerSmartStar::compute ( Block *block) throw()
{

  // Loop through the grid and check star formation criteria
  // There are a number of star formation criteria that must be
  // fulfilled before star formation is triggered

  // TODO: Should check if we are on the maximum refinement level
  if (! block->is_leaf()){
    block->compute_done();
    return;
  }

  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();
  
  Particle particle = enzo_block->data()->particle();
  Field field = enzo_block->data()->field();

  double dx, dy, dz;
  block->cell_width(&dx, &dy, &dz);
  
  double lx, ly, lz;
  double ux, uy, uz;
  block->lower(&lx,&ly,&lz);
  block->upper(&ux,&uy,&uz);

  // declare particle position arrays
  //  default particle type is "star", but this will default
  //  to subclass particle_type
  const int it   = particle.type_index (this->particle_type());
  const int ia_m = particle.attribute_index (it, "mass");
  const int ia_pm = particle.attribute_index (it, "prevmass");
  const int ia_x = particle.attribute_index (it, "x");
  const int ia_y = particle.attribute_index (it, "y");
  const int ia_z = particle.attribute_index (it, "z");
  const int ia_vx = particle.attribute_index (it, "vx");
  const int ia_vy = particle.attribute_index (it, "vy");
  const int ia_vz = particle.attribute_index (it, "vz");
  const int ia_ax = particle.attribute_index(it,"ax");
  const int ia_metal = particle.attribute_index (it, "metal_fraction");
  const int ia_to    = particle.attribute_index (it, "creation_time");
  const int ia_l     = particle.attribute_index (it, "lifetime");
  const int ia_timeindex = particle.attribute_index (it, "timeindex");
  const int ia_class = particle.attribute_index (it, "class");
  const int ia_accrate = particle.attribute_index (it, "accretion_rate");
  const int ia_accrate_time = particle.attribute_index (it, "accretion_rate_time");

  /// pointers for particle attribute arrays
  enzo_float * pmass = 0;
  enzo_float * prevmass = 0;
  enzo_float * px   = 0;
  enzo_float * py   = 0;
  enzo_float * pz   = 0;
  enzo_float * pvx  = 0;
  enzo_float * pvy  = 0;
  enzo_float * pvz  = 0;
  enzo_float * pmetal = 0;
  enzo_float * pform  = 0;
  enzo_float * plifetime = 0;
  int        * ptimeindex = 0;
  int        * pclass     = 0;
  enzo_float * paccrate = 0;
  enzo_float * paccrate_time = 0;

  // obtain the particle stride length
  // (Not sure if this is right, isn't there a separate stride length for each attribrute?)
  const int ps = particle.stride(it, ia_m);
  int ipp; // Particle index
  int ib; // Batch index
  int rank = cello::rank();

  int gx,gy,gz;
  field.ghost_depth (0, &gx, &gy, &gz);

  int nx, ny, nz;
  field.size ( &nx, &ny, &nz);

  int mx = nx + 2*gx;
  int my = ny + 2*gy;
  int mz = nz + 2*gz;

  // get pointers to field values
  enzo_float * density     = (enzo_float *) field.values("density");
  enzo_float * temperature = (enzo_float *) field.values("temperature");
  enzo_float * potential   = (enzo_float *) field.values("potential");
  enzo_float * velocity_x = (rank >= 1) ?
    (enzo_float *)field.values("velocity_x") : NULL;
  enzo_float * velocity_y = (rank >= 2) ?
    (enzo_float *)field.values("velocity_y") : NULL;
  enzo_float * velocity_z = (rank >= 3) ?
    (enzo_float *)field.values("velocity_z") : NULL;
  enzo_float * metal = field.is_field("metal_density") ?
    (enzo_float *) field.values("metal_density") : NULL;
  
  EnzoComputeTemperature compute_temperature
    (enzo_config->ppm_density_floor,
     enzo_config->ppm_temperature_floor,
     enzo_config->ppm_mol_weight,
     enzo_config->physics_cosmology);

  compute_temperature.compute(enzo_block);
  int count = 0; // Number of particles formed
  // iterate over all cells (not including ghost zones)
  for (int iz=gz; iz<nz+gz; iz++){
    for (int iy=gy; iy<ny+gy; iy++){
      for (int ix=gx; ix<nx+gx; ix++){

        int i = ix + mx*(iy + my*iz);

        // need to compute this better for Grackle fields (on to-do list)
        double mean_particle_mass = enzo_config->ppm_mol_weight * cello::mass_hydrogen;
	//        double ndens = rho_cgs / mean_particle_mass;
	double rho_cgs = density[i] * enzo_units->density();
        double mass_in_solar_masses  = density[i] *dx*dy*dz * enzo_units->mass() / cello::mass_solar;
	double jeans_density;
        //
        // Apply the criteria for star formation
        //

	//(i) The first criteria is that the local gas density must
	//    exceed the jeans density
	if(! this->check_jeans_density(temperature[i],
				       dx*enzo_units->length(),
				       density[i],&jeans_density)) continue;

	// (ii) The second criteria is that the velocity divergence is negative
	if (! this->check_velocity_divergence(velocity_x, velocity_y,
					      velocity_z, i,
                                              1, my, my*mz)) continue;

	//(iii) Is cell the gravitational minimum over a jeans length
	double cellpos[3];
	cellpos[0] = lx + (ix - gx + 0.5) * dx;
	cellpos[1] = ly + (iy - gy + 0.5) * dy;
	cellpos[2] = lz + (iz - gz + 0.5) * dz;
	
	if(! this->check_potential_minimum(enzo_block, cellpos, potential[i],
					   temperature[i], rho_cgs)) continue;

        

        // now create a star particle
	count++; 
        int my_particle = particle.insert_particles(it, 1);
        particle.index(my_particle, &ib, &ipp);
        int io = ipp; // ipp*ps Stefan: Not sure about this
        pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
	prevmass = (enzo_float *) particle.attribute_array(it, ia_pm, ib);
	// particle mass is mass here
        pmass[io] = (density[i] - jeans_density) * dx * dy * dz;

	// set particle position to be at centre of cell
        px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
        py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
        pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);
        px[io] = lx + (ix - gx + 0.5) * dx;
        py[io] = ly + (iy - gy + 0.5) * dy;
        pz[io] = lz + (iz - gz + 0.5) * dz;

	// set particle velocity to be cell velocity
        pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
        pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
        pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib);
        pvx[io] = velocity_x[i];
        if (velocity_y) pvy[io] = velocity_y[i];
        if (velocity_z) pvz[io] = velocity_z[i];

        // finalize attributes
        plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib);
        pform     = (enzo_float *) particle.attribute_array(it, ia_to, ib);
        pform[io]     =  enzo_block->time();   // formation time
	// lifetime is hard-coded here, maybe should be a parameter
        plifetime[io] =  10.0 * cello::Myr_s / enzo_units->time() ;
	paccrate = (enzo_float *) particle.attribute_array(it, ia_accrate, ib);
        paccrate_time = (enzo_float *) particle.attribute_array(it, ia_accrate_time,
								ib);
	ptimeindex = (int *) particle.attribute_array(it, ia_timeindex, ib);
	ptimeindex[io] = 0;
	pclass = (int *) particle.attribute_array(it, ia_class, ib);
	pclass[io] = SMS;  //see star_type enum in _enzo.hpp
	*(&paccrate[io]) = 0.0; // I don't count formation as accretion per se
	*(&paccrate_time[io]) = pform[io];

	/* Starting point of accretion rate pointer for this particle */
	paccrate = &paccrate[io];
	paccrate_time = &paccrate_time[io];
	paccrate += sizeof(enzo_float);
	paccrate_time += sizeof(enzo_float);

	/* Now initialise the array from this offset point onwards */
	/* I don't think this is necessary. I think we can set attribrutes
	   to be arrays */
	for(int k = 1; k < SS_NTIMES; k++, paccrate++, paccrate_time++) {
	  *paccrate = k;
	  *paccrate_time = k;
	}
	
        if (metal){
          pmetal     = (enzo_float *) particle.attribute_array(it, ia_metal, ib);
          pmetal[io] = metal[i] / density[i];
        }

	/* Specify that particle is local to the block */
	//is_local = (int64_t *) particle.attribute_array(it,ia_local,ib);
	//is_local[io] = 1;
	
        // Remove mass from grid and rescale fraction fields
	density[i] = jeans_density;
        if (density[i] < 0){
          CkPrintf("Smartstar: density index mass: %g %i %g\n",
                   density[i],i,mass_in_solar_masses);
          ERROR("EnzoMethodStarMakerSmartStar::compute()",
                "Negative densities in star formation");
        }

        /* rescale tracer fields to maintain constant mass fraction
	   with the corresponding new density... */
        rescale_densities(enzo_block, i, 1.0 - jeans_density/density[i]);
      }
    }
  } // end loop iz

#ifdef DEBUG_SMARTSTARS
  if (count > 0){
    CkPrintf("SmartStars: Number of particles formed on Block %s = %i \n",
	     block->name().c_str(),count);
  }
#endif

  block->compute_done();
  return;
} // end compute


/*
 * Compute the gravitational minimum within 1 Jeans length of the cell
 * In order to pass this test the potential of the local cell
 * must be at the minimum of the potential.
 */
int EnzoMethodStarMakerSmartStar::check_potential_minimum(EnzoBlock * enzo_block,
						     const double * cellpos,
						     const double local_potential,
						     const double temperature,
						     const double rho_cgs
)	     
{
  
  if (!check_potential_minimum_)
    return 1;
  Field field = enzo_block->data()->field();
  const EnzoConfig * enzo_config = enzo::config();
  const EnzoUnits * enzo_units = enzo::units();
  int gx,gy,gz;
  field.ghost_depth (0, &gx, &gy, &gz);

  int mx, my, mz;
  field.dimensions (0, &mx, &my, &mz);

  int nx, ny, nz;
  field.size ( &nx, &ny, &nz);

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;
  double dx, dy, dz;
  enzo_block->cell_width(&dx, &dy, &dz);
  double lx, ly, lz;
  enzo_block->lower(&lx,&ly,&lz);
  const double mean_particle_mass = cello::mass_hydrogen*enzo_config->ppm_mol_weight;
  const double jeans_length = sqrt(gamma_ * cello::kboltz * temperature /
			      (cello::pi * cello::grav_constant
			       * mean_particle_mass * rho_cgs)) /
                                 enzo_units->length();
  
  double pot_min = 1e10;
  enzo_float * g = (enzo_float *) field.values("potential");
  /*  
   * Now need to determine the gravitational potential in all 
   * cells within a jeans_length
   */

  /* Here we only need to loop over a subset of the cells, 
     and can break the loop when a value larger than the test
     value is found. */

  /* ISSUE: Jeans Length can 'extend' beyond the ghost zone */

  /* However, if we check for whether Jeans length is resolved by less 
     than N cells, and require the ghost depth to be greater than or equal
     to N, and check this condition after the Jeans density condition,
     we are fine */

  /* According to Fedderath, we check within a 'control volume' */
   for (int iz = 0; iz < ngz; iz++){
      for (int iy = 0; iy < ngy; iy++){
        for (int ix = 0; ix < ngx; ix++){
          int i = INDEX(ix,iy,iz,ngx,ngy);

	  double posx = lx + (ix - gx + 0.5) * dx;
	  double posy = ly + (iy - gy + 0.5) * dy;
	  double posz = lz + (iz - gz + 0.5) * dz;
	  double dist = sqrt(pow((cellpos[0] - posx), 2.0) +
			     pow((cellpos[1] - posy), 2.0) +
			     pow((cellpos[2] - posz), 2.0));
	  if(dist <= jeans_length) {
	    if(g[i] < pot_min)
	      pot_min = g[i];
	  }
	}
      }
   }
   return (local_potential <= pot_min);
   
}


  
/*
   Defaults to parent class timestep if nothing declared here
double EnzoMethodStarMakerSmartStar::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}
*/

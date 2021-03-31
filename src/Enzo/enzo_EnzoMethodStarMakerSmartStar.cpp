/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodStarMakerSmartStar.cpp
/// @author     John Regan (john.regan@mu.ie)
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date
/// @brief      Functionality to allow for SMS, PopIII, PopII and
///             Black Hole Formation
///
///     Derived star maker class that actually makes stars. This is partially
///     adapted after the star_maker_ssn method from Enzo and from the active
///     particle type SmartStars from Enzo

#include "cello.hpp"
#include "enzo.hpp"
#include <time.h>

//#define DEBUG_SMARTSTARS
#define SCALAR 1

//-------------------------------------------------------------------

EnzoMethodStarMakerSmartStar::EnzoMethodStarMakerSmartStar
()
  : EnzoMethodStarMaker() 
{
  /// Check we are using sensible parameters
  const EnzoConfig * enzo_config = enzo::config();

  /// At the moment, can only run with SmartStars in unigrid mode
  ASSERT("EnzoMethodStarMakerSmartStar::EnzoMethodStarMakerSmartStar()",
	 "SmartStars requires unigrid mode (Adapt : max_level = 0) since " 
         "the Jeans length refinement criterion is not yet implemented",
	 enzo_config->mesh_max_level == 0);  
  ASSERT("EnzoMethodStarMakerSmartStar::EnzoMethodStarMakerSmartStar()",
	 "control_volume_cells_min_ must be at least 1",
	 control_volume_cells_min_ > 0);
  ASSERT("EnzoMethodStarMakerSmartStar::EnzoMethodStarMakerSmartStar()",
	 "control_volume_cells_max_ must be greater than or equal to "
	 "control_volume_cells_min_",
	 control_volume_cells_max_ >= control_volume_cells_min_);

  /// Control volume must lie within ghost zones
  const int * ghost_depth = enzo_config->field_ghost_depth;
  const int min_ghost_depth = std::min(ghost_depth[0],
				       std::min(ghost_depth[1],ghost_depth[2]));

  ASSERT("EnzoMethodStarMakerSmartStar::EnzoMethodStarMakerSmartStar()",
	 "control_volume_cells_max_ must be less than or equal to the "
	 "ghost zone depth",
	 control_volume_cells_max_ <= min_ghost_depth);

  /// We should generalise these checks to the case where we have different
  /// number of resultion elements in each dimension

  srand(time(NULL)); // need random number generator for later


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
  double cell_volume = dx * dy * dz;
  double mean_cell_width = cbrt(cell_volume);
  
  double lx, ly, lz;
  double ux, uy, uz;
  block->lower(&lx,&ly,&lz);
  block->upper(&ux,&uy,&uz);

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
  const int ia_mf = particle.attribute_index (it, "metal_fraction");
  const int ia_c    = particle.attribute_index (it, "creation_time");
  const int ia_l     = particle.attribute_index (it, "lifetime");
  const int ia_timeindex = particle.attribute_index (it, "timeindex");
  const int ia_class = particle.attribute_index (it, "class");
  const int ia_accrate = particle.attribute_index (it, "accretion_rate");
  const int ia_accrate_time = particle.attribute_index (it, "accretion_rate_time");
  const int ia_loc = particle.attribute_index (it, "is_local");

  // Attribrute stride lengths
  const int dm   = particle.stride(it, ia_m);
  const int dp   = particle.stride(it, ia_x);
  const int dv   = particle.stride(it, ia_vx);
  const int dl   = particle.stride(it, ia_l);
  const int dc   = particle.stride(it, ia_c);
  const int dmf  = particle.stride(it, ia_mf);
  const int dloc = particle.stride(it, ia_loc);
  
  /// Initialise pointers for particle attribute arrays
  enzo_float * pmass = 0;
  enzo_float * prevmass = 0;
  enzo_float * px   = 0;
  enzo_float * py   = 0;
  enzo_float * pz   = 0;
  enzo_float * pvx  = 0;
  enzo_float * pvy  = 0;
  enzo_float * pvz  = 0;
  enzo_float * pmetal = 0;
  enzo_float * pcreation  = 0;
  enzo_float * plifetime = 0;
  int        * ptimeindex = 0;
  int        * pclass     = 0;
  enzo_float * paccrate = 0;
  enzo_float * paccrate_time = 0;
  int64_t * is_local = 0;

  int ipp; // Particle index
  int ib; // Batch index

  int gx,gy,gz;
  field.ghost_depth (0, &gx, &gy, &gz);

  int nx, ny, nz;
  field.size ( &nx, &ny, &nz);

  int mx = nx + 2*gx;
  int my = ny + 2*gy;
  int mz = nz + 2*gz;

  // get pointers to field values
  enzo_float * density     = (enzo_float *) field.values("density");
  enzo_float * specific_internal_energy =
    (enzo_float *) field.values("internal_energy");
  enzo_float * potential   = (enzo_float *) field.values("potential");
  enzo_float * velocity_x  = (enzo_float *) field.values("velocity_x");
  enzo_float * velocity_y  = (enzo_float *) field.values("velocity_y");
  enzo_float * velocity_z  = (enzo_float *) field.values("velocity_z");
  enzo_float * metal = field.is_field("metal_density") ?
    (enzo_float *) field.values("metal_density") : NULL;

  int count = 0; // Number of particles formed
  // iterate over all cells (not including ghost zones)
  for (int iz=gz; iz<nz+gz; iz++){
    for (int iy=gy; iy<ny+gy; iy++){
      for (int ix=gx; ix<nx+gx; ix++){

        int i = INDEX(ix,iy,iz,mx,my);
        // need to compute this better for Grackle fields (on to-do list)

        // Apply the criteria for star formation

	/*(i) The first criteria is that the local gas density must
	      exceed the jeans density */
	
	double jeans_density;
	if(!check_jeans_density(specific_internal_energy[i],
				mean_cell_width,
				density[i],&jeans_density)) continue;

	// (ii) The second criteria is that the flow is converging
	if (!check_converging_flow(velocity_x, velocity_y,
				   velocity_z, i,
				   1, my, my*mz,dx,dy,dz)) continue;

	//(iii) Is cell the density maximum within the control volume
	if (!check_density_maximum(enzo_block,ix,iy,iz)) continue;

        

        // now create a star particle
	count++; 
        int my_particle = particle.insert_particles(it, 1);
        particle.index(my_particle, &ib, &ipp);
        pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
	prevmass = (enzo_float *) particle.attribute_array(it, ia_pm, ib);
	// particle mass is mass here
        pmass[ipp * dm] = (density[i] - jeans_density) * cell_volume;

	// set particle position to be at centre of cell, plus some
	// randomness
        px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
        py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
        pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);

	// generate random numbers from uniform distribution between 0 and 1
	double rnum1 = (double(rand())) / (double(RAND_MAX));
	double rnum2 = (double(rand())) / (double(RAND_MAX));
	double rnum3 = (double(rand())) / (double(RAND_MAX));
        px[ipp * dp] = lx + (ix - gx + 0.5) * dx + 0.1 * (rnum1-0.5) * dx;
        py[ipp * dp] = ly + (iy - gy + 0.5) * dy + 0.1 * (rnum2-0.5) * dy;
        pz[ipp * dp] = lz + (iz - gz + 0.5) * dz + 0.1 * (rnum3-0.5) * dz;

	// set particle velocity to be cell velocity
        pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
        pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
        pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib);
        pvx[ipp * dv] = velocity_x[i];
        pvy[ipp * dv] = velocity_y[i];
        pvz[ipp * dv] = velocity_z[i];

        // finalize attributes
        plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib);
        pcreation     = (enzo_float *) particle.attribute_array(it, ia_c, ib);
        pcreation[ipp * dc]     =  enzo_block->time();   // formation time
	// lifetime is hard-coded here, maybe should be a parameter
        plifetime[ipp * dl] =  10.0 * cello::Myr_s / enzo_units->time() ;
	//paccrate = (enzo_float *) particle.attribute_array(it, ia_accrate, ib);
        //paccrate_time = (enzo_float *) particle.attribute_array(it, ia_accrate_time,
	//	ib);
      //ptimeindex = (int *) particle.attribute_array(it, ia_timeindex, ib);
      //ptimeindex[io] = 0;
      //pclass = (int *) particle.attribute_array(it, ia_class, ib);
      //pclass[io] = SMS;  //see star_type enum in _enzo.hpp
      //*(&paccrate[io]) = 0.0; // I don't count formation as accretion per se
       //*(&paccrate_time[io]) = pform[io];

	/* Starting point of accretion rate pointer for this particle */
	//paccrate = &paccrate[io];
	//paccrate_time = &paccrate_time[io];
	//paccrate += sizeof(enzo_float);
	//paccrate_time += sizeof(enzo_float);

	/* Now initialise the array from this offset point onwards */
	/* STEFAN: I don't think this is necessary. I think we can set attribrutes
	   to be arrays */
	//for(int k = 1; k < SS_NTIMES; k++, paccrate++, paccrate_time++) {
      // *paccrate = k;
      //  *paccrate_time = k;
      //}
	
        if (metal){
          pmetal     = (enzo_float *) particle.attribute_array(it, ia_mf, ib);
          pmetal[ipp * dmf] = metal[i] / density[i];
        }

	/* Specify that particle is local to the block */
	is_local = (int64_t *) particle.attribute_array(it,ia_loc,ib);
	is_local[ipp * dloc] = 1;
	
        // Remove mass from grid and rescale fraction fields
	density[i] = jeans_density;
        if (density[i] < 0){
          CkPrintf("Smartstar: density index mass: %g %i %g\n",
                   density[i],i,density[i]*cell_volume);
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
   Defaults to parent class timestep if nothing declared here
double EnzoMethodStarMakerSmartStar::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}
*/

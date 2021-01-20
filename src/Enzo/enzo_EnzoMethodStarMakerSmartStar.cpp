/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodStarMakerSmartStar.cpp
/// @author     John Regan (john.regan@mu.ie)
/// @date
/// @brief      Functionality to allow for SMS, PopIII, PopII and Black Hole Formation
///
///     Derived star maker class that actually makes stars. This is partially
///     adapted after the star_maker_ssn method from Enzo and from the active
///     particle type SmartStars from Enzo

#include "cello.hpp"
#include "enzo.hpp"
#include "FofLib.hpp"
#include <time.h>

int FofList(int, enzo_float *, enzo_float, int *, int **, int ***);
// #define DEBUG_SF
#define SCALAR 1

//-------------------------------------------------------------------

EnzoMethodStarMakerSmartStar::EnzoMethodStarMakerSmartStar
()
  : EnzoMethodStarMaker()
{
  // To Do: Make the seed an input parameter
  srand(time(NULL)); // need randum number generator for later
  return;
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

  int count = 0;
  //CkPrintf("OK we at least check for smart star formation!\n");
  
  // Are we at the highest level?
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
  block->lower(&lx,&ly,&lz);
  accretion_radius_cells_ = enzo_config->method_star_maker_accretion_radius_cells;
  CkPrintf("Accretion radius cells: %d \n",accretion_radius_cells_);
  // declare particle position arrays
  //  default particle type is "star", but this will default
  //  to subclass particle_type
  const int it   = particle.type_index (this->particle_type());
  CkPrintf("Particle type index = %d \n",it);
  const int ia_m = particle.attribute_index (it, "mass");
  CkPrintf("ia_m = %d \n",ia_m);
  const int ia_pm = particle.attribute_index (it, "prevmass");
  CkPrintf("ia_pm = %d \n",ia_pm);
  const int ia_x = particle.attribute_index (it, "x");
  CkPrintf("ia_x = %d \n",ia_x);
  const int ia_y = particle.attribute_index (it, "y");
  CkPrintf("ia_y = %d \n",ia_y);
  const int ia_z = particle.attribute_index (it, "z");
  CkPrintf("ia_z = %d \n",ia_z);
  const int ia_vx = particle.attribute_index (it, "vx");
  CkPrintf("ia_vx = %d \n",ia_vx);
  const int ia_vy = particle.attribute_index (it, "vy");
  CkPrintf("ia_vy = %d \n",ia_vy);
  const int ia_vz = particle.attribute_index (it, "vz");
  CkPrintf("ia_vz = %d \n",ia_vz);
  const int ia_ax = particle.attribute_index(it,"ax");
  CkPrintf("ia_ax = %d \n",ia_ax);
  // additional particle attributes
  const int ia_metal = particle.attribute_index (it, "metal_fraction");
  CkPrintf("ia_metal = %d \n",ia_metal);
  const int ia_to    = particle.attribute_index (it, "creation_time");
  CkPrintf("ia_to = %d \n",ia_to);
  const int ia_l     = particle.attribute_index (it, "lifetime");
  CkPrintf("ia_l = %d \n",ia_l);
  const int ia_timeindex = particle.attribute_index (it, "timeindex");
  CkPrintf("ia_timeindex = %d \n",ia_timeindex);
  const int ia_class = particle.attribute_index (it, "class");
  CkPrintf("ia_class = %d \n",ia_class);
  const int ia_accrate = particle.attribute_index (it, "accretion_rate");
  CkPrintf("ia_accrate = %d \n",ia_accrate);
  const int ia_accrate_time = particle.attribute_index (it, "accretion_rate_time");
  CkPrintf("ia_accrate_time = %d \n",ia_accrate_time);
  const int ia_foo = particle.attribute_index (it, "foo");
  CkPrintf("ia_foo = %d \n",ia_foo);
  int ib  = 0; // batch counter
  int ipp = 0; // counter
  
  /// pointers for particle attribute arrays (later)
  enzo_float * pmass = 0;
  enzo_float * prevmass = 0;
  enzo_float * px   = 0;
  enzo_float * py   = 0;
  enzo_float * pz   = 0;
  enzo_float * pvx  = 0;
  enzo_float * pvy  = 0;
  enzo_float * pvz  = 0;
   ///
  enzo_float * pmetal = 0;
  enzo_float * pform  = 0;
  enzo_float * plifetime = 0;
  int        * ptimeindex = 0;
  int        * pclass     = 0;
  enzo_float * paccrate = 0;
  enzo_float * paccrate_time = 0;

  // obtain the particle stride length
  const int ps = particle.stride(it, ia_m);

  int rank = cello::rank();

  int gx,gy,gz;
  field.ghost_depth (0, &gx, &gy, &gz);

  int nx, ny, nz;
  field.size ( &nx, &ny, &nz);

  int mx = nx + 2*gx;
  int my = ny + 2*gy;
  int mz = nz + 2*gz;


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

  const double Zsolar = 0.02;  // TODO: Update to more accurate value

  // Idea for multi-metal species - group these using 'group'
  // class in IC parameter file and in SF / Feedback routines simply
  // check if this group exists, and if it does, loop over all of these
  // fields to assign particle chemical tags and deposit yields

  // compute the temperature (we need it here)
  // NB: @JR This doesn't give the correct temperature. I need to come
  // and look at this. 
  EnzoComputeTemperature compute_temperature
    (enzo_config->ppm_density_floor,
     enzo_config->ppm_temperature_floor,
     enzo_config->ppm_mol_weight,
     enzo_config->physics_cosmology);

  compute_temperature.compute(enzo_block);
  // iterate over all cells (not including ghost zones)
  //
  //   To Do: Allow for multi-zone star formation by adding mass in
  //          surrounding cells if needed to accumulte enough mass
  //          to hit target star particle mass ()
  for (int iz=gz; iz<nz+gz; iz++){
    for (int iy=gy; iy<ny+gy; iy++){
      for (int ix=gx; ix<nx+gx; ix++){

        int i = ix + mx*(iy + my*iz);

        // need to compute this better for Grackle fields (on to-do list)
        double mean_particle_mass = enzo_config->ppm_mol_weight * cello::mass_hydrogen;
	//        double ndens = rho_cgs / mean_particle_mass;
	double rho_cgs = density[i] * enzo_units->density();
        double mass_in_solar_masses  = density[i] *dx*dy*dz * enzo_units->mass() / cello::mass_solar;
        double metallicity = (metal) ? metal[i]/density[i]/Zsolar : 0.0;
	double jeans_density;
        //
        // Apply the criteria for star formation
        //

	//(i) The first criteria is that the local gas density must
	//    exceed the jeans density
	if(! this->check_jeans_density(temperature[i],
				       dx*enzo_units->length(),
				       density[i],&jeans_density)) continue;
	CkPrintf("First criteria passed!\n");
	
	//CkExit(-99);
	// (ii) The second criteria is that the velocity divergence is negative
	if (! this->check_velocity_divergence(velocity_x, velocity_y,
					      velocity_z, i,
                                              1, my, my*mz)) continue;

	CkPrintf("Second criteria passed!\n");
	//CkExit(-99);
	//(iii) Is cell the gravitational minimum over a jeans length
	double cellpos[3];
	cellpos[0] = lx + (ix - gx + 0.5) * dx;
	cellpos[1] = ly + (iy - gy + 0.5) * dy;
	cellpos[2] = lz + (iz - gz + 0.5) * dz;
	if(! this->check_gravitational_minimum(enzo_block, cellpos, potential[i],
					       temperature[i],
					       rho_cgs, enzo_units->length())) continue;
	CkPrintf("Final criteria passed!\n");
	//CkExit(-99);
	// (iv) 
	
#ifdef DONOTCOMPILE
	    
	  //if (! this->check_number_density_threshold(ndens)) continue;
        //if (! this->check_self_gravitating( mean_particle_mass, rho_cgs, temperature[i],
        //                                    velocity_x, velocity_y, velocity_z,
        //                                    enzo_units->length(), enzo_units->density(),
        //                                    i, 1, my, my*mz, dx, dy, dz)) continue;

        // AJE: TO DO ---
        //      If Grackle is used, check for this and use the H2
        //      fraction from there instead if h2 is used. Maybe could
        //      do this in the self shielding factor function

        // Only allow star formation out of the H2 shielding component (if used)
        const double f_h2 = this->h2_self_shielding_factor(density,
                                                           metallicity,
                                                           enzo_units->density(),
                                                           enzo_units->length(),
                                                           i, 1, my, my*mz,
                                                           dx, dy, dz);
        mass *= f_h2; // apply correction (f_h2 = 1 if not used)

       
        // Check whether mass in [min_mass, max_range] range and if specified, Jeans unstable
        //if (! this->check_mass(mass)) continue;

        double tdyn = sqrt(3.0 * cello::pi / 32.0 / cello::grav_constant /
                      (density[i] * enzo_units->density()));

        //
        // compute fraction that can / will be converted to stars this step
        // (just set to efficiency if dynamical time is ignored)
        //
        double star_fraction =  this->use_dynamical_time_ ?
                                std::min(this->efficiency_ * enzo_block->dt * enzo_units->time() / tdyn, 1.0) :
                                         this->efficiency_ ;
        // if this is less than the mass of a single particle,
        // use a random number draw to generate the particle
	CkPrintf("star_fraction = %g \n",star_fraction);
	CkPrintf("star_particle_min_mass_ = %g \n",this->star_particle_min_mass_);
	CkPrintf("star_particle_max_mass_ = %g \n",this->star_particle_max_mass_);
	CkPrintf("star_fraction*mass = %g \n",star_fraction*mass);
        if ( star_fraction * mass < this->star_particle_min_mass_){
          // get a random number
          double rnum = (double(rand())) / (double(RAND_MAX));
          double probability = this->efficiency_ * mass / this->star_particle_min_mass_;
	  CkPrintf("rnum = %g, probability = %g \n",rnum,probability);
          if (rnum > probability){
	      CkPrintf("Not forming stars");
              continue; // do not form stars
          } else{
            star_fraction = this->star_particle_min_mass_ / mass;
          }
        } else {
          // else allow the total mass of stars to form to be up to the
          // maximum particle mass OR the maximum gas->stars conversion fraction.
          // AJE: Note, this forces there to be at most one particle formed per
          //      cell per timestep. In principle, this could be bad if
          //      the computed gas->stars mass (above) is >> than maximum particle
          //      mass b/c it would artificially extend the lifetime of the SF
          //      region and presumably increase the amount of SF and burstiness
          //      of the SF and feedback cycle. Check this!!!!

          if (star_fraction * mass > this->star_particle_max_mass_){
#ifdef DEBUG_SF
            CkPrintf( "DEBUG_SF: SmartStar - SF mass = %g ; max particle mass = %g\n",
                                         star_fraction*mass, this->star_particle_max_mass_);
#endif
            star_fraction = this->star_particle_max_mass_ / mass;
          }

          star_fraction = std::min(star_fraction, this->maximum_star_fraction_);
        }
#endif
        count++; //

        // now create a star particle
        //    insert_particles( particle_type, number_of_particles )
        int my_particle = particle.insert_particles(it, 1);
	CkPrintf("My particle = %d. \t it = %d \n",my_particle,it);

        // For the inserted particle, obtain the batch number (ib)
        //  and the particle index (ipp)
        particle.index(my_particle, &ib, &ipp);
	CkPrintf("ib = %d, ipp = %d\n",ib,ipp);

        int io = ipp; // ipp*ps
	
        // pointer to mass array in block
        pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
	CkPrintf("ia_m = %d \n",ia_m);
	prevmass = (enzo_float *) particle.attribute_array(it, ia_pm, ib);
	CkPrintf("ia_pm = %d \n",ia_pm);
	CkPrintf("prevmass address = %p \n",prevmass);
	CkPrintf("prevmass[%d] = %g \n",io,prevmass[io]);
        pmass[io] = (density[i] - jeans_density) * dx * dy * dz;
	CkPrintf("io = %d \n",io);
	CkPrintf("jeans_density = %g \n",jeans_density);
	CkPrintf("Particle mass = %g \n",pmass[io]);
	//CkPrintf("prevmass address = %s \n",prevmass);
	//prevmass[io] = (density[i] - jeans_density) * dx * dy * dz;
	//CkPrintf("prevmass[io] = %g \n",prevmass[io]);
        px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
        py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
        pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);
	//CkPrintf("px address = %s \n",px);
	//CkPrintf("py address = %s \n",py);
	//CkPrintf("pz address = %s \n",pz);
        // need to double check that these are correctly handling ghost zones
        //   I believe lx is lower coordinates of active region, but
        //   ix is integer index of whole grid (active + ghost)
        //
        px[io] = lx + (ix - gx + 0.5) * dx;
        py[io] = ly + (iy - gy + 0.5) * dy;
        pz[io] = lz + (iz - gz + 0.5) * dz;
	

        pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
        pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
        pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib);

        pvx[io] = velocity_x[i];
        if (velocity_y) pvy[io] = velocity_y[i];
        if (velocity_z) pvz[io] = velocity_z[i];

        // finalize attributes
        plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib);
        pform     = (enzo_float *) particle.attribute_array(it, ia_to, ib);
	CkPrintf("ia_l = %d \n",ia_l);
	CkPrintf("pform = %d \n",ia_to);
        pform[io]     =  enzo_block->time();   // formation time
        plifetime[io] =  10.0 * cello::Myr_s / enzo_units->time() ; // lifetime

	paccrate = (enzo_float *) particle.attribute_array(it, ia_accrate, ib);
        paccrate_time = (enzo_float *) particle.attribute_array(it, ia_accrate_time, ib);
	CkPrintf("ia_accrate = %d \n",ia_accrate);
	CkPrintf("ia_accrate_time = %d \n",ia_accrate_time);
	ptimeindex = (int *) particle.attribute_array(it, ia_timeindex, ib);
	ptimeindex[io] = 0;
	pclass = (int *) particle.attribute_array(it, ia_class, ib);
	pclass[io] = SMS;  //see star_type enum in _enzo.hpp
	*(&paccrate[io]) = 0.0; // I don't count formation as accretion per se
	*(&paccrate_time[io]) = pform[io];
	
	CkPrintf("it = %d\t ia_accrate = %d\t ib = %d\n", it, ia_accrate, ib);
	CkPrintf("offset = 1: accrate + offset  = %p \n", &paccrate[io] + sizeof(enzo_float));
	CkPrintf("io = %d: accrate address = %p %p %p\n", io,
		 paccrate, &paccrate[io], &(paccrate[io]));
	CkPrintf("%d: accrate address = %p (*accrate = %f)\n", 0, &(paccrate[io]), *(&paccrate[io]));
	/* Starting point of accretion rate pointer for this particle */
	paccrate = &paccrate[io];
	paccrate_time = &paccrate_time[io];
	paccrate += sizeof(enzo_float);
	paccrate_time += sizeof(enzo_float);
	CkPrintf("%d: accrate address = %p (*accrate = %f)\n", 0, paccrate, *paccrate);
	/* Now initialise the array from this offset point onwards */
	for(int k = 1; k < SS_NTIMES; k++, paccrate++, paccrate_time++) {
	  *paccrate = k;
	  *paccrate_time = k;
	  CkPrintf("%d: accrate address = %p (*accrate = %f)\n", k, paccrate, *paccrate);
	 
	}
        if (metal){
          pmetal     = (enzo_float *) particle.attribute_array(it, ia_metal, ib);
          pmetal[io] = metal[i] / density[i];
        }

        // Remove mass from grid and rescale fraction fields
	density[i] = jeans_density;
        if (density[i] < 0){
          CkPrintf("Smartstar: density index mass: %g %i %g\n",
                   density[i],i,mass_in_solar_masses);
          ERROR("EnzoMethodStarMakerSmartStar::compute()",
                "Negative densities in star formation");
        }

          // rescale tracer fields to maintain constant mass fraction
        // with the corresponding new density...
        //enzo_float scale = density[i] / old_density;
        rescale_densities(enzo_block, i, 1.0 - jeans_density/density[i]);
      }
    }
  } // end loop iz

  if (count > 0){
      CkPrintf("SmartStar: Number of particles formed = %i \n", count);
      //CkExit(-1);
  }


  /* We should now merge particles which come within the mergin radius of each other */
  /* Particle merging can be disabled by setting SmartStarMerging = False */
  /* Find mergeable groups using an FOF search */

  // apply accretion depending on particle type
  // some particles do not accrete.
  // for now, just do this for all star particles 

  
  int numparticles = particle.num_particles(it);
  if(numparticles > 1 && count > 0) {
    CkPrintf("numparticles = %d\n", numparticles);
    int GroupNumberAssignment[numparticles];
    int *groupsize = NULL;
    int **grouplist = NULL;
    bool delete_mask[numparticles];
    enzo_float ParticleCoordinates[3*numparticles];
    /* Particles merge once they come within 3 accretion radii of one another */
    enzo_float MergingRadius = dx*accretion_radius_cells_*3; 
    px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
    py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
    pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);
    pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
    pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
    pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib);
    plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib);
    pform     = (enzo_float *) particle.attribute_array(it, ia_to, ib);
    pmetal    = (enzo_float *) particle.attribute_array(it, ia_metal, ib);
    
    
    for (int i=0; i<numparticles; i++) {
      particle.index(i, &ib, &ipp);
      int j = i*3;
      ParticleCoordinates[j]   = px[ipp];
      ParticleCoordinates[j+1] = py[ipp];
      ParticleCoordinates[j+2] = pz[ipp];
      CkPrintf("ipp = %d\n", ipp);
      delete_mask[i] = 0;
    }
    int ngroups = 0;
    
    /* 
     * Group particles into groups using a FoF algorithm
     * Need to merge particles in each group
     * GroupNumberAssignment is an array of size numparticles indexing each group to which 
     * particle belongs
     * groupsize is a pointer to an array with size equal to number of groups found. Each element 
     * gives number of particles in that group. groupsize allocated in routine. 
     * grouplist is a pointer to an array with size equal to number of groups found. 
     * Each element of the array is
     * itself an array listing the indices of all particles in that
     * group.
     */
    ngroups = FofList(numparticles, ParticleCoordinates, MergingRadius, 
		      GroupNumberAssignment, &groupsize, &grouplist);
    for (int i=0; i<numparticles; i++) {
      int j = i*3;
      CkPrintf("ParticleCoordinates = %e %e %e\n", ParticleCoordinates[j],
	       ParticleCoordinates[j+1],     ParticleCoordinates[j+2]);
      CkPrintf("MergingRadius = %e\n", MergingRadius);
    }
    
    for(int i = 0; i < ngroups; i++) {
      if (groupsize[i] != 1) {
	/* Particle 0 */
	int ippa = -1;
	particle.index(grouplist[i][0], &ib, &ippa);
	CkPrintf("P0 has mass %f msolar", pmass[ippa]/cello::mass_solar);
	for (int j=1; j < groupsize[i]; j++) {
	  int ippb = -1;
	  particle.index(grouplist[i][j], &ib, &ippb);
	  CkPrintf("P%d has mass %f msolar", j, pmass[ippb]/cello::mass_solar);
	  /* Merge everything into first particle in group */
	  enzo_float ratio1 = pmass[ippa] / (pmass[ippa] + pmass[ippb]);
	  enzo_float ratio2 = 1.0 - ratio1;

	  px[ippa] = ratio1 * px[ippa] + ratio2 * px[ippb];
	  py[ippa] = ratio1 * py[ippa] + ratio2 * py[ippb];
	  pz[ippa] = ratio1 * pz[ippa] + ratio2 * pz[ippb];
	  pvx[ippa] = ratio1 * pvx[ippa] + ratio2 * pvx[ippb];
	  pvy[ippa] = ratio1 * pvy[ippa] + ratio2 * pvy[ippb];
	  pvz[ippa] = ratio1 * pvz[ippa] + ratio2 * pvz[ippb];
	  plifetime[ippa] = std::min(plifetime[ippa], plifetime[ippb]);
	  pform[ippa] = std::min(pform[ippa], pform[ippb]);
	  pmetal[ippa] = ratio1 * pmetal[ippa] + ratio2 * pmetal[ippb];
	  pmass[ippa] += pmass[ippb];
	  delete_mask[grouplist[i][j]] = 1;
	  
	  
	}
	CkPrintf("groupsize[%d] = %d.\t Grouplist[%d][0] = %d\t GroupList[%d][1] = %d\n",
		 i, groupsize[i], i, grouplist[i][0], i, grouplist[i][1]);

      }
	
    }
    /* Delete redundant particles */
    int numdeleted = particle.delete_particles(it, ib, delete_mask);
    CkPrintf("FoF: ngroups = %d\n", ngroups);
    CkPrintf("Number of Particles after merging = %d\n", particle.num_particles(it));
    CkPrintf("Number of Particles Deleted = %d\n", numdeleted);
  }
  block->compute_done();

  return;
}

  
/*
   Defaults to parent class timestep if nothing declared here
double EnzoMethodStarMakerSmartStar::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}
*/

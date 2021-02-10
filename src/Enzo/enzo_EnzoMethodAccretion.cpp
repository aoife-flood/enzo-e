/// See LICENSE_CELLO file for license and copyright information
/// @file	enzo_EnzoMethodAccretion.cpp

///
///
///
///

#include "cello.hpp"
#include "enzo.hpp"
//#define UPDATE_STATS 1
EnzoMethodAccretion::EnzoMethodAccretion
()
  : Method()
{
  FieldDescr * field_descr = cello::field_descr();
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();


  cello::simulation()->new_refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();
  CkPrintf("Printing EnzoMethodAccretion's refresh object: \n");
  refresh->print();

  dual_energy_         = field_descr->is_field("internal_energy") &&
                         field_descr->is_field("total_energy");
  prescription_   = enzo_config->method_accretion_prescription;
  //Creating global gamma here for testing. 
  gamma_ = 5.0/3.0;
  return;
}

void EnzoMethodAccretion::pup (PUP::er &p)
{
  /// NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | prescription_;
  p | dual_energy_;

  return;
}

void EnzoMethodAccretion::compute (Block * block) throw()
{

  if (block->is_leaf()){
    this->compute_(block);
  }

  block->compute_done();

  return;
}

void EnzoMethodAccretion::compute_ (Block * block) throw()
{

  //CkPrintf("%s: Fantastic we will now attempt to accrete onto a star!\n", __FUNCTION__);
  //CkExit(-99);
  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();
  const int rank = cello::rank();
  Field field = block->data()->field();
  // Obtain grid sizes and ghost sizes

  enzo_float * d           = (enzo_float *) field.values ("density");
  enzo_float * te          = (enzo_float *) field.values("total_energy");
  enzo_float * ge          = (enzo_float *) field.values("internal_energy");
  enzo_float * temperature = (enzo_float *) field.values("temperature");
  enzo_float * metal = field.is_field("metal_density") ?
                       (enzo_float *) field.values("metal_density") : NULL;
  enzo_float * velocity_x = (rank >= 1) ?
    (enzo_float *)field.values("velocity_x") : NULL;
  enzo_float * velocity_y = (rank >= 2) ?
    (enzo_float *)field.values("velocity_y") : NULL;
  enzo_float * velocity_z = (rank >= 3) ?
    (enzo_float *)field.values("velocity_z") : NULL;
  int mx, my, mz, gx, gy, gz;
  double xm, ym, zm, xp, yp, zp, hx, hy, hz;
  int nx, ny, nz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);
  field.size ( &nx, &ny, &nz);
  block->data()->lower(&xm,&ym,&zm);
  block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);
  double dx, dy, dz;
  block->cell_width(&dx, &dy, &dz);
  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;
  // We will probably never be in the situation of constant acceleration
  // and cosmology, but just in case.....
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();
  enzo_float cosmo_a = 1.0;
  enzo_float accretion_rate = 0.0, accreted_mass = 0.0;
  
  double current_time  = block->time();
  if (cosmology) {
    enzo_float cosmo_dadt = 0.0;
    double dt    = block->dt();
    cosmology->compute_expansion_factor(&cosmo_a,&cosmo_dadt,current_time+0.5*dt);
    if (rank >= 1) hx *= cosmo_a;
    if (rank >= 2) hy *= cosmo_a;
    if (rank >= 3) hz *= cosmo_a;
  }

  double inv_vol = 1.0 / (hx*hy*hz), vol = hx*hy*hz;
  
  Particle particle = enzo_block->data()->particle();

  // apply accretion depending on particle type
  // some particles do not accrete.
  // for now, just do this for all star particles 

  int it = particle.type_index("star");
  int feedback_count = 0;
  CkPrintf("Number of Star Particles = %d \n",particle.num_particles(it));
  if (particle.num_particles(it) > 0 ){

    CkPrintf("Found a particle\n");


    
    const int ia_m = particle.attribute_index (it, "mass");
    const int ia_pm = particle.attribute_index (it, "prevmass");
    const int ia_x = (rank >= 1) ? particle.attribute_index (it, "x") : -1;
    const int ia_y = (rank >= 2) ? particle.attribute_index (it, "y") : -1;
    const int ia_z = (rank >= 3) ? particle.attribute_index (it, "z") : -1;

    const int ia_vx = (rank >= 1) ? particle.attribute_index (it, "vx") : -1;
    const int ia_vy = (rank >= 2) ? particle.attribute_index (it, "vy") : -1;
    const int ia_vz = (rank >= 3) ? particle.attribute_index (it, "vz") : -1;

    const int ia_l = particle.attribute_index (it, "lifetime");
    const int ia_c = particle.attribute_index (it, "creation_time");
    const int ia_timeindex = particle.attribute_index (it, "timeindex");
    const int ia_class = particle.attribute_index (it, "class");
    const int ia_acc = particle.attribute_index (it, "accretion_rate");
    const int ia_acc_time = particle.attribute_index (it, "accretion_rate_time");
    const int dm = particle.stride(it, ia_m);
    const int dp = particle.stride(it, ia_x);
    const int dl = particle.stride(it, ia_l);
    const int dc = particle.stride(it, ia_c);
    const int dacc = particle.stride(it, ia_acc);
    const int nb = particle.num_batches(it);
    //CkPrintf("Num Batches = %d\n", nb);

    for (int ib=0; ib<nb; ib++){
      enzo_float *px=0, *py=0, *pz=0;
      enzo_float *pvx=0, *pvy=0, *pvz=0;
      enzo_float *plifetime=0, *pcreation=0, *pmass=0, *paccrate=0, *paccrate_time=0;
      enzo_float *prevmass=0;
      int *ptimeindex=0, *pclass=0;
      prevmass = (enzo_float *) particle.attribute_array(it, ia_pm, ib);
      pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
      ptimeindex = (int *)  particle.attribute_array(it, ia_timeindex, ib);
      pclass = (int *)  particle.attribute_array(it, ia_class, ib);
      paccrate = (enzo_float *) particle.attribute_array(it, ia_acc, ib);
      paccrate_time = (enzo_float *) particle.attribute_array(it, ia_acc_time, ib);
      px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
      py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
      pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);

      pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
      pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
      pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib);

      plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib);
      pcreation = (enzo_float *) particle.attribute_array(it, ia_c, ib);
      const int paccradius = 4; // need to make this read in a parameter
      
      int np = particle.num_particles(it,ib);
      CkPrintf("Num Particles in this Batch= %d\n", np);
      
      /*Loop over all particles in this batch */
      for (int ip=0; ip<np; ip++){

        int ipdp = ip*dp;

	enzo_float pvel[3];
	pvel[0] = pvx[ipdp];
	pvel[1] = pvy[ipdp];
	pvel[2] = pvz[ipdp];
	// get corresponding grid position
        double xp = (px[ipdp] - xm) / hx;
        double yp = (py[ipdp] - ym) / hy;
        double zp = (pz[ipdp] - zm) / hz;

	
        // get 3D grid index for particle - account for ghost zones!!
        int ix = ((int) std::floor(xp))  + gx;
        int iy = ((int) std::floor(yp))  + gy;
        int iz = ((int) std::floor(zp))  + gz;

        // now deposit feedback in this cell
        int i = INDEX(ix,iy,iz,mx,my);
	enzo_float accretion_radius = paccradius*dx;
	/* Things to do:
	 * 1. Calculate the Bondi-Hoyle radius.
	 * 2. Compute Kernel Radius (see Krumholtz et al. (2004))
	 * 3. Compute weights (Krumholtz et al. (2004))
	 * 4. Compute Accretion rate 
	 */

	const enzo_float cell_temp = temperature[i];
	enzo_float cell_vel[3];
	cell_vel[0] = velocity_x[i];
	cell_vel[1] = velocity_y[i];
	cell_vel[2] = velocity_z[i];
	//CkPrintf("Pos = %e %e %e\n", px[ipdp], py[ipdp], pz[ipdp]);
	//CkPrintf("Vel = %e %e %e\n", pvel[0]/1e5, pvel[1]/1e5, pvel[2]/1e5);
	//CkPrintf("Cell temperature = %f\n", cell_temp);
	const enzo_float bondi_hoyle_radius =
	  calculate_bondi_hoyle_radius(pmass[ipdp], pvel, cell_temp, cell_vel); 
	//CkPrintf("Bondi-Hoyle Radius = %e cm\n", bondi_hoyle_radius);

	/* Compute kernel radius */
	enzo_float kernel_radius = 0.0;
	
	/* Compute the kernel radius */
	if (bondi_hoyle_radius < dx/4.0) {  /* For BHs whose Bondi radius is not resolved */
	  //printf("%s: Setting kernel radius to CellWidth, BH not resolved\n", __FUNCTION__);
	  kernel_radius = dx*2.0;
	}
	else { /*Accrete out to the BH radius */
	  //printf("%s: Setting kernel radius to BondiHoyleRadius\n", __FUNCTION__);
	  kernel_radius = std::max(bondi_hoyle_radius, accretion_radius);
	}

	int numcells = 0, num_ghost_cells = 0;
	double weighted_sum = 0.0, avg_temp = 0.0, total_gas_mass = 0.0, sum_of_weights = 0.0;
	/* Calculate cell weights within the accretion zone */
	/* Also identify if ghost zones fall within the accretion zone. */
	for (int iz = 0; iz < ngz; iz++){

	  for (int iy = 0; iy < ngy; iy++){
	    for (int ix = 0; ix < ngx; ix++){
	      int cellindex = INDEX(ix,iy,iz,ngx,ngy);
	      /* Get the positions of the cells */
	      double posx = xm + (ix - gx + 0.5) * dx;
	      double posy = ym + (iy - gy + 0.5) * dy;
	      double posz = zm + (iz - gz + 0.5) * dz;
	      
	      /* Calculate distance from cells to the particle. I need 
	       * to do this to weight each cells contribution */
	      double dist2 = pow((px[ipdp] - posx), 2.0) +
		             pow((py[ipdp] - posy), 2.0) +
		             pow((pz[ipdp] - posz), 2.0);
	      if ((accretion_radius*accretion_radius) > dist2) {
		weighted_sum += d[cellindex]*
		  exp(-dist2/((kernel_radius)*(kernel_radius)));
		sum_of_weights += exp(-dist2/((kernel_radius)*(kernel_radius)));
		avg_temp += temperature[cellindex];
		total_gas_mass += d[cellindex]*vol;
		numcells++;
		/*cells in ghost zone*/
		if(ix < gx || ix >(nx+gx) ||
		   iy < gy || iy >(ny+gy) ||
		   iz < gz || iz >(nz+gz))
		  {
		    //CkPrintf("So Cell %d %d %d is inside Accretion zone\n", ix, iy, iz);
		    num_ghost_cells++;
		  }
		
		
	      }
	    }
	  }
	}
	CkPrintf("numcells = %d\n", numcells);
	CkPrintf("total_gas_mass = %e Msun\n", total_gas_mass/cello::mass_solar);
	if(num_ghost_cells > 0) {
	  CkPrintf("num_ghost_cells = %d\n", num_ghost_cells);
	  //CkExit(-99);
	}
	double weight = 1.0/numcells;
	double avg_density = weighted_sum / sum_of_weights;
	const double mean_particle_mass = cello::mass_hydrogen*enzo_config->ppm_mol_weight;
	const double cInfinity = sqrt(gamma_ * cello::kboltz * cell_temp / mean_particle_mass);
	/* Estimate the relative velocity */
	enzo_float vInfinity = sqrt(pow(pvel[0] - cell_vel[0],2) +
				    pow(pvel[1] - cell_vel[1],2) +
				    pow(pvel[2] - cell_vel[2],2));
	/* 
	 * Traditional Bondi-Hoyle Prescription using the formalism
	 * presented in Krumholtz et al. (2004)
	 * The particle accretes according to the Bondi-Hoyle formula
	 * with a density given by the average density within the 
	 * accretion radius.
	 */

	if(prescription_ == SPHERICAL_BONDI_HOYLE_FORMALISM) {

	  CkPrintf("Doing SPHERICAL_BONDI_HOYLE_FORMALISM, Accretion Prescription = %d\n",
		   prescription_);
	  double rho_infinity = avg_density /
	    bondi_alpha((float)1.2*dx / bondi_hoyle_radius);
	  double lambda_c = 0.25*exp(1.5); // Only valid for isothermal gas
	  /* Bondi Hoyle */
	  //CkPrintf("rho_infinity = %e\n", rho_infinity);
	  //CkPrintf("cInfinity = %e\n", cInfinity);
	  //CkPrintf("vInfinity = %e\n", vInfinity);
	  
	  accretion_rate = (4*cello::pi*rho_infinity*pow(bondi_hoyle_radius,2)*
			   sqrt(pow(lambda_c*cInfinity,2) + pow(vInfinity,2)));
	  //CkPrintf("accretion_rate = %e\n", accretion_rate);
	}
	//CkPrintf("!!!!!!!!!!!!!accretion_rate = %e Msolar/yr\n", accretion_rate*cello::yr_s/cello::mass_solar);
	//CkExit(-99);

	enzo_float DeltaV[3] = {0.0, 0.0, 0.0}; 
	enzo_float ppos[3] = {0.0, 0.0, 0.0};
	ppos[0] = px[ipdp]; ppos[1] = py[ipdp]; ppos[2] = pz[ipdp];

	/* Remove mass accreted from the grid */
	if (!remove_accreted_mass(enzo_block, ppos, pvel, pmass[ipdp], kernel_radius,
				  sum_of_weights, accretion_radius,
				  accretion_rate, &accreted_mass,
				  DeltaV))
	  {
	    CkPrintf("%s: Failed to remove mass from grid after accretion. Cannot continue.\n",
		     __FUNCTION__);
	    CkExit(-99);
	  }

	/* Mass and energy of grid updated. */

	/*Now update particle */
	pvx[ipdp] += DeltaV[0];
	pvy[ipdp] += DeltaV[1];
	pvz[ipdp] += DeltaV[2];
	pmass[ipdp] += accreted_mass;
	//paccrate[ipdp] = accreted_mass/block->dt(); /* We calculate the accretion rate below */
        feedback_count++;

#ifdef UPDATE_STATS
	/* 
	 * Every X number of years we need to calculate the actual accretion rate 
	 * onto the star. This is necessary because the accretion rate 
	 * calculated using the accretion presciption (e.g. Bondi-Hoyle or Mass Flux)
	 * is a lower limit as it misses accretion due to mergers which can dominate
	 * at certain times. 
	 * To Calculate the actual accretion rate we simply look at the mass difference
	 * at fixed intervals and calculate accretion rate from that value. 
	 * The fixed interval could be user defined or perhaps determined by being
	 * a fixed number of timesteps of the max refinement level grid. 
	 */
	
	if(block->cycle() % 10 ==  0) {
	  CkPrintf("OK so we need to calculate the accretion rate here\n");
	  //enzo_float omass = prevmass[ipdp];
	  enzo_float cmass = pmass[ipdp];
	  if(cmass - omass < 0.0) {
	    CkPrintf("Weird - the previous mass is greater than the current mass....\n");
	    CkExit(-99);
	  }
	  int offset=ptimeindex[ipdp];
	  enzo_float otime = *(&(paccrate_time[ipdp])+sizeof(enzo_float)*offset);
	  enzo_float ctime = current_time;
	  //int ctimeindex = (ptimeindex[ipdp]++)%SS_NTIMES;
	  //int otimeindex = ctimeindex - 1;
	  //if(otimeindex == -1) //loop back
	  //  otimeindex = SS_NTIMES -1;

	  enzo_float deltatime = ctime - otime;
	  enzo_float accrate = (cmass - omass)/deltatime;
	  /* Which timeindex do we want to update */
	  

	  CkPrintf("Setting with accrate = %e\t accrate_time = %e\n", accrate, ctime);
	  /* Do the updates */
	  *(&paccrate[ipdp] + sizeof(enzo_float)*ptimeindex[ipdp]) = accrate;
	  *(&paccrate_time[ipdp] + sizeof(enzo_float)*ptimeindex[ipdp]) = ctime;
	  ptimeindex[ipdp]++;
	  for(int i = 0; i <= offset; i++) {
	    CkPrintf("%d: paccrate =  %p\t *paccrate = %e\n", i,
		     &paccrate[ipdp] + sizeof(enzo_float)*i,
		     (*(&paccrate[ipdp] + sizeof(enzo_float)*i))*cello::yr_s/cello::mass_solar);
	    CkPrintf("%d: accrate_time = %p\t *paccrate_time = %e\n", i, &paccrate_time[ipdp] +
		     sizeof(enzo_float)*i,
		     (*(&paccrate_time[ipdp] + sizeof(enzo_float)*i))/cello::kyr_s);
	  }
	 

	  //prevmass[ipdp] = cmass;

	  CkPrintf("old_mass = %e Msolar\t cmass = %e Msolar\t DeltaM = %f\n",
		   omass/cello::mass_solar,
		   cmass/cello::mass_solar,
		   (cmass - omass)/cello::mass_solar);
	  CkPrintf("old_time = %e kyr\t ctime = %e kyr\t DeltaT = %e kyr\n",
		   otime/cello::kyr_s,
		   ctime/cello::kyr_s,
		   (ctime - otime)/cello::kyr_s);
	  CkPrintf("accrate = %e Msolar/yr\t accratetime = %e kyr \t " \
		   "deltatime = %f kyr\t index = %d\t Particle Mass = %e Msolar\t Class = %d\n",
		   accrate*cello::yr_s/cello::mass_solar,
		   ctime/cello::kyr_s,
		   deltatime/cello::kyr_s,
		   offset,
		   pmass[ipdp]/cello::mass_solar,
		   pclass[ipdp]);
	  if(accrate*cello::yr_s/cello::mass_solar > CRITICAL_ACCRETION_RATE) {
	    pclass[ipdp] = SMS;
	  }
	  else {
	    float Age = block->time() - pcreation[ipdp];
	    if(Age/cello::yr_s > 1e4) { /* Don't do this at very start */
	      CkPrintf("%s: WARNING: ParticleClass switching from SMS to POPIII " \
		       " (deltatime = %f kyrs)\n", __FUNCTION__,
		       deltatime*(cello::yr_s*1e3));
	      CkPrintf("%s: WARNING: Accretion Rate = %f Msolar/yr. Critical rate = %f Msolar/yr\n",
		       __FUNCTION__,
		       (accrate*cello::yr_s/cello::mass_solar),
		       CRITICAL_ACCRETION_RATE);
	      pclass[ipdp] = popIII;
	    }
	    CkPrintf("Age = %e kyr\t Time = %e kyr", Age/cello::kyr_s, block->time()/cello::kyr_s);
	  }
	}
	
#endif

      } // end loop over particles

    } // end loop over batches


    

    
  }

  if (feedback_count > 0){
      CkPrintf("Number of feedback particles:  %i \n", feedback_count);
  }

 
  
  return;
}



int EnzoMethodAccretion::remove_accreted_mass (Block * block, enzo_float ppos[3],
					       enzo_float pvel[3], enzo_float pmass,
					       enzo_float kernel_radius,
					       enzo_float sum_of_weights,
					       enzo_float particle_accretion_radius,
					       enzo_float accretion_rate,
					       enzo_float *accreted_mass,
					       enzo_float *DeltaV)
{
  double SmallRho = 1e-20, SmallEint = 1e-20;
  double tiny_number = 1e-20;
  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();
  const double mean_particle_mass = cello::mass_hydrogen*enzo_config->ppm_mol_weight;
  Particle particle = enzo_block->data()->particle();
  Field field = enzo_block->data()->field();
  double xm, ym, zm, xp, yp, zp;
  double dx, dy, dz;
  block->cell_width(&dx, &dy, &dz);
  double dt = block->dt();
  int rank = cello::rank();
  int numcells = 0;
  int gx,gy,gz;
  field.ghost_depth (0, &gx, &gy, &gz);
  block->data()->lower(&xm,&ym,&zm);
  block->data()->upper(&xp,&yp,&zp);
  double cumulative_accreted_mass = 0.0, totalmass_before = 0.0, totalmass_after = 0.0;
  enzo_float AveragedVelocity[3] = {0.0, 0.0, 0.0};

  int nx, ny, nz;
  field.size ( &nx, &ny, &nz);
  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;
  enzo_float * density     = (enzo_float *) field.values("density");
  enzo_float * te          = (enzo_float *) field.values("total_energy");
  enzo_float * ge          = (enzo_float *) field.values("internal_energy");
  enzo_float * velocity_x = (rank >= 1) ?
    (enzo_float *)field.values("velocity_x") : NULL;
  enzo_float * velocity_y = (rank >= 2) ?
    (enzo_float *)field.values("velocity_y") : NULL;
  enzo_float * velocity_z = (rank >= 3) ?
    (enzo_float *)field.values("velocity_z") : NULL;
  enzo_float cellvolume = dx*dy*dz;
  for (int iz = 0; iz < ngz; iz++){
    for (int iy = 0; iy < ngy; iy++){
      for (int ix = 0; ix < ngx; ix++){
	int cellindex = INDEX(ix,iy,iz,ngx,ngy);
	enzo_float rho_cell = density[cellindex];
	enzo_float mcell = rho_cell*cellvolume;
	if(mcell < tiny_number || rho_cell < tiny_number) {
	  continue;
	}

	/* Get the positions of the cells */
	double cellposx = xm + (ix - gx + 0.5) * dx;
	double cellposy = ym + (iy - gy + 0.5) * dy;
	double cellposz = zm + (iz - gz + 0.5) * dz;
	/* Calculate distance from cells to the particle. I need 
	 * to do this to weight each cells contribution */
	double dist2 = pow((ppos[0] - cellposx), 2.0) +
		       pow((ppos[1] - cellposy), 2.0) +
		       pow((ppos[2] - cellposz), 2.0);

	double radius = sqrt(dist2);
	double vgas[3] = {0.0, 0.0, 0.0};
	vgas[0] = velocity_x[cellindex];
	if(rank > 1)
	  vgas[1] = velocity_y[cellindex];
	if(rank > 2)
	  vgas[2] = velocity_z[cellindex];

	double weight = exp(-dist2/(kernel_radius*kernel_radius))/sum_of_weights;
	if ( (particle_accretion_radius < radius) || (weight < tiny_number)) {
	  // outside the accretion radius
	  ;
	}
	else {  //Inside accretion radius
	  
	  // TE and GE are stored per unit mass - is this true in enzo-e - I think so (JR)....
	  // double etot = te[cellindex]*mcell;
	  // double eint = 0.0;
	  // if (dual_energy_)  /*Look into this and see where it gets set */
	  //   eint = ge[cellindex]*mcell;
	  // else
	  //   eint = etot - 0.5*mcell*
	  //     (vgas[0]*vgas[0] + vgas[1]*vgas[1] + vgas[2]*vgas[2]);
	  // double ke = 0.5*mcell*(vgas[0]*vgas[0] + vgas[1]*vgas[1] + vgas[2]*vgas[2]);
	
	  // Calculate mass we need to subtract from this cell
	  double maccreted =  dt * accretion_rate * weight;
	  if (maccreted > ACCRETION_LIMIT*mcell)
	    maccreted = ACCRETION_LIMIT*mcell;

	  double mnew = 0.0;
	  if ((mcell - maccreted)/cellvolume > SmallRho) {
	    mnew = mcell - maccreted;
	  }
	  else {
	    mnew = SmallRho*cellvolume;
	    maccreted = mcell - mnew;
	  }
	  mnew = mcell - maccreted;
	  maccreted = mcell - mnew;
	  double rhonew = mnew/cellvolume;
	  

	  numcells++;
	  // Compute new total internal energy. By construction,
	  // this keeps the specific internal energy constant after
	  // accretion
	  //double eintnew = eint * (1.0 - maccreted/mcell);
	  //
	  // Compute new total kinetic energy
	  //double kenew = ke * (1.0 - maccreted/mcell);
	  
	  // Compute the new total energy
	  //double etotnew = eintnew + kenew;
	  
	  // Update the densities
	  density[cellindex] -= maccreted/cellvolume;
	  
	  // Update the energies. ge is unchanged in this framework. 
	  //te[cellindex] = etotnew/mnew;
	  
	  // Check if mass or energy is too small, correct if necessary
	  if (density[cellindex] < SmallRho) {
	    density[cellindex] = SmallRho;
	    velocity_x[cellindex] = vgas[0];
	    velocity_y[cellindex] = vgas[1];
	    velocity_z[cellindex] = vgas[2];
	  }
	  // if (dual_energy_) {
	  //   if (ge[cellindex] < SmallEint) {
	  //     ge[cellindex] = SmallEint;
	  //   }
	  // }
	  // else {
	  //   if (te[cellindex] -
	  // 	0.5 * (pow(velocity_x[cellindex],2) +
	  // 	       pow(velocity_y[cellindex],2) +
	  // 	       pow(velocity_z[cellindex],2))
	  // 	< SmallEint) {
	  //     te[cellindex] = SmallEint +
	  // 	0.5 * (pow(velocity_x[cellindex],2) +
	  // 	       pow(velocity_y[cellindex],2) +
	  // 		   pow(velocity_z[cellindex],2));
	  //   }
	  // }
	  // Everything is OK we can update the particle
	  // Mass first
	  
	  *accreted_mass += maccreted;
	  cumulative_accreted_mass += maccreted;
	 
	  AveragedVelocity[0] += mcell*vgas[0];
	  AveragedVelocity[1] += mcell*vgas[1];
	  AveragedVelocity[2] += mcell*vgas[2];
	  totalmass_before += mcell;
	  totalmass_after += mnew;
#ifdef MAX_ACCRETION_RATE
	  /* Max Accretion Rate criteria */	  
	  if(*accreted_mass*cellvolume > MaxAccretionRate*this->dtFixed) {
	    printf("%s: We have removed the maximum allowed mass from the grid", __FUNCTION__);
	    printf("%s: Accreted Mass = %e Msolar\t Max Allowed = %e\n", __FUNCTION__,
		   *accreted_mass*cellvolume*MassUnits/SolarMass,  
		   MaxAccretionRate*this->dtFixed*MassUnits/SolarMass);
	    return ENZO_SUCCESS;
	  }
#endif
	}

      }
    }
  }
  if(numcells == 0) { //Nothing to do
    DeltaV[0] = 0.0, DeltaV[1] = 0.0; DeltaV[2] = 0.0;
    *accreted_mass= 0.0;
    return ENZO_SUCCESS;
  }

   /* Calculate mass weighted average velocity inside accretion sphere. */
  for(int i = 0; i < 3 ; i++)
    AveragedVelocity[i] /= totalmass_before;
  enzo_float NewVelocity[3] = {0.0, 0.0, 0.0};


 
  NewVelocity[0] = (pmass*pvel[0] + cumulative_accreted_mass*AveragedVelocity[0])/(pmass + cumulative_accreted_mass);
  NewVelocity[1] = (pmass*pvel[1] + cumulative_accreted_mass*AveragedVelocity[1])/(pmass + cumulative_accreted_mass);
  NewVelocity[2] = (pmass*pvel[2] + cumulative_accreted_mass*AveragedVelocity[2])/(pmass + cumulative_accreted_mass);
  for(int i = 0; i < 3; i++)
    DeltaV[i] = NewVelocity[i] - pvel[i];
  



  
  
  return ENZO_SUCCESS;
  
}

double EnzoMethodAccretion::timestep (Block * block) const throw()
{
  // In general this is not needed, but could imagine putting timestep
  // limiters in situations where, for example, one would want
  // dt < star_lifetime (or something like that), especially if
  // important things happen throughout the star's lifetime.

  return std::numeric_limits<double>::max();
}

double EnzoMethodAccretion::calculate_bondi_hoyle_radius(enzo_float pmass,
							 enzo_float *pvel,
							 enzo_float cell_temp,
							 enzo_float *cell_vel)
{
  double BHradius = 0.0;
  const EnzoConfig * enzo_config = enzo::config();
  /* Estimate the relative velocity */
  enzo_float vInfinity = sqrt(pow(pvel[0] - cell_vel[0],2) +
			      pow(pvel[1] - cell_vel[1],2) +
			      pow(pvel[2] - cell_vel[2],2));
  const double mean_particle_mass = cello::mass_hydrogen*enzo_config->ppm_mol_weight;
  const double cInfinity_squared = cello::kboltz * cell_temp / mean_particle_mass;
  BHradius = cello::grav_constant*pmass/(pow(vInfinity,2) + cInfinity_squared);
  //CkPrintf("%s: vInfinity = %f km/s\n", __FUNCTION__,  (vInfinity)/1e5);
  //CkPrintf("%s: cInfinity = %f km/s\n", __FUNCTION__,  sqrt(cInfinity_squared)/1e5);
  //CkPrintf("%s: CellTemperature = %f K\n", __FUNCTION__, cell_temp);
  //CkPrintf("%s: ParticleMass = %e Msolar\n", __FUNCTION__, pmass/cello::mass_solar);
  return BHradius;
}


/* Ported (copy and paste) direct from enzo-dev */
double EnzoMethodAccretion::bondi_alpha(float x) {

#define XMIN 0.01
#define XMAX 2.0
#define NTABLE 51

  float lambda_c, xtable, xtablep1, alpha_exp;
  int idx;

  /* This is a precomputed table of alpha values.  These correspond to x values
     that run from 0.01 to 2.0 with uniform logarithmic spacing.  The reason for
     this choice of range is that the asymptotic expressions are accurate to
     better than 2% outside this range */

  float alphatable[NTABLE] = {820.254, 701.882, 600.752, 514.341, 440.497, 377.381, 323.427,
			      277.295, 237.845, 204.1, 175.23, 150.524, 129.377, 111.27, 95.7613,
			      82.4745, 71.0869, 61.3237, 52.9498, 45.7644, 39.5963, 34.2989,
			      29.7471, 25.8338, 22.4676, 19.5705, 17.0755, 14.9254, 13.0714,
			      11.4717, 10.0903, 8.89675, 7.86467, 6.97159, 6.19825, 5.52812,
			      4.94699, 4.44279, 4.00497, 3.6246, 3.29395, 3.00637, 2.75612,
			      2.53827, 2.34854, 2.18322, 2.03912, 1.91344, 1.80378, 1.70804,
			      1.62439};

  // A constant that appears in the following formulae.  This hardcoded value is
  // valid for an isothermal gas.
  lambda_c = 0.25*exp(1.5);

  // deal with the off-the-table cases
  if (x < XMIN) 
    return lambda_c / sqrt(2.*x*x);
  else if (x >= XMAX)
    return exp(1./x);
  else {
    // we are on the table
    
    idx = floor((NTABLE-1)*log(x/XMIN)/log(XMAX/XMIN));
    xtable = exp(log(XMIN) + idx*log(XMAX/XMIN)/(NTABLE-1));
    xtablep1 = exp(log(XMIN) + (idx+1)*log(XMAX/XMIN)/(NTABLE-1));
    alpha_exp = log(x/xtable) / log(xtablep1/xtable);

    return alphatable[idx] * pow(alphatable[idx+1]/alphatable[idx],alpha_exp);
  }

#undef NTABLE
#undef XMIN
#undef XMAX

}


 

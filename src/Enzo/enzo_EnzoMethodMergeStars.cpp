/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodMergeStars.cpp
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date
/// @brief  Implements a star merging method class
///
///         This method is designed to be run after a method which forms stars,
///         and will not work properly otherwise. After copying particles from
///         neighbouring blocks during the refresh, this method then uses an
///         iterative procedure to merge particles together, so that no
///         two particles are within a 'merging radius' of each other. This
///         method also requires that blocks containing a star particle, and
///         all its neighbouring blocks, are on the maximum refinement level
///         (yet to be implemented).

#include "cello.hpp"
#include "enzo.hpp"
#include "FofLib.hpp"
#include <time.h>

#define DEBUG_MERGESTARS
int FofList(int, enzo_float *, enzo_float, int *, int **, int ***);

EnzoMethodMergeStars::EnzoMethodMergeStars()
  : Method()
{
  const EnzoConfig * enzo_config = enzo::config();
  ASSERT("EnzoMethodMergeStars::EnzoMethodMergeStars()",
	 "EnzoMethodMergeStars requires unigrid mode (Adapt : max_level = 0). "
	 "In future, we may put in a refinement condition that blocks containing "
	 "star particles or neighbouring such a block is at highest refinement "
	 "level", enzo_config->mesh_max_level == 0);

  // Merging radius is in units of the cell width. Must be at least twice the
  // accretion kernel radius (assumes accretion method is activated, may change this
  // in the future)
  merging_radius_cells_ = enzo_config->method_merge_stars_merging_radius_cells;

  ASSERT("EnzoMethodMergeStars::EnzoMethodMergeStars() ",
	 "merging_radius_cells_ must be at least twice the accretion radius "
	 "in order to ensure that accretion zones do not overlap",
         merging_radius_cells_ >=
	 2 * enzo_config->method_accretion_kernel_radius_cells);

  // Refresh copies all star particles from neighbouring blocks
  cello::simulation()->new_refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  ParticleDescr * particle_descr = cello::particle_descr();
  refresh->add_particle(particle_descr->type_index("star"),true);
  

}

void EnzoMethodMergeStars::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | merging_radius_cells_;
 
  return;
}


void EnzoMethodMergeStars::compute ( Block *block) throw()
{

  
  if (block->is_leaf()){
    this->compute_(block);
  }
  block->compute_done();

  return;
}

// Required
double EnzoMethodMergeStars::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}

void EnzoMethodMergeStars::compute_(Block * block)
{
  Hierarchy * hierarchy = cello::hierarchy();  
  EnzoBlock * enzo_block = enzo::block(block);
  Particle particle = enzo_block->data()->particle();
  int it = particle.type_index("star");
  int num_particles = particle.num_particles(it);
#ifdef DEBUG_MERGESTARS
  CkPrintf("In total, there are %d particles on Block %s \n",num_particles,
	   block->name().c_str());
#endif
  int ngroups = 1;
  //  if(num_particles > 1) {
    // Declare pointers to particle attributes
    enzo_float *px, *py, *pz, *pvx, *pvy, *pvz;
    enzo_float *plifetime, *pcreation, *pmass, *pmetal;
    int64_t *is_local;
    
    // Get attribute indices
    const int ia_m   = particle.attribute_index (it, "mass");
    const int ia_x   = particle.attribute_index (it, "x");
    const int ia_y   = particle.attribute_index (it, "y");
    const int ia_z   = particle.attribute_index (it, "z");
    const int ia_vx  = particle.attribute_index (it, "vx");
    const int ia_vy  = particle.attribute_index (it, "vy");
    const int ia_vz  = particle.attribute_index (it, "vz");
    const int ia_l   = particle.attribute_index (it, "lifetime");
    const int ia_c   = particle.attribute_index (it, "creation_time");
    const int ia_mf  = particle.attribute_index (it, "metal_fraction");
    const int ia_loc = particle.attribute_index (it, "is_local");
	
    // Attribrute stride lengths
    const int dm   = particle.stride(it, ia_m);
    const int dp   = particle.stride(it, ia_x);
    const int dv   = particle.stride(it, ia_vx);
    const int dl   = particle.stride(it, ia_l);
    const int dc   = particle.stride(it, ia_c);
    const int dmf  = particle.stride(it, ia_mf);
    const int dloc = particle.stride(it, ia_loc);
      
    // Array giving the FoF group number of each particle
    int * group_index = new int[num_particles];
      
    // group_sizes will be an array giving the number of particles in each
    // group
    int *  group_size;

    // group_lists will be an 'array of arrays'. Each element will be an
    // array containing the indices of particles belonging to a particular
    // group
    int ** group_list;
      
    // Array containing particle positions in 'block units'
    double * particle_coordinates_block_units = new double[3 * num_particles];
      
    // Merging radius in 'block units'
    double merging_radius_block_units;
      
    // Fill in particle coordinates array
    get_particle_coordinates_block_units_(enzo_block, it,
					  particle_coordinates_block_units,
					  &merging_radius_block_units);

    ngroups = FofList(num_particles, particle_coordinates_block_units,
			  merging_radius_block_units, 
			  group_index, &group_size, &group_list);

#ifdef DEBUG_MERGESTARS
    CkPrintf("The %d particles on Block %s are in %d FoF groups \n",num_particles,
	   block->name().c_str(),ngroups);
#endif 
      
    for (int i = 0; i < ngroups; i++){
#ifdef DEBUG_MERGESTARS
      CkPrintf("Group %d out of %d on block %s: Group size = %d \n",i, ngroups,
	       block->name().c_str(),group_size[i]);
#endif 
      if (group_size[i] > 1){
	  
	ASSERT("EnzoMethodMergeStars::compute_()",
	       "There is a FoF group containing a pair of star particles "
	       "on non-neighbouring blocks. Since this cannot be properly "
	       "dealt with we exit the programme here. This has likely "
	       "happened because the merging radius is too large in "
	       "comparison to the block size",
	       particles_in_neighbouring_blocks_(particle_coordinates_block_units,
						   group_list,group_size,i));
	// ib1 and ip1 correspond to the first particle in group i
	// ib2 and ip2 will be used to index the other particles in the
	// group
	int ib1, ib2, ip1, ip2;

	particle.index(group_list[i][0],&ib1,&ip1);

	// We get the attribrutes of this particle and store them in a new
	// variable
	pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib1);
	px = (enzo_float *) particle.attribute_array(it, ia_x, ib1);
	py = (enzo_float *) particle.attribute_array(it, ia_y, ib1);
	pz = (enzo_float *) particle.attribute_array(it, ia_z, ib1);
	pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib1);
	pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib1);
	pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib1);
	plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib1);
	pcreation = (enzo_float *) particle.attribute_array(it, ia_c, ib1);
	pmetal    = (enzo_float *) particle.attribute_array(it, ia_mf, ib1);

	enzo_float pmass1 = pmass[ip1*dm];
	enzo_float px1 = px[ip1*dp];
	enzo_float py1 = py[ip1*dp];
	enzo_float pz1 = pz[ip1*dp];
	double pos1[3] = {px1,py1,pz1};
	enzo_float pvx1 = pvx[ip1*dv];
	enzo_float pvy1 = pvy[ip1*dv];
	enzo_float pvz1 = pvz[ip1*dv];
	enzo_float plifetime1 = plifetime[ip1*dl];
	enzo_float pcreation1 = pcreation[ip1*dc];
	enzo_float pmetal1 = pmetal[ip1*dmf];

#ifdef DEBUG_MERGESTARS
	CkPrintf("0th particle in group %d out of %d on Block %s"
                 "Particle block index = %d, batch_index = %d, "
		 "particle batch index = %d. Mass = %g. Position = (%g,%g,%g)\n",
		 i,ngroups, block->name().c_str(),group_list[i][0],ib1,ip1,
		 pmass1,px1,py1,pz1);
#endif 
	  
	// now loop over the rest of the particles in this group, and merge
	// them in to the first particle

	for (int j = 1; j < group_size[i]; j++){

	  particle.index(group_list[i][j],&ib2,&ip2);

	  // get attributes of this particle
	  pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib2);
	  px = (enzo_float *) particle.attribute_array(it, ia_x, ib2);
	  py = (enzo_float *) particle.attribute_array(it, ia_y, ib2);
	  pz = (enzo_float *) particle.attribute_array(it, ia_z, ib2);
	  pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib2);
	  pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib2);
	  pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib2);
	  plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib2);
	  pcreation = (enzo_float *) particle.attribute_array(it, ia_c, ib2);
	  pmetal    = (enzo_float *) particle.attribute_array(it, ia_mf, ib2);
	  is_local = (int64_t *) particle.attribute_array(it, ia_loc, ib2);

	  // Mark this particle as a non-local particle so that it gets
	  // deleted later... a bit of a hack
	  is_local[ip2 * dloc] = 0;

	  enzo_float pmass2 = pmass[ip2*dm];
	  enzo_float px2 = px[ip2*dp];
	  enzo_float py2 = py[ip2*dp];
	  enzo_float pz2 = pz[ip2*dp];
	  enzo_float pvx2 = pvx[ip2*dv];
	  enzo_float pvy2 = pvy[ip2*dv];
	  enzo_float pvz2 = pvz[ip2*dv];
	  const double pos2[3] = {px2,py2,pz2};
	  enzo_float plifetime2 = plifetime[ip2*dl];
	  enzo_float pcreation2 = pcreation[ip2*dc];
	  enzo_float pmetal2 = pmetal[ip2*dmf];
	    
	  enzo_float f1 = pmass1 / (pmass1 + pmass2);
	  enzo_float f2 = 1.0 - f1;	  
	    
	  // Get the nearest periodic image of particle 2 to particle 1
	  double npi[3];
	  hierarchy->get_nearest_periodic_image(pos2,pos1,npi);

	  // Compute new properties of 'particle 1'
	  px1 = f1 * px1 + f2 * npi[0];
	  py1 = f1 * py1 + f2 * npi[1];
	  pz1 = f1 * pz1 + f2 * npi[2];
	  pvx1 = f1 * pvx1 + f2 * pvx2;
	  pvy1 = f1 * pvy1 + f2 * pvy2;
	  pvz1 = f1 * pvz1 + f2 * pvz2;
	  plifetime1 = std::min(plifetime1, plifetime2);
	  pcreation1 = std::min(pcreation1, pcreation2);
	  pmetal1 = f1 * pmetal1 + f2 * pmetal2;
	  pmass1 += pmass2;

#ifdef DEBUG_MERGESTARS
	  CkPrintf("%dth particle in group %d out of %d on Block %s "
		   "is merged into 0th particle. New properties of 0th particle: "
		   "Mass = %g. Position = (%g,%g,%g)\n",
		   j,i,ngroups, block->name().c_str(),pmass1,px1,py1,pz1);
#endif 
	  
	} // Loop over particles in group

	// Set new properties of 'particle 1'

	pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib1);
	px = (enzo_float *) particle.attribute_array(it, ia_x, ib1);
	py = (enzo_float *) particle.attribute_array(it, ia_y, ib1);
	pz = (enzo_float *) particle.attribute_array(it, ia_z, ib1);
	pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib1);
	pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib1);
	pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib1);
	plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib1);
	pcreation = (enzo_float *) particle.attribute_array(it, ia_c, ib1);
	pmetal    = (enzo_float *) particle.attribute_array(it, ia_mf, ib1);
	  
	pmass[ip1*dm] = pmass1;

	 
	double folded_pos[3];
	pos1[0] = px1;
	pos1[1] = py1;
	pos1[2] = pz1;
	// Fold position within the domain if periodic boundary positions
	hierarchy->get_folded_position(pos1,folded_pos);
	px[ip1*dp] = folded_pos[0];
	py[ip1*dp] = folded_pos[1];
	pz[ip1*dp] = folded_pos[2];

#ifdef DEBUG_MERGESTARS
	CkPrintf("After merging all particles in group %d out of %d on Block %s: "
	         "New properties of 0th particle: "
	         "Mass = %g. Position = (%g,%g,%g)\n",
	          i,ngroups, block->name().c_str(),pmass[ip1*dm],px[ip1*dp],
	          py[ip1*dp],pz[ip1*dp]);
#endif 
	
	pvx[ip1*dv] = pvx1;
	pvy[ip1*dv] = pvy1;
	pvz[ip1*dv] = pvz1;
	plifetime[ip1*dl] = plifetime1;
	pcreation[ip1*dc] = pcreation1;
	pmetal[ip1*dmf] = pmetal1;
	    
        }// if (group_size[i] > 1)


      // Mark particle as local /non-local if it is / is not in the block
#ifdef DEBUG_MERGESTARS
	CkPrintf("After merging all particles in group %d out of %d on Block %s: "
	    "particle_in_block = %d \n",i,ngroups,block->name().c_str(),
	    particle_in_block_(group_list[i][0],enzo_block, it));
#endif
	
	int ib, ip;
	particle.index(group_list[i][0],&ib,&ip);
	is_local = (int64_t *) particle.attribute_array(it, ia_loc, ib);
	is_local[ip*dloc] = particle_in_block_(group_list[i][0],enzo_block,it);
	
	}// Loop over Fof groups

      // Delete the dynamically allocated arrays

      delete [] group_index;
      free(group_size);
      free(group_list);
      delete [] particle_coordinates_block_units;

      //     }// if num_particles > 1
	
// Now we delete particles we marked as non-local.
  int delete_count = enzo_block->delete_particle_copies_(it);
  cello::simulation()->data_delete_particles(delete_count);

#ifdef DEBUG_MERGESTARS
  CkPrintf("After merging, ngroups = %d, num_particles = %d \n", ngroups,
	    particle.num_particles(it));
#endif
  return;
  
}

// This function loops over all particles on the block, and returns 1 if any
// of them are local
bool EnzoMethodMergeStars::any_local_particles_(EnzoBlock * enzo_block, int it)
{
  bool return_val = 0;
  Particle particle = enzo_block->data()->particle();
  const int ia_loc = particle.attribute_index (it, "is_local");
  const int dloc = particle.stride(it, ia_loc);
  const int nb = particle.num_batches(it);

  // Loop over batches
  for (int ib = 0; ib < nb; ib++){
    int64_t * is_local;
    is_local = (int64_t *) particle.attribute_array(it,ia_loc,ib);
    const int np = particle.num_particles(it,ib);

    // Loop over particles in batch
    for (int ip = 0; ip < np; ip++){
      CkPrintf("is_local[%d] = %d \n",ip,is_local[ip*dloc]);
      if (is_local[ip * dloc] == 1){
	return_val = 1;
	break; // breaks out of ip loop
      }
    }
    if (return_val == 1) break; // breaks out of ib loop
  }
  return return_val;
}

// This fills an array, which must have already been allocated with length
// 3 * num_particles, with particle coordinates in 'block units'. This means
// that the centre of the block has coordinates (0,0,0), and a block has unit
// side length. This also takes care of periodic boundary conditions.
void EnzoMethodMergeStars::get_particle_coordinates_block_units_
  (EnzoBlock * enzo_block, int it,
   double * particle_coordinates_block_units,
   double * merging_radius_block_units)
{
  Hierarchy * hierarchy = cello::hierarchy();
  
  // Get dimensions of block
  double block_xm, block_ym, block_zm, block_xp, block_yp, block_zp;
  enzo_block->lower(&block_xm,&block_ym,&block_zm);
  enzo_block->upper(&block_xp,&block_yp,&block_zp);

  // Get the cell width, must be the same in all dimensions
  double cell_width_x, cell_width_y, cell_width_z;
  enzo_block->cell_width(&cell_width_x,&cell_width_y,&cell_width_z);

  // We assume that blocks are always cubes, will put in a check somewhere else
  // to make sure this is always the case
  const double block_width = block_xp - block_xm;

  // Get the merging radius in block units
  *merging_radius_block_units = merging_radius_cells_ * cell_width_x / block_width;

  // Get the coordinates of the centre of the block
  const double block_centre_x = 0.5 * (block_xm + block_xp);
  const double block_centre_y = 0.5 * (block_ym + block_yp);
  const double block_centre_z = 0.5 * (block_zm + block_zp);
  const double block_centre[3] = {block_centre_x,block_centre_y,block_centre_z};
  
  Particle particle = enzo_block->data()->particle();
  const int ia_x = particle.attribute_index (it, "x");
  const int ia_y = particle.attribute_index (it, "y");
  const int ia_z = particle.attribute_index (it, "z");
  const int dp = particle.stride(it, ia_x);
  enzo_float *px, *py, *pz;

  // ip_block is particle index within the block
  // ib is batch index
  // ip_batch is particle index within a batch
  int ib, ip_batch;

  const int num_particles = particle.num_particles(it);

  // Loop over all particles in block
  for (int ip_block = 0; ip_block < num_particles; ip_block++){

    // Get the particle's batch index and its index within the batch
    particle.index(ip_block,&ib,&ip_batch);
    
    // Get pointers to the attribute arrays
    px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
    py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
    pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);

    // Get the nearest periodic image to the block centre. If
    // boundary conditions are non-periodic, this just returns the particle
    // coordinates
    double npi[3];
    const double pos[3] = {px[ip_batch*dp],py[ip_batch*dp],pz[ip_batch*dp]};
    hierarchy->get_nearest_periodic_image(pos,block_centre,npi);

    // Now we can set particle coordinates in block units

    particle_coordinates_block_units[3*ip_block]
                             = (npi[0] - block_centre[0]) / block_width;
    particle_coordinates_block_units[3*ip_block + 1]
                             = (npi[1] - block_centre[1]) / block_width; 
    particle_coordinates_block_units[3*ip_block + 2]
                             = (npi[2] - block_centre[2]) / block_width; 
    
  } // ip loop
  
  return;
}

// Checks if all the particles within group i are in neighbouring blocks
bool EnzoMethodMergeStars::particles_in_neighbouring_blocks_
(double * particle_coordinates,
 int ** group_list,int * group_size,
 int i)
{
  bool return_val = 1;
  // Loop over all pairs of particles
  for (int j = 0; j < group_size[i]; j++){
    const int ind_1 = group_list[i][j];
    const double px1 = particle_coordinates[3*ind_1];
    const double py1 = particle_coordinates[3*ind_1 + 1];
    const double pz1 = particle_coordinates[3*ind_1 + 2];
    for (int k = j; k < group_size[i]; k++){
      const int ind_2 = group_list[i][k];
      const double px2 = particle_coordinates[3*ind_2];
      const double py2 = particle_coordinates[3*ind_2 + 1];
      const double pz2 = particle_coordinates[3*ind_2 + 2];

      // in block units, just need to check if one of the pair is less than -0.5, and
      // the other greater than 0.5, for each coordinate
      if (px1 < -0.5 && px2 > 0.5){
	return_val = 0;
	break; // break out of the k loop
      }

      if (px2 < -0.5 && px1 > 0.5){
	return_val = 0;
	break; // break out of the k loop
      }

      if (py1 < -0.5 && py2 > 0.5){
	return_val = 0;
	break; // break out of the k loop
      }

      if (py2 < -0.5 && py1 > 0.5){
	return_val = 0;
	break; // break out of the k loop
      }

      if (pz1 < -0.5 && pz2 > 0.5){
	return_val = 0;
	break; // break out of the k loop
      }

      if (pz2 < -0.5 && pz1 > 0.5){
	return_val = 0;
	break; // break out of the k loop
      }  
    } // k loop

    if (return_val == 0) break; // break out of the j loop
  }
  
  return return_val;
}

 bool EnzoMethodMergeStars::particle_in_block_(int i,  EnzoBlock * enzo_block, int it)
{
  bool return_val = 1;

  double block_xm, block_ym, block_zm, block_xp, block_yp, block_zp;
  enzo_block->lower(&block_xm,&block_ym,&block_zm);
  enzo_block->upper(&block_xp,&block_yp,&block_zp);
  
  Particle particle = enzo_block->data()->particle();
  enzo_float *px, *py, *pz;
  const int ia_x   = particle.attribute_index (it, "x");
  const int ia_y   = particle.attribute_index (it, "y");
  const int ia_z   = particle.attribute_index (it, "z");
  const int dp     = particle.stride(it, ia_x);

  int ib, ip;
  particle.index(i,&ib,&ip);
  px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
  py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
  pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);

  // Check if particle is outside the bounds of the block
  if (px[ip * dp] < block_xm || px[ip * dp] > block_xp ||
      py[ip * dp] < block_ym || py[ip * dp] > block_yp ||
      pz[ip * dp] < block_zm || pz[ip * dp] > block_zp)
    return_val = 0;
  
  return return_val;
}

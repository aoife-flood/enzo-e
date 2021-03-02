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

int FofList(int, enzo_float *, enzo_float, int *, int **, int ***);

EnzoMethodMergeStars::EnzoMethodMergeStars()
  : Method()
{

  const EnzoConfig * enzo_config = enzo::config();

  // Refresh copies all star particles from neighbouring blocks
  cello::simulation()->new_refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  ParticleDescr * particle_descr = cello::particle_descr();
  refresh->add_particle(particle_descr->type_index("star"),true);
  
  // Merging radius is in units of the cell width. Must be less than or equal to
  // the field ghost depth and equal to the accretion radius if an accretion method
  // is used.
  merging_radius_            = enzo_config->method_merge_stars_merging_radius;
}

void EnzoMethodMergeStars::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | merging_radius_;
 
  return;
}


void EnzoMethodMergeStars::compute ( Block *block) throw()
{

  //  if (! block->is_leaf()) return;
  CkPrintf("%s : %s : %d: Message from EnzoMethodMergeStars \n" , __FILE__,__FUNCTION__,__LINE__); 
  block->compute_done();

  return;
}

// Required
double EnzoMethodMergeStars::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}

#ifdef DONOTCOMPILE
  /* We should now merge particles which come within the mergin radius of each other */
  /* Particle merging can be disabled by setting SmartStarMerging = False */
  /* Find mergeable groups using an FOF search */

  // apply accretion depending on particle type
  // some particles do not accrete.
  // for now, just do this for all star particles 

  
  int numparticles = particle.num_particles(it);
  if(numparticles > 1 && count > 0) {
    //CkPrintf("numparticles = %d\n", numparticles);
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
      //CkPrintf("ipp = %d\n", ipp);
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
    // for (int i=0; i<numparticles; i++) {
    //   int j = i*3;
    //   //CkPrintf("ParticleCoordinates = %e %e %e\n", ParticleCoordinates[j],
    // 	       ParticleCoordinates[j+1],     ParticleCoordinates[j+2]);
    //   //CkPrintf("MergingRadius = %e\n", MergingRadius);
    // }
    // //CkExit(-1);
    
    for(int i = 0; i < ngroups; i++) {
      if (groupsize[i] != 1) {
	
	/* Particle 0 */
	int ippa = -1;
	particle.index(grouplist[i][0], &ib, &ippa);
	//CkPrintf("ippa = %d \n",ippa);
	//CkPrintf("P0 has mass %f msolar\n", pmass[ippa]/cello::mass_solar);
	for (int j=1; j < groupsize[i]; j++) {
	  int ippb = -1;
	  particle.index(grouplist[i][j], &ib, &ippb);

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
    //CkExit(-1);
    /* Delete redundant particles */
    int numdeleted = particle.delete_particles(it, ib, delete_mask);
    CkPrintf("FoF: ngroups = %d\n", ngroups);
    CkPrintf("Number of Particles after merging = %d\n", particle.num_particles(it));
    CkPrintf("Number of Particles Deleted = %d\n", numdeleted);
    
    
    
    px = (enzo_float *) particle.attribute_array(it, ia_x, ib);
    py = (enzo_float *) particle.attribute_array(it, ia_y, ib);
    pz = (enzo_float *) particle.attribute_array(it, ia_z, ib);
    for (int i2=0; i2<particle.num_particles(it); i2++) {
      particle.index(i2, &ib, &ipp);
      int j2 = i2*3;
      ParticleCoordinates[j2]   = px[ipp];
      ParticleCoordinates[j2+1] = py[ipp];
      ParticleCoordinates[j2+2] = pz[ipp];
      CkPrintf("After merging, particle %d: x,y,z = %g,%g,%g\n", i2,px[ipp],py[ipp],pz[ipp]);
    }
      
  
  } 
  block->compute_done();

  return;
}
#endif //DONOTCOMPILE

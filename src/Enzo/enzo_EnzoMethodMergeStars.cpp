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

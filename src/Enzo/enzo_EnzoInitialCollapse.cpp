// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialCollapse.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:30:35 PST 2011
/// @brief    Implementation of an array of Collapse problems, one per Block
///
/// This problem is designed for scaling studies of the Gravity
/// solver.  Each block contains a spherical collapse problem.

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

// #define DEBUG_PERFORMANCE
#define UNIFORM_DENSITY_PROFILE 1
#define R2_PROFILE              2
//----------------------------------------------------------------------

void EnzoInitialCollapse::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  PUParray(p,array_,3);
  p | radius_relative_;
  
}

//----------------------------------------------------------------------
void EnzoInitialCollapse::enforce_block
( Block * block, const Hierarchy  * hierarchy ) throw()

{

  if (!block->is_leaf()) return;

  Timer timer;
  timer.start();

  ASSERT("EnzoInitialCollapse",
	 "Block does not exist",
	 block != NULL);

  Field field = block->data()->field();

  // Get Field parameters
  
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  //CkPrintf("%s: (nx,ny,nz) = (%d,%d,%d)\n", __FUNCTION__, nx,ny,nz);
  // domain extents
  double dxm,dym,dzm;
  double dxp,dyp,dzp;
  hierarchy->lower(&dxm,&dym,&dzm);
  hierarchy->upper(&dxp,&dyp,&dzp);
  //CkPrintf("%s: (dxm,dym,dzm) = (%g,%g,%g)\n", __FUNCTION__, dxm,dym,dzm);
  //CkPrintf("%s: (dxp,dyp,dzp) = (%g,%g,%g)\n", __FUNCTION__, dxp,dyp,dzp);
  // Block extents
  double bxm,bym,bzm;
  double bxp,byp,bzp;
  block->data()->lower(&bxm,&bym,&bzm);
  block->data()->upper(&bxp,&byp,&bzp);
  //CkPrintf("%s: (bxm,bym,bzm) = (%g,%g,%g)\n", __FUNCTION__, bxm,bym,bzm);
  //CkPrintf("%s: (bxp,byp,bzp) = (%g,%g,%g)\n", __FUNCTION__, bxp,byp,bzp);
  const int rank = cello::rank();
  //CkPrintf("%s: rank = %d\n", __FUNCTION__, rank);
  double hx,hy,hz;
  field.cell_width(bxm,bxp,&hx,
		   bym,byp,&hy,
		   bzm,bzp,&hz);
  if (rank < 2) hy = 0.0;
  if (rank < 3) hz = 0.0;

  //CkPrintf("%s: (hx,hy,hz) = (%g,%g,%g)\n", __FUNCTION__, hx,hy,hz);

  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);
  //CkPrintf("%s: (gx,gy,gz) = (%d,%d,%d)\n", __FUNCTION__, gx,gy,gz);
  const int mx = nx + 2*gx;
  const int my = ny + 2*gy;
  const int mz = nz + 2*gz;

  const int m = mx*my*mz;
  //CkPrintf("%s: (mx,my,mz) = (%d,%d,%d)\n", __FUNCTION__, mx,my,mz);
  //CkPrintf("%s: m = %d \n", __FUNCTION__, m);
  
  // Get Fields
  bool RotatingSphere = true; 
  enzo_float *  d = (enzo_float *) field.values ("density");
  enzo_float * dt = (enzo_float *) field.values ("density_total");
  enzo_float *  p = (enzo_float *) field.values ("pressure");
  enzo_float * po = (enzo_float *) field.values ("potential");
  enzo_float *  t = (enzo_float *) field.values ("temperature");
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

  // Initialize Fields

  const int in = cello::index_static();
  //CkPrintf("%s: in = %d\n", __FUNCTION__, in);
  const double gamma = EnzoBlock::Gamma[in];
  //CkPrintf("%s: gamma = %g \n", __FUNCTION__,gamma);
  double mu = 3.0; // To match Federrath 2010
  const double energy = cello::kboltz*temperature_ / ((gamma - 1.0) * (mu * cello::mass_hydrogen));
  //CkPrintf("%s: energy = %g\n", __FUNCTION__, energy);
  //CkPrintf("%s: temperature_ = %g\n", __FUNCTION__, temperature_);
  
  // ...compute ellipsoid density

  const double rx = (dxp - dxm)/2. * radius_relative_ / array_[0] ;
  const double ry = (dyp - dym)/2. * radius_relative_ / array_[1] ;
  const double rz = (dzp - dzm)/2. * radius_relative_ / array_[2] ;
  //CkPrintf("%s: array_[0] = %d\n", __FUNCTION__, array_[0]);
  //CkPrintf("%s: array_[1] = %d\n", __FUNCTION__, array_[1]);
  //CkPrintf("%s: array_[2] = %d\n", __FUNCTION__, array_[2]);
  //CkPrintf("%s: radius_relative = %g\n", __FUNCTION__, radius_relative_);
  //CkPrintf("%s: dxp-dxm = %g\n", __FUNCTION__, dxp-dxm);
  //CkPrintf("%s: dyp-dym = %g\n", __FUNCTION__, dyp-dym);
  //CkPrintf("%s: dzp-dzm = %g\n", __FUNCTION__, dzp-dzm);
  //CkPrintf("%s: (rx,ry,rz) = (%g,%g,%g)\n", __FUNCTION__, rx,ry,rz);
  const double rx2i = 1.0/(rx*rx); 
  const double ry2i = 1.0/(ry*ry);
  const double rz2i = 1.0/(rz*rz);
  double density = 0.0;
  // This is the density at the trucation radius
  if(density_profile_ == R2_PROFILE) {
    if(truncation_density_ != 0.0){
      //CkPrintf("Truncation Density set. Overriding mass.\n");
      density = truncation_density_;
    }
    else
      density = mass_/(4.0*(cello::pi)*rx*ry*rz);
  }
  else { //This is the mean density
    density = mass_ / (4.0/3.0*(cello::pi)*rx*ry*rz);
  }
  
  double soundspeed = sqrt(temperature_*(cello::kboltz)/((mu)*(cello::mass_hydrogen))); 
  //CkPrintf("%s: Profile = %d\n", __FUNCTION__, density_profile_);
  //CkPrintf("%s: Density = %e\n", __FUNCTION__, density);
  //CkPrintf("%s: mass = %e\n", __FUNCTION__, mass_);
  //CkPrintf("%s: rx = %e\n", __FUNCTION__, rx);
  //CkPrintf("%s: Sound Speed = %e km/s\n", __FUNCTION__, soundspeed/1e5);
  //CkPrintf("%s: CellWidth = %e cm\n", __FUNCTION__, hx);
  //if(density_profile_ == UNIFORM_DENSITY_PROFILE)
    //CkPrintf("%s: calculated mass (assuming uniform density) = %e\n",
  //	     __FUNCTION__, (density*(4.0/3.0*(cello::pi)*rx*ry*rz))/cello::mass_solar);
  //if(density_profile_ == R2_PROFILE)
    //CkPrintf("%s: calculated mass (assuming non-uniform density) = %g solar masses\n",
  //	     __FUNCTION__, (density*(4.0*(cello::pi)*rx*ry*rz))/cello::mass_solar);
  //CkPrintf("%s: Instability Parameter = %f\n", __FUNCTION__,
  //	   4*(cello::pi)*(cello::grav_constant)*density*rx*rx/(soundspeed*soundspeed));
  // bounds of possible explosions intersecting this Block
  //CkExit(-99);
  int kxm = MAX((int)floor((bxm-dxm-rx)/(dxp-dxm)*array_[0])-1,0);
  int kym = MAX((int)floor((bym-dym-ry)/(dyp-dym)*array_[1])-1,0);
  int kzm = MAX((int)floor((bzm-dzm-rz)/(dzp-dzm)*array_[2])-1,0);
  int kxp = MIN( (int)ceil((bxp-dxm+rx)/(dxp-dxm)*array_[0])+1,array_[0]);
  int kyp = MIN( (int)ceil((byp-dym+ry)/(dyp-dym)*array_[1])+1,array_[1]);
  int kzp = MIN( (int)ceil((bzp-dzm+rz)/(dzp-dzm)*array_[2])+1,array_[2]);
  //CkPrintf("%s: (kxm,kym,kzm) = (%d,%d,%d)\n", __FUNCTION__, kxm,kym,kzm);
  //CkPrintf("%s: (kxp,kyp,kzp) = (%d,%d,%d)\n", __FUNCTION__, kxp,kyp,kzp);
  double hxa = (dxp-dxm) / array_[0];
  double hya = (rank >= 2) ? (dyp-dym) / array_[1] : 0.0;
  double hza = (rank >= 3) ? (dzp-dzm) / array_[2] : 0.0;

  //CkPrintf("%s: (hxa,hya,hza) = (%g,%g,%g)\n", __FUNCTION__, hxa,hya,hza);
  // (kx,ky,kz) index bounds of collapse in domain

  // Initialize background 

  // ratio of density inside and outside the cloud 
  const double density_ratio = 228.33;
  
  std::fill_n(d,m,density / density_ratio);
  std::fill_n(te,m,energy);
  std::fill_n(ie,m,energy);
  std::fill_n(t,m,temperature_);
  std::fill_n(dt,m,0.0);
  std::fill_n(p,m,0.0);
  std::fill_n(po,m,0.0);
  std::fill_n(ax,m,0.0);
  std::fill_n(vx,m,0.0);
  if (rank >= 2) std::fill_n(ay,m,0.0);
  if (rank >= 2) std::fill_n(vy,m,0.0);
  if (rank >= 3) std::fill_n(az,m,0.0);
  if (rank >= 3) std::fill_n(vz,m,0.0);
  std::fill_n(x,m,0.0);
  std::fill_n(b,m,0.0);

  // Initialize sphere (ellipsoid)

  for (int kz=kzm; kz<kzp; kz++) {
    double zc = dzm + hza*(0.5+kz);
    for (int ky=kym; ky<kyp; ky++) {
      double yc = dym + hya*(0.5+ky);
      for (int kx=kxm; kx<kxp; kx++) {
	double xc = dxm + hxa*(0.5+kx);

	// (collapse center xc,yc,zc)
	//CkPrintf("%s: (kx,ky,kz) = (%d,%d,%d)\n", __FUNCTION__, kx,ky,kz);
	//CkPrintf("%s: (xc,yc,zc) = (%g,%g,%g)\n", __FUNCTION__, xc,yc,zc);
	for (int iz=0; iz<mz; iz++) {
	  double z = bzm + (iz - gz + 0.5)*hz - zc;
	  for (int iy=0; iy<my; iy++) {
	    double y = bym + (iy - gy + 0.5)*hy - yc;
	    for (int ix=0; ix<mx; ix++) {
	      double x = bxm + (ix - gx + 0.5)*hx - xc;
	      int i = INDEX(ix,iy,iz,mx,my);
	      double R2 = x*x + y*y + z*z;
              double r2 = x*x*rx2i + y*y*ry2i + z*z*rz2i;
	      bool in_sphere = (r2 <= 1.0);
	      if (in_sphere) {
		double radius = sqrt(R2);
		if(R2_PROFILE == density_profile_) //1/r^2 density profile
		  d[i]  = density*rx*rx/(R2);
		else if(UNIFORM_DENSITY_PROFILE == density_profile_)
		  d[i] = density;
		else {
		  CkPrintf ("%s:%d %s Unknown density_profile selected\n",
			    __FILE__,__LINE__,block->name().c_str());
		  CkExit(-99);
		}
                t[i]  = temperature_;
	      }
	    }
	  }
	}
      }
    }
  }

  static int counter = 0;
  CkPrintf("Grand. Block done %d\n", counter++); 
#ifdef DEBUG_PERFORMANCE  
  if (CkMyPe()==0) {
    CkPrintf ("%s:%d %s DEBUG_PERFORMANCE %f\n",
	      __FILE__,__LINE__,block->name().c_str(),
	      timer.value());
  }
#endif  
  // Initialize particles

  Particle particle = block->data()->particle();
  
}


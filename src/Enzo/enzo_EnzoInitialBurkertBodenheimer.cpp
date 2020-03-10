// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialBurkertBodenheimer.cpp
/// @author   John Regan (john.regan@mu.ie)
/// @date     2020-04-17
/// @brief    Port of the original collapse PM infrastructure. Now used to test
///           the Burkert-Bodenheimer setup. 

#include "enzo.hpp"
#include <random>
#define BB_RADIUS 5e16  //Radius of Could Core in cm 
//----------------------------------------------------------------------

void EnzoInitialBurkertBodenheimer::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  PUParray(p,array_,3);
  p | radius_relative_;
  
}

//----------------------------------------------------------------------

void EnzoInitialBurkertBodenheimer::enforce_block
( Block * block, const Hierarchy  * hierarchy ) throw()

{

  if (!block->is_leaf()) return;

  Timer timer;
  timer.start();

  ASSERT("EnzoInitialBurkertBodenheimer",
	 "Block does not exist",
	 block != NULL);

  Field field = block->data()->field();

  // Get Field parameters
  
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  // domain extents
  double dxm,dym,dzm;
  double dxp,dyp,dzp;
  hierarchy->lower(&dxm,&dym,&dzm);
  hierarchy->upper(&dxp,&dyp,&dzp);

  // Block extents
  double bxm,bym,bzm;
  double bxp,byp,bzp;
  block->data()->lower(&bxm,&bym,&bzm);
  block->data()->upper(&bxp,&byp,&bzp);

  const int rank = cello::rank();
  
  double hx,hy,hz;
  field.cell_width(bxm,bxp,&hx,
		   bym,byp,&hy,
		   bzm,bzp,&hz);
  if (rank < 2) hy = 0.0;
  if (rank < 3) hz = 0.0;


  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  const int mx = nx + 2*gx;
  const int my = ny + 2*gy;
  const int mz = nz + 2*gz;

  const int m = mx*my*mz;

  // Get Fields

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

  const double gamma = EnzoBlock::Gamma[in];
  const double energy = 1e-3*(cello::kboltz)*temperature_ / ((gamma - 1.0) * (1.0 * cello::mass_hydrogen));


  printf("%s: gamma = %1.20lf\n", __FUNCTION__, gamma);
  printf("%s: radius_relative = %f\n", __FUNCTION__, radius_relative_); fflush(stdout);
  printf("%s: mass = %e\n", __FUNCTION__, mass_);
  // ...compute ellipsoid density

  printf("%s: dxp = %e\n", __FUNCTION__, dxp);
  printf("%s: array_[0] = %d\n", __FUNCTION__, array_[0]);
  //const double rx = (dxp - dxm) * radius_relative_ / array_[0] ;
  //const double ry = (dyp - dym) * radius_relative_ / array_[1] ;
  //const double rz = (dzp - dzm) * radius_relative_ / array_[2] ;
  
  const double rx = BB_RADIUS;
  const double ry = BB_RADIUS;
  const double rz = BB_RADIUS;
  
  const double rx2i = 1.0/(rx*rx);
  const double ry2i = 1.0/(ry*ry);
  const double rz2i = 1.0/(rz*rz);

  const double R = BB_RADIUS;
  const double DomainWidth = dxp-dxm;
  const double density = mass_ / (4.0/3.0*(cello::pi)*rx*ry*rz);
  printf("%s: Radius = %e\n", __FUNCTION__, rx);
  printf("%s: Volume = %e\n", __FUNCTION__, (4.0/3.0*(cello::pi)*rx*ry*rz));
  printf("%s: density = %e\n", __FUNCTION__, density);
  // bounds of possible explosions intersecting this Block

  printf("rx = %e\n", rx);
  printf("DomainWidth = %e\n", DomainWidth);
  printf("rx/DomainWidth = %e\n", rx/DomainWidth);
  printf("%s: array_[0] = %d\n", __FUNCTION__, array_[0]);
  printf("rx/DomainWidth*array_[0] = %e\n", rx/DomainWidth*array_[0]);
  int kxm = MAX((int)floor(rx/DomainWidth*array_[0])-1,0);
  int kym = MAX((int)floor(ry/DomainWidth*array_[1])-1,0);
  int kzm = MAX((int)floor(rz/DomainWidth*array_[2])-1,0);
  int kxp = MIN( (int)ceil(rx/DomainWidth*array_[0])+1,array_[0]);
  int kyp = MIN( (int)ceil(ry/DomainWidth*array_[1])+1,array_[1]);
  int kzp = MIN( (int)ceil(rz/DomainWidth*array_[2])+1,array_[2]);
  
  printf("%s: kxm = %d\n", __FUNCTION__, kxm);
  printf("%s: kxp = %d\n", __FUNCTION__, kxp);
  printf("DomainWidth = %e\n", DomainWidth);
  printf("%s: array_[0] = %d\n", __FUNCTION__, array_[0]);
  double hxa = DomainWidth;
  double hya = (rank >= 2) ? DomainWidth : 0.0;
  double hza = (rank >= 3) ? DomainWidth : 0.0;
  printf("hxa = %e\n", hxa);
  printf("hya = %e\n", hya);
  printf("hza = %e\n", hza);
  // (kx,ky,kz) index bounds of collapse in domain
 
  // Initialize background 

  // ratio of density inside and outside the cloud 
  const double density_ratio = 100.0;
  
  std::fill_n(d,m,density / density_ratio);
  std::fill_n(te,m,energy);
  std::fill_n(ie,m,energy);
  std::fill_n(t,m,temperature_*density_ratio);
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
  printf("kzm = %d\n", kzm);
  printf("kzp = %d\n", kzp);
  double xc = 0.0, yc = 0.0, zc = 0.0;

  printf("xc,yc,zc = %e %e %e\n", xc,yc,zc);
  // (explosion center xc,yc,zc)
  // mx,my,mz are the lengths of the blocks including ghost zones

  printf("mx,my,mz = %d %d %d\n", mx,my,mz);
  int numinsphere = 0;
  for (int iz=0; iz<mz; iz++) {   
    double z = bzm + (iz - gz + 0.5)*hz;
    for (int iy=0; iy<my; iy++) {
      double y = bym + (iy - gy + 0.5)*hy;
      for (int ix=0; ix<mx; ix++) {
	double x = bxm + (ix - gx + 0.5)*hx;
	int i = INDEX(ix,iy,iz,mx,my);
	//printf("ix,iy,iz = %d %d %d\n", ix, iy, iz);
	//printf("x,y,z = %e %e %e\n\n", x, y, z);
	double r2 = x*x*rx2i + y*y*ry2i + z*z*rz2i;
	//printf("r2 = %lf\n", r2);
	bool in_sphere = (r2 < 1.0);

	if (in_sphere) {
	  //printf("r2 = %e\n", r2);
	  //printf("in_sphere = %d\n", in_sphere);
	  //printf("Inside!!!!! temperature = %f\n", temperature_);
	  d[i]  = density;
	  t[i]  = temperature_;
	  numinsphere++;
	}
      }
    }
  }
  printf("numinsphere = %d\t total = %d\t Precentage = %lf\%\n", numinsphere, m, float(numinsphere)/float(m)*100.0);
  
  //exit(-99);
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


/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodStarMaker.cpp
/// @author     Andrew Emerick (aemerick11@gmail.com)
/// @author     Stefan Arridge (stefan.arridge@gmail.com)
/// @date
/// @brief  Implements a star maker class
///
///       This is supposed to be a general class for star formation routines
///       where individual SF routines can be created as derived classes.
///       The intention is to try and improve upon the clutter / mess
///       in Enzo's Grid_StarParticleHandler. This will do this, but
///       will still require quite a bit of repeated code across
///       individual SF (the derived classes) routines... so not perfect...

#include "cello.hpp"
#include "enzo.hpp"

// #define SHU_COLLAPSE
// #define DEBUG_SF

//-------------------------------------------------------------------

EnzoMethodStarMaker::EnzoMethodStarMaker
()
  : Method()
{
  ASSERT("EnzoMethodStarMaker::EnzoMethodStarMaker()",
	 "Cannot use star_maker method with rank < 3",cello::rank() == 3);
  const EnzoConfig * enzo_config = enzo::config();
  const EnzoUnits * enzo_units = enzo::units();
  
  // AJE: This was the old way this was done
  // Initialize default Refresh object
  // const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
  //                           enzo_sync_id_method_star_maker);
  // refresh(ir)->add_all_fields();

  cello::simulation()->new_refresh_set_name(ir_post_,name());

  Refresh * refresh = cello::refresh(ir_post_);
  ParticleDescr * particle_descr = cello::particle_descr();
  refresh->add_particle(particle_descr->type_index("star"));
  refresh->add_all_fields();

  // Copy over parameters from config to local names here for convenience
  check_number_density_threshold_ =
    enzo_config->method_star_maker_check_number_density_threshold;
  check_converging_flow_=
    enzo_config->method_star_maker_check_converging_flow;
  check_jeans_density_ =
    enzo_config->method_star_maker_check_jeans_density;
  jeans_density_factor_ =
    enzo_config->method_star_maker_jeans_density_factor;
  use_dynamical_time_        = enzo_config->method_star_maker_use_dynamical_time;
  check_self_gravitating_      =
    enzo_config->method_star_maker_check_self_gravitating;
  use_h2_self_shielding_     = enzo_config->method_star_maker_use_h2_self_shielding;
  check_jeans_mass_            = enzo_config->method_star_maker_check_jeans_mass;
  check_potential_minimum_ = enzo_config->method_star_maker_check_potential_minimum;
  check_density_maximum_ = enzo_config->method_star_maker_check_density_maximum;
  control_volume_cells_min_ =
    enzo_config->method_star_maker_control_volume_cells_min;
  control_volume_cells_max_ =
    enzo_config->method_star_maker_control_volume_cells_max;
  number_density_threshold_  =
    enzo_config->method_star_maker_number_density_threshold;
  efficiency_                = enzo_config->method_star_maker_efficiency;
  maximum_star_fraction_     = enzo_config->method_star_maker_maximum_mass_fraction;
  star_particle_min_mass_    = enzo_config->method_star_maker_minimum_star_mass;
  star_particle_max_mass_    = enzo_config->method_star_maker_maximum_star_mass;
  gamma_ = enzo_config->field_gamma;
  ggm1_ = gamma_ * (gamma_ - 1.0);
  grav_constant_internal_units_ = cello::grav_constant * enzo_units->mass()
                                  * enzo_units->time() * enzo_units->time() /
                                  (enzo_units->length() * enzo_units->length() *
				   enzo_units->length());
}

//-------------------------------------------------------------------

void EnzoMethodStarMaker::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | check_number_density_threshold_;
  p | check_jeans_density_;
  p | jeans_density_factor_;
  p | check_converging_flow_;
  p | use_dynamical_time_;
  p | check_number_density_threshold_;
  p | check_potential_minimum_;
  p | check_density_maximum_;
  p | control_volume_cells_min_;
  p | control_volume_cells_max_;
  p | efficiency_;
  p | maximum_star_fraction_;
  p | star_particle_min_mass_;
  p | star_particle_max_mass_;
  p | check_self_gravitating_;
  p | use_h2_self_shielding_;
  p | check_jeans_mass_;

  return;
}

//------------------------------------------------------------------
//   This does nothing at the moment - business is done in derived
//   class (Currently EnzoMethodStarMakerStochasticSF)
void EnzoMethodStarMaker::compute ( Block *block) throw()
{

  if (! block->is_leaf()) return;

  block->compute_done();

  return;
}

// Required
double EnzoMethodStarMaker::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}

void EnzoMethodStarMaker::rescale_densities(EnzoBlock * enzo_block,
                                            const int index,
                                            const double density_ratio) throw() {

  // Loop through all passive scalars (color fields)
  // which are mass fractions stored as densities, and rescale
  // to the new density after star formation.
  //
  // AE: NOTE: Change this routine if / whenever there needs to be
  //           fraction fields that are not labelled as color
  //           Obviously requires these fields to be declared as color
  //           in input file to work.
  //           This can / should likely get moved to be a more general
  //           function of the block / data / field class (one of those)
  //
  //    density_ratio = new_density / old_density
  //

  Field field = enzo_block->data()->field();

  Grouping * field_groups = field.groups();
  int nc = field_groups->size("color");

  for (int ic = 0; ic < nc; ic++){
    enzo_float * cfield = (enzo_float *)
      field.values(field_groups->item("color",ic));

    cfield[index] *= density_ratio;

  }

  return;
}

/*
void EnzoMethodStarMaker::convert_densities_to_fraction(EnzoBlock * enzo_block,
                                                        int direction) throw() {

  // Actually, I don't really think we need this...
  //   this only needs to be done with cells that either have
  //   star formation or get gas removed for star formation. This is
  //   likely to be a small number of cells on a given grid, so
  //   really there is no need to do this conversion for every cell...


  Field field = enzo_block->data()->field();

  int gx,gy,gz;
  field.ghost_depth (0, &gx, &gy, &gz);

  int mx, my, mz;
  field.dimensions (0, &mx, &my, &mz);

  int nx, ny, nz;
  field.size ( &nx, &ny, &nz);

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;

  enzo_float * d = (enzo_float *) field.values("density");

  if (direction == 0){ // convert density to fraction

    for (int iz = 0; iz < ngz; iz++){
      for (int iy = 0; iy < ngy; iy++){
        for (int ix = 0; ix < ngx; ix++){
          int i = INDEX(ix,iy,iz,ngx,ngy);

          double inv_dens = 1.0 / d[i];

          if (metal) metal[i] = metal[i] * inv_dens;
        }
      }
    }


  } else { // convert fraction to density

    if (metal) metal[i] *= d[i];
  }

  return;
}
*/

// ---------------------------------------------------------

int EnzoMethodStarMaker::check_number_density_threshold(
                                                       const double &d
                                                        ){

  ///  Apply the criteria that the local number density be greater
  ///  than the provided number density if check_density_threshold_ is
  ///  desired by the user.

  return !(check_number_density_threshold_) +
          (d >= number_density_threshold_);
}

int EnzoMethodStarMaker::check_self_gravitating(
                const double mean_particle_mass, const double rho_cgs,
		const enzo_float temperature,
                enzo_float *vx, enzo_float *vy, enzo_float *vz,
                const double lunit, const double vunit,
                const int &index, const int &dix, const int &diy, const int &diz,
                const double dx, const double dy, const double dz)
{

  if (!check_self_gravitating_)
    return 1;

  // Hopkins et al. (2013). Virial parameter: alpha < 1 -> self-gravitating

  double div_v_norm2, cs2, alpha;
  double dx2 = dx*dx * lunit*lunit;
  double dy2 = dy*dy * lunit*lunit;
  double dz2 = dz*dz * lunit*lunit;

  // Frobenius norm of the velocity gradient tensor
  div_v_norm2 = (pow(vx[index+dix] - vx[index-dix], 2) +
                 pow(vy[index+dix] - vy[index-dix], 2) +
                 pow(vz[index+dix] - vz[index-dix], 2)) / dx2 +
                (pow(vx[index+diy] - vx[index-diy], 2) +
                 pow(vy[index+diy] - vy[index-diy], 2) +
                 pow(vz[index+diy] - vz[index-diy], 2)) / dy2 +
                (pow(vx[index+diz] - vx[index-diz], 2) +
                 pow(vy[index+diz] - vy[index-diz], 2) +
                 pow(vz[index+diz] - vz[index-diz], 2)) / dz2;
  div_v_norm2 *= (vunit*vunit);

  // constant for testing. TODO: change to variable  
  cs2 = (gamma_ * cello::kboltz * temperature) / mean_particle_mass;

  alpha = (div_v_norm2 + cs2/dx2) / (8 * cello::pi * cello::grav_constant * rho_cgs);
  return (alpha < 1);

}

double EnzoMethodStarMaker::h2_self_shielding_factor(
                enzo_float *rho, const double metallicity,
                const double dunit, const double lunit,
                const int &index, const int &dix, const int &diy, const int &diz,
                const double dx, const double dy, const double dz)
{

  if (!use_h2_self_shielding_)
    return 1;

  // Hopkins et al. (2017) and Krumholz & Gnedin (2011). Constant numbers come from their models and fits.
  // Mass fraction that is self-shielded and able to cool. f_shield > 0

  double tau, phi, psi, f_shield, grad_rho;

  const double rho_cgs = rho[index] * dunit;

  grad_rho = sqrt(pow((rho[index+dix] - rho[index-dix]) / dx, 2) +
                  pow((rho[index+diy] - rho[index-diy]) / dy, 2) +
                  pow((rho[index+diz] - rho[index-diz]) / dz, 2));
  grad_rho *= dunit / lunit;
  tau = 434.8 * rho_cgs * (dx + rho_cgs / grad_rho);  // 434.8 cm^2 / g
  phi = 0.756 * pow(1.0 + 3.1 * metallicity, 0.365);
  psi = (0.6 * tau * (0.01 + metallicity)) / (log(1.0 + 0.6*phi + 0.01*phi*phi));
  f_shield = 1.0 - 3.0 / (1.0 + 4.0*psi);
  return f_shield;

}

/*
 * Check if the local Jeans is resolved by 1/jeans_density_factor cell widths, which
 * is equivalent to checking if cell density is larger than some threshold
 * density (see Krumholz+ 2004, ApJ, 611, 399).
 * Modifies the value of jeans_density
 */
int EnzoMethodStarMaker::check_jeans_density(const double specific_internal_energy,
					     const double mean_cell_width,
					     const double cell_density,
					     double* jeans_density)
{

  if (!check_jeans_density_)
    return 1;

  const EnzoUnits * enzo_units = enzo::units();
  const double cs2 = ggm1_ * specific_internal_energy;  

  *jeans_density  = jeans_density_factor_ * jeans_density_factor_ *
                    cello::pi * cs2 /
                    (grav_constant_internal_units_ *
                    mean_cell_width * mean_cell_width);
  
  return (cell_density > *jeans_density);
}



 
int EnzoMethodStarMaker::check_jeans_mass(
  const double temperature, const double mean_particle_mass,
  const double rho_cgs, const double mass
)
{
  if (!check_jeans_mass_)
    return 1;


  const double minimum_jeans_mass = 1000 * cello::mass_solar;
  double cs2 = (gamma_ * cello::kboltz * temperature) / mean_particle_mass;
  double m_jeans = (cello::pi/6) * pow(cs2, 1.5) / (pow(cello::grav_constant, 1.5) * sqrt(rho_cgs));
  double m_jcrit = MAX(minimum_jeans_mass, m_jeans);
  return (mass < m_jcrit);
}

// This function implements the converging flow condition for star formation.
// This is done by computing the symmetrised grad velocity tensor (or strain tensor)
// a_{ij} = 0.5*(dv_i/dx_j + dv_j/dx_i), then first checking its trace is negative
// (i.e. the velocity divergence is negative), then if this is satisfied we check
// if all the eigenvalues are negative

int EnzoMethodStarMaker::check_converging_flow(
                enzo_float *vx, enzo_float *vy, enzo_float *vz,
                const int &index, const int &dix, const int &diy,
                const int &diz, const double &dx, const double &dy,
		const double &dz){

    if (!check_converging_flow_){
      return 1;
    }
    
    const double a_11 = (vx[index+dix] - vx[index-dix])/dx;
    const double a_22 = (vy[index+diy] - vy[index-diy])/dy;
    const double a_33 = (vz[index+diz] - vz[index-diz])/dz;

    /// If trace is positive, at least one of the eigenvalues is
    /// positive, so cell fails the test
    
    if (a_11 + a_22 + a_33 > 0) return 1;

    const double a_12 = 0.5 * ( (vx[index+diy] - vx[index-diy]) / dy
			       +(vy[index+dix] - vy[index-dix]) / dx );
    const double a_13 = 0.5 * ( (vx[index+diz] - vx[index-diz]) / dz
			       +(vz[index+dix] - vz[index-dix]) / dx );
    const double a_23 = 0.5 * ( (vy[index+diz] - vy[index-diz]) / dz
			       +(vz[index+diy] - vz[index-diy]) / dy );

    // Set the coefficients of the cubic equation which gives the
    // eigenvalues, i.e. lambda^3 + A*lambda^2 + B*lambda + C = 0

    const double A = -1.0 * (a_11 + a_22 + a_33);
    const double B = -1.0 * (  a_11 * a_22 + a_11 * a_33
			     + a_22 * a_33 + a_12 * a_12
			     + a_23 * a_23 + a_13 * a_13 );
    const double C = a_11 * a_23 * a_23 + a_22 * a_13 * a_13
                   + a_33 * a_12 * a_12 - a_11 * a_22 * a_33;

    // Equation can be transformed to the form t^3 + beta*t + gamma = 0
    // where t = lambda - alpha, with alpha, beta, gamma defined as
    // follows

    const double alpha = -3.0 * A;
    const double beta  = B - A * A / 3.0;
    const double gamma = C + 2.0 * A * A * A / 27.0 - A * B / 3.0;

    /// Can transform this to a trigonometric equation by taking
    /// t = 2*sqrt(-beta/3) * cos(theta). See
    /// https://en.wikipedia.org/wiki/Cubic_equation for derivation

    /// We check if any of the eigenvalues are positive

    for (int i = 0; i < 3; i++){
      const double lambda_k = alpha + 2.0 * sqrt(-1.0 * beta / 3.0) *
	cos( acos(1.5 * gamma / beta * sqrt(-3.0 / beta)) / 3.0
	     + i * 2.0 * cello::pi / 3.0);
      if (lambda_k > 0) return 0;
    }
    
    return 1;
}


/*
 * This checks whether the given cell is at a potential minimum
 * within a spherical control volume. The radius of the control volume is the Jeans
 * length, but bounded below by min_control_volume_radius_ (which must be at least
 * one cell width) and above by  max_control_volume_radius_ (which must be no greater
 * than the ghost zone depth). Currently assumes that the spatial resolution is 
 * the same in each dimension.  
 */

/*
 * Bleuler and Teyssier 2014, MNRAS, 445, 4015 have this to say:
 * '' A local minimum in the gravitational potential is not a pre-requisite 
 *    for local gravitational collapse. This can be seen in a thought experiment
 *    where a constant force field is applied to the region of interest. The 
 *    addition of a constant force term corresponds to adding a linear term in 
 *    the gravitational potential. This changes the position and/or existence of
 *    local extrema in the potential without changing the local dynamics. This
 *    demonstrates why the tidal tensor, which is not affected by the addition of
 *    a linear term, is the right quantity for the evaluation of local 
 *    gravitational collapse (see Section 2.2). ''
*/
int EnzoMethodStarMaker::check_potential_minimum(
		       			EnzoBlock * enzo_block,
			       		const int ix, const int iy, const int iz)
{

  if (!check_potential_minimum_)
    return 1;
  Field field = enzo_block->data()->field();
  int mx, my, mz;
  field.dimensions (0, &mx, &my, &mz);
  double dx, dy, dz;
  enzo_block->cell_width(&dx, &dy, &dz);

  enzo_float * phi = (enzo_float *) field.values("potential");
  enzo_float * u = (enzo_float *) field.values("internal_energy");
  enzo_float * d = (enzo_float *) field.values("density");
  const int cell_index = INDEX(ix,iy,iz,mx,my);
  const double jeans_length = sqrt(cello::pi * ggm1_ * u[cell_index] /
			          (grav_constant_internal_units_ * d[cell_index]));
  
  /* The potential field is defined so that it has positive values, so to find
     the potential 'minimum' we actually find the positive value
   */

  /* Note: This assumes spatial resolution is the same in all dimensions */
  const double control_volume_radius = std::min(control_volume_cells_max_ * dx,
					  std::max(jeans_length,
					       control_volume_cells_min_ * dx));

  /* Convert this to a number of cells */
  const int control_volume_cells = floor(control_volume_radius / dx);

  /* We loop over cells within a cube, centred on the 'test cell' with side 
     length 2 * control_volume_cells. For each cell, we check whether it lies within
     the control volume, and if so, we check if the potential has a value larger than
     the local potential, and return 0 if so */
  int result = 1;
   for (int jz = iz - control_volume_cells; jz < iz + control_volume_cells; jz++){
     for (int jy = iy - control_volume_cells; jy < iy + control_volume_cells; jy++){
       for (int jx = ix - control_volume_cells; jx < ix + control_volume_cells; jx++){
	 const int j = INDEX(jx,jy,jz,mx,my);
	  
	 const double distance2 = dx * dx * (jx - ix) * (jx - ix) +
	                          dy * dy * (jy - iy) * (jy - iy) +
	                          dz * dz * (jz - iz) * (jz - iz);
	  if(distance2 <= control_volume_radius * control_volume_radius) {
	    if (phi[j] > phi[cell_index]){
	      result = 0;
	      break; //break out of the jx loop
	    }
	  }
       } // jx
       if (!result) break; // break out of the jy loop
     } // jy
     if (!result) break; // break out of the jz loop
   } // jz
		    
   return result;
   
}



/*
 * This checks whether the value of the density in a given cell is at a minimum
 * within a spherical control volume. The radius of the control volume is the Jeans
 * length, but bounded below by min_control_volume_radius_ (which must be at least
 * one cell width) and above by  max_control_volume_radius_ (which must be no greater
 * than the ghost zone depth). Currently assumes that the spatial resolution is 
 * the same in each dimension.  
 */

int EnzoMethodStarMaker::check_density_maximum(
		       			EnzoBlock * enzo_block,
			       		const int ix, const int iy, const int iz)
{

  if (!check_density_maximum_)
    return 1;

  // Get some data about the fields
  Field field = enzo_block->data()->field();
  int mx, my, mz;
  field.dimensions (0, &mx, &my, &mz);
  double dx, dy, dz;
  enzo_block->cell_width(&dx, &dy, &dz);
  const int cell_index = INDEX(ix,iy,iz,mx,my);
  
  // Get pointers to relevant fields
  enzo_float * u = (enzo_float *) field.values("internal_energy");
  enzo_float * d = (enzo_float *) field.values("density");

  const double jeans_length = sqrt(cello::pi * ggm1_ * u[cell_index] /
			          (grav_constant_internal_units_ * d[cell_index]));

  /* Note: This assumes spatial resolution is the same in all dimensions */
  const double control_volume_radius = std::min(control_volume_cells_max_ * dx,
					  std::max(jeans_length,
					       control_volume_cells_min_ * dx));

  /* Convert this to a number of cells */
  const int control_volume_cells = floor(control_volume_radius / dx);

  int result = 1;
  /* We loop over cells within a cube, centred on the 'test cell' with side 
     length 2 * control_volume_cells. For each cell, we check whether it lies within
     the control volume, and if so, we check if the potential has a value larger than
     the local potential, and return 0 if so */
   for (int jz = iz - control_volume_cells; jz < iz + control_volume_cells; jz++){
     for (int jy = iy - control_volume_cells; jy < iy + control_volume_cells; jy++){
       for (int jx = ix - control_volume_cells; jx < ix + control_volume_cells; jx++){
	 const int j = INDEX(jx,jy,jz,mx,my);

	 const double distance2 = dx * dx * (jx - ix) * (jx - ix) +
	                          dy * dy * (jy - iy) * (jy - iy) +
	                          dz * dz * (jz - iz) * (jz - iz);
	  if(distance2 <= control_volume_radius * control_volume_radius) {
	    if (d[j] > d[cell_index]){
	      result = 0;
	      break; // break out of the jx loop
	    }
	  }
       } // jx
       if (!result) break; // break out of the jy loop
     } // jy
     if (!result) break; // break out of the jz loop
   }
   return result;
   
}


int EnzoMethodStarMaker::check_mass(const double &m){
  /// Apply the condition that the mass of gas converted into
  /// stars in a single cell cannot exceed a certain fraction
  /// of that cell's mass. There does not need to be a check on
  /// the maximum particle mass.

  int minlimit = ((maximum_star_fraction_ * m) > star_particle_min_mass_);
  return minlimit;

}

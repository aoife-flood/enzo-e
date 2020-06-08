// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodAccretion.hpp
/// @author     John Regan (john.regan@mu.ie)
/// @date
/// @brief

#ifndef ENZO_ENZO_METHOD_ACCRETION
#define ENZO_ENZO_METHOD_ACCRETION

#define ACCRETION_LIMIT                       0.1
#define SPHERICAL_BONDI_HOYLE_FORMALISM       1
class EnzoMethodAccretion : public Method {

  /// @class   EnzoMethodAccretion
  /// @ingroup Enzo
  /// @btief   [\ref Enzo] Encapsulate Accretion Routines

public:

  EnzoMethodAccretion();

  /// Destructor
  virtual ~EnzoMethodAccretion() throw() {};

  /// CHarm++ Pup::able declarations
  PUPable_decl(EnzoMethodAccretion);

  /// Charm++ Pup::able migration Constructor
  EnzoMethodAccretion (CkMigrateMessage *m)
    : Method (m)
    {  }

  /// Charm++ Pack / Unpack function
  void pup(PUP::er &p);

  /// Apply the method
  virtual void compute (Block * block) throw();

  void compute_ (Block * block) throw();

  /// name
  virtual std::string name() throw()
  { return "accretion"; }

  // Compute the maximum timestep for this method
  virtual double timestep (Block * block) const throw();


protected: //methods
  double calculate_bondi_hoyle_radius(enzo_float pmass,
				      enzo_float *pvel,
				      enzo_float cell_temp,
				      enzo_float *cell_vel);
  double bondi_alpha(float x);

  int remove_accreted_mass(Block * block,
			   enzo_float ppos[3],
			   enzo_float pvel[3],
			   enzo_float pmass,
			   enzo_float kernel_radius,
			   enzo_float sum_of_weights,
			   enzo_float particle_accretion_radius,
			   enzo_float accretion_rate,
			   enzo_float *accreted_mass,
			   enzo_float *DeltaV);
protected: //variables

  int prescription_;
  int dual_energy_;
  // variables to be passsed here
  double gamma_;
};


#endif /* EnzoMethodAccretion */

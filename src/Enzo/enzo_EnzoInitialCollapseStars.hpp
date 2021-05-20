// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialShuCollapse.hpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     2021-05-20
/// @brief    [\ref Enzo] Initial conditions for a uniform density sphere of stars

#ifndef ENZO_ENZO_INITIAL_COLLAPSE_STARS_HPP
#define ENZO_ENZO_INITIAL_COLLAPSE_STARS_HPP

class EnzoInitialCollapseStars : public Initial {

  /// @class    EnzoInitialCollapseStars
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] initial conditions for singular isothermal sphere problem

public: // interface

  /// CHARM++ constructor
  EnzoInitialCollapseStars() throw()
  : Initial (),
    truncation_radius_(0.0),
    sound_speed_(0.0),
    instability_parameter_(0.0),
    central_particle_(false),
    central_particle_mass_(0.0)
  {
  }    

    /// Constructor
  EnzoInitialCollapseStars(const EnzoConfig * enzo_config) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialCollapseStars);

  /// CHARM++ migration constructor
  EnzoInitialCollapseStars(CkMigrateMessage *m)
    : Initial (m),
      truncation_radius_(0.0),
      sound_speed_(0.0),
      instability_parameter_(0.0),
      central_particle_(false),
      central_particle_mass_(0.0)
  {
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  virtual void enforce_block
  ( Block * block, Hierarchy * hierarchy ) throw();

  private: // attributes
  
  /// Location of the centre of collapse
  double centre_[3];

  /// Drift velocity of the whole system
  double drift_velocity_[3];

  /// Truncation radius - must be less than half the total domain width
  double truncation_radius_;

  /// Sound speed of the gas - internal energy calculated from this
  double sound_speed_;

  /// Instability parameter as defined in Shu 1977
  double instability_parameter_;

  /// Is there are star particle initialised at the centre of collapse?
  bool central_particle_;

  /// Mass of the central particle
  double central_particle_mass_;

  /// Physics variables we store as attribrutes for convenience
  double gamma_;
  double ggm1_;
  double grav_constant_internal_units_;
};

#endif /* ENZO_ENZO_INITIAL_COLLAPSE_HPP */


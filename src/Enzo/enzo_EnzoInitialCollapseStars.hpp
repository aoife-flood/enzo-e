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
    truncation_radius_(0.0)
  {
  }    

    /// Constructor
  EnzoInitialCollapseStars(const EnzoConfig * enzo_config) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialCollapseStars);

  /// CHARM++ migration constructor
  EnzoInitialCollapseStars(CkMigrateMessage *m)
    : Initial (m),
      truncation_radius_(0.0)
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

  /// Mass density of the sphere
  double density_;

  /// The random seed used to choose the cells hosting a star particle
  double random_seed_;

  /// Particle positions are offset by a random number, uniformly distributed
  /// between +/- 0.5 * offset_factor_ * cell_width
  double offset_factor_;
};

#endif /* ENZO_ENZO_INITIAL_COLLAPSE_HPP */


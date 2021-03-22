// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialShuCollapse.hpp
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     2021-03-21
/// @brief    [\ref Enzo] Singular isothermal sphere problem, as described by
///                       Shu (1997)

#ifndef ENZO_ENZO_INITIAL_SHU_COLLAPSE_HPP
#define ENZO_ENZO_INITIAL_SHU_COLLAPSE_HPP

class EnzoInitialShuCollapse : public Initial {

  /// @class    EnzoInitialShuCollapse
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] initial conditions for singular isothermal sphere problem

public: // interface

  /// CHARM++ constructor
  EnzoInitialShuCollapse
  (int cycle, double time,
   int rank,
   const double centre[3],
   double truncation_radius,
   double sound_speed,
   double instability_parameter,
   const double drift_velocity[3],
   bool central_particle
   ) throw()
    : Initial (cycle,time),
      rank_(rank),
      truncation_radius_(truncation_radius_),
      sound_speed_(sound_speed),
      instability_parameter_(instability_parameter),
      central_particle_(central_particle)
  {
    centre_[0] = centre[0];
    centre_[1] = centre[1];
    centre_[2] = centre[2];

    drift_velocity_[0] = drift_velocity[0];
    drift_velocity_[1] = drift_velocity[1];
    drift_velocity_[2] = drift_velocity[2];
  }    

    /// Constructor
  EnzoInitialShuCollapse(const EnzoConfig * enzo_config) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialShuCollapse);

  /// CHARM++ migration constructor
  EnzoInitialShuCollapse(CkMigrateMessage *m)
    : Initial (m),
      rank_(0),
      truncation_radius_(0.0),
      sound_speed_(0.0),
      instability_parameter_(0.0),
      central_particle_(false),
  {

    centre_[0] = centre[0];
    centre_[1] = centre[1];
    centre_[2] = centre[2];

    drift_velocity_[0] = drift_velocity[0];
    drift_velocity_[1] = drift_velocity[1];
    drift_velocity_[2] = drift_velocity[2];
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

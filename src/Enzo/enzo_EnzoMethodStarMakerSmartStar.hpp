///

///
///
///
///

#ifndef ENZO_ENZO_METHOD_STARMAKER_SMARTSTAR
#define ENZO_ENZO_METHOD_STARMAKER_SMARTSTAR

#define SS_NTIMES 5
#define CRITICAL_ACCRETION_RATE 0.04    //Msun/yr see Sakurai et al. (2016)
///
///
///
///

class EnzoMethodStarMakerSmartStar : public EnzoMethodStarMaker {

public:
  // Create new EnzoStarMakerSmartStar object
  EnzoMethodStarMakerSmartStar();

  // Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodStarMakerSmartStar);

  // Charm++ PUP::able declarations
  EnzoMethodStarMakerSmartStar (CkMigrateMessage *m)
   : EnzoMethodStarMaker (m)
   {  }

  /// Charm++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply method
  virtual void compute ( Block * block) throw();

  //virtual double timestep (Block * block) const throw();

  virtual std::string particle_type () throw()
  { return "star";}

  /// Name
  virtual std::string name () throw()
   { return "star_maker";}

  virtual ~EnzoMethodStarMakerSmartStar() throw() {};

protected:

  //int FofList(int, enzo_float *, enzo_float, int *, int **, int ***);
  
  int accretion_radius_cells_;

  
};

#endif /* EnzoMethodStarMakerSmartStar */

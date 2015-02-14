
#ifndef TUBE_CHARACTERISTICS
#define TUBE_CHARACTERISTICS

  //! systemic properties of the microtubule network
  struct tubeCharacteristcs_t
  {
    //! the volume of space in which tubes originate (nm^3)
    double volume;
	//! the mean tube length in nm
    double mean_tube_length;
	//! the mean angle between segments within a tube (deg)
    double mean_intra_tube_angle;
	//! the mean angle between tubes (deg)
    double mean_inter_tube_angle;
	//! the density of tube segments within the volume (segments/nm^3)
    double mean_tube_density;
	//! segment length (nm)
    double segLength;
	//! the total number of segments for all tubes
    size_t numSegments;
	//! the persistence length (nm)
    double persistenceLength;
	//! the number of segments per tube
    size_t segPerTube;
	//! the total number of tubes
	size_t numTubes;
	//! the structural entropy of all tubes in bits
    double se;
  };
  
  #endif /* TUBE_CHARACTERISTICS */
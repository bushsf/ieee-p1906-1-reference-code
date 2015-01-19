/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 *  Copyright � 2015 by IEEE.
 *
 *  This source file is an essential part of IEEE Std 1906.1,
 *  Recommended Practice for Nanoscale and Molecular
 *  Communication Framework.
 *  Verbatim copies of this source file may be used and
 *  distributed without restriction. Modifications to this source
 *  file as permitted in IEEE Std 1906.1 may also be made and
 *  distributed. All other uses require permission from the IEEE
 *  Standards Department (stds-ipr@ieee.org). All other rights
 *  reserved.
 *
 *  This source file is provided on an AS IS basis.
 *  The IEEE disclaims ANY WARRANTY EXPRESS OR IMPLIED INCLUDING
 *  ANY WARRANTY OF MERCHANTABILITY AND FITNESS FOR USE FOR A
 *  PARTICULAR PURPOSE.
 *  The user of the source file shall indemnify and hold
 *  IEEE harmless from any damages or liability arising out of
 *  the use thereof.
 *
 * Author: Stephen F Bush - GE Global Research
 *                      bushsf@research.ge.com
 *                      http://www.amazon.com/author/stephenbush
 */


#ifndef P1906_MOL_EXTENDED_FIELD
#define P1906_MOL_EXTENDED_FIELD

#include <iostream>
#include <fstream>
using namespace std;

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

#include "ns3/object.h"
#include "ns3/nstime.h"
#include "ns3/ptr.h"
#include "ns3/p1906-mol-field.h"
#include "ns3/p1906-mol-pos.h"

namespace ns3 {

/**
 * \ingroup P1906 framework
 *
 * \class P1906MOL_ExtendedField
 *
 * \brief Class extends the molecular Field component of the P1906 framework towards implementing a vector field
 *
 * This class implements persistence length as described in:
 *	  Bush, S. F., & Goel, S. (2013). Persistence Length as a Metric for Modeling and 
 *	    Simulation of Nanoscale Communication Networks, 31(12), 815-824. 
 *		http://dx.doi.org/10.1109/JSAC.2013.SUP2.12130014.
 *	  
 * The Gnu Scientific Library (GSL) is required for these methods.
 *  All points and positions are in three dimensions comprised of a gsl_vector * of length three (x, y, z).
 *  All lines and segments are comprised of two points within a gsl_vector * of length size (x1, y1, z1), (x2, y2, z2).
 *  Lists of points or positions are in a gsl_matrix * of size n x 3.
 *  Each tube is comprised of a list of segments within a gsl_matrix * of size s x 6 -> s x ((x1, y1, z1), (x2, y2, z2)).
 *  A set of tubes is also a gsl_matrix * of size (s * t) x 6, where s is the number of segments and t the number of tubes.
 *  All random number are derived from gsl_rng *.
 */

class P1906MOL_ExtendedField : public P1906MOLField
{
public:
  static TypeId GetTypeId (void);
    
  //! properties of the microtubule network
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
  
  //! Methods related to tube properties
  //! set the volume in which tubes will be generated
  //! \todo the setTube* functions could be specified before ns-3 starts and exist as ns-3 nodes
  void setTubeVolume(struct tubeCharacteristcs_t * ts, double volume = 25);
  //! set the mean tube length
  void setTubeLength(struct tubeCharacteristcs_t * ts, double mean_tube_length = 100);
  //! set the mean angle between tube segments
  void setTubeIntraAngle(struct tubeCharacteristcs_t * ts, double mean_intra_tube_angle = 30);
  //! set the mean angle between entire tubes
  void setTubeInterAngle(struct tubeCharacteristcs_t * ts, double mean_inter_tube_angle = 10);
  //! set the density of tube segments
  void setTubeDensity(struct tubeCharacteristcs_t * ts, double mean_tube_density = 10);
  //! set the persistence length of the tubes
  void setTubePersistenceLength(struct tubeCharacteristcs_t * ts, double persistenceLength = 50);
  //! set the number of segments per tube
  void setTubeSegments(struct tubeCharacteristcs_t * ts, size_t segPerTube = 10);
  //! display all the microtubule network properties
  void displayTubeChars(struct tubeCharacteristcs_t * ts);
  
  //! Methods related to vector fields
  //! return closest point within threshold
  void findClosestPoint(gsl_vector * pt, gsl_matrix * vf, gsl_vector *result);
  //! convert the tube structures to a vector field of the same dimensions as the tubeMatrix
  void tubes2VectorField(gsl_matrix * tubeMatrix, gsl_matrix * vf);
 
  //! Methods implementing unit tests
  //! test tube overlaps
  bool unitTestgetOverlap();
  //! test the findClosestPoint function
  bool unitTestfindClosestPoint();

    //! Methods related to creating and displaying points and lines in 3D
  //! return a pt vector comprised of x, y, z
  void point(gsl_vector * pt, double x, double y, double z);
  //! return a line comprised of two points
  void line(gsl_vector * line, gsl_vector * pt1, gsl_vector * pt2);
  //! place a line comprised of pt1 and pt2 into a list of lines at position mp
  void line(gsl_matrix * line, int mp, gsl_vector * pt1, gsl_vector * pt2);
  //! retrieve segment mp from tubeMatrix and place it in segment
  void line(gsl_vector * segment, gsl_matrix * tubeMatrix, int mp);
  //! simply print the list of points
  void displayPoints(gsl_matrix *pts);
  //! print only the first numPts
  void displayPoints(gsl_matrix *pts, size_t numPts);
  //! simply print the value of a single point pt
  void displayPoint(gsl_vector *pt);
  
  //! true if pt insects segment, false otherwise
  //! \todo this should be in Motion class
  bool isPointOverlap(gsl_vector * pt, gsl_vector * segment);
   
  //! Methods related to computing structural entropy
  //! return the structural entropy give be a list of angles
  double sEntropy(gsl_matrix *segAngle);
  //! return the cross product of u and v in product, all vectors are the same size
  void cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *product);
  //! return the shortest distance in 3D space between the point pt and the segment
  double distance(gsl_vector *pt, gsl_vector *segment);

  //! return all the points where tubes overlap with one another in pts
  void getAllOverlaps3D(gsl_matrix *tubeMatrix, vector<P1906MOL_Pos> & pts);
  //! return all the points where a segment overlaps with a list of tubes in pts
  int getOverlap3D(gsl_vector *segment, gsl_matrix *tubeMatrix, gsl_matrix *pts, gsl_vector *tubeSegments);
  //! return the nearest segment in tubeMatrix to the point pt that falls within radius, otherwise -1
  int findNearestTube(gsl_vector *pt, gsl_matrix *tubeMatrix, double radius);  
   
  P1906MOL_ExtendedField ();
  virtual ~P1906MOL_ExtendedField ();

};

}

#endif /* P1906_MOL_EXTEND_FIELD */

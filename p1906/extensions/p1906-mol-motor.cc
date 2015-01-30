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
 
/* \details
 * <pre>
 *                                                upper-right corner
 *                               +--------------------+             
 *                               |                    |             
 *                               |                    |             
 *                   ++          |    DESTINATION     |             
 *    +-->           ++          |       VOLUME       |             
 *    |           MOLECULAR      |                    |             
 *    |             MOTOR        |                    |             
 *    +                          |                    |             
 * STARTING                      +--------------------+             
 *  POINT                  lower-left corner                                             
 * </pre>
 */

#include "ns3/log.h"

#include "ns3/p1906-mol-extended-motion.h"
#include "ns3/p1906-mol-motor.h"

namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("P1906MOL_Motor");

TypeId P1906MOL_Motor::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::P1906MOL_Motor")
    .SetParent<P1906MOL_ExtendedMotion> ();
  return tid;
}

P1906MOL_Motor::P1906MOL_Motor ()
{
  /** This class implements persistence length as described in:
	  Bush, S. F., & Goel, S. (2013). Persistence Length as a Metric for Modeling and 
	    Simulation of Nanoscale Communication Networks, 31(12), 815-824. 
		http://dx.doi.org/10.1109/JSAC.2013.SUP2.12130014.
		
    The Gnu Scientific Library (GSL) is required for these methods.
    All points and positions are in three dimensions comprised of a gsl_vector * of length three (x, y, z).
    All lines and segments are comprised of two points within a gsl_vector * of length size (x1, y1, z1), (x2, y2, z2).
    Lists of points or positions are in a gsl_matrix * of size n x 3.
    Each tube is comprised of a list of segments within a gsl_matrix * of size s x 6 -> s x ((x1, y1, z1), (x2, y2, z2)).
    A set of tubes is also a gsl_matrix * of size (s * t) x 6, where s is the number of segments and t the number of tubes.
    All random number are derived from gsl_rng *.
  */
  
  //! keep track of the motor's current position
  current_location = gsl_vector_alloc (3);
  //! location of opposite corners of the destination volume
  destination_volume = gsl_vector_alloc (6);
  
  //! start with an empty record of for tracking position
  pos_history.clear();
  
  //! random number generation structures and initialization
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_env_setup();
  
}

//! this is where the motor starts, for example, the location of the transmitter
void P1906MOL_Motor::setStartingPoint(gsl_vector * pt)
{
  gsl_vector_memcpy(current_location, pt);
}

//! this is where the motor ends defined by the lower left corner and upper right corner of a cube, for example, the location of a receiver
void P1906MOL_Motor::setDestinationVolume(gsl_vector * lower_left, gsl_vector * upper_right)
{
    gsl_vector_set(destination_volume, 0, gsl_vector_get (lower_left, 0));
	gsl_vector_set(destination_volume, 1, gsl_vector_get (lower_left, 1));
	gsl_vector_set(destination_volume, 2, gsl_vector_get (lower_left, 2));
	gsl_vector_set(destination_volume, 3, gsl_vector_get (upper_right, 0));
	gsl_vector_set(destination_volume, 4, gsl_vector_get (upper_right, 1));
	gsl_vector_set(destination_volume, 5, gsl_vector_get (upper_right, 2));	
}

//! return true if motor is in the destination volume, false otherwise
bool P1906MOL_Motor::inDestination()
{
  bool inDest = false;
  
  if ((gsl_vector_get (current_location, 0) >= gsl_vector_get (destination_volume, 0)) &&
      (gsl_vector_get (current_location, 0) <= gsl_vector_get (destination_volume, 3)) &&
      (gsl_vector_get (current_location, 1) <= gsl_vector_get (destination_volume, 1)) &&
      (gsl_vector_get (current_location, 1) <= gsl_vector_get (destination_volume, 4)) &&
      (gsl_vector_get (current_location, 2) <= gsl_vector_get (destination_volume, 2)) &&
      (gsl_vector_get (current_location, 2) <= gsl_vector_get (destination_volume, 5)))
        inDest = true;
		
  return inDest;
}

//! use microtubules, if available, Brownian motion otherwise until destination is reached
void P1906MOL_Motor::move2Destination(gsl_matrix * tubeMatrix, size_t segPerTube, double timePeriod, vector<P1906MOL_Pos> & pts)
{
  int timeout = 100; //! in case motor never reaches destination
  int loops = 0; //! keep track of iterations
  
  //! \todo consider a bounding box (defined in extended Motion) around the system in which motors are reflected back into the medium
  while (!inDestination() && (loops < timeout))
  {
    //! returns the index of the segment in tubeMatrix to which the motor is bound 
    float2Tube(r, current_location, pts, tubeMatrix, timePeriod);
	setLocation(pts.back());
    //printf ("(move2Destination) current location after float2Tube\n");
	//displayLocation();
	pts.back().displayPos();
	//! walk along tube until end of tube or unbound
	motorWalk(r, current_location, pts, tubeMatrix, segPerTube);
	setLocation(pts.back());
    //printf ("(move2Destination) current location after motorWalk\n");
	//pts.back().displayPos();
	//displayLocation();
	loops++;
  }
}

//! set the current motor location to pt
void P1906MOL_Motor::setLocation(P1906MOL_Pos pt)
{
  double x, y, z;
  
  pt.getPos (&x, &y, &z);
  gsl_vector_set( current_location, 0, x);
  gsl_vector_set( current_location, 1, y);
  gsl_vector_set( current_location, 2, z);
}

//! set the current location to the 3D Cartesian coordinates
void P1906MOL_Motor::setLocation(double x, double y, double z)
{
  gsl_vector_set( current_location, 0, x);
  gsl_vector_set( current_location, 1, y);
  gsl_vector_set( current_location, 2, z);
}

//! print the current motor location
void P1906MOL_Motor::displayLocation()
{
  printf ("current location: %f %f %f\n", 
    gsl_vector_get( current_location, 0),
	gsl_vector_get( current_location, 1),
	gsl_vector_get( current_location, 2));
}

//! motor uses Brownian motion until the destination volume is reached
void P1906MOL_Motor::float2Destination(double timePeriod)
{
  gsl_vector * newPos = gsl_vector_alloc (3);
  P1906MOL_Pos Pos;
  
  Pos.setPos ( gsl_vector_get (current_location, 0),
               gsl_vector_get (current_location, 1),
			   gsl_vector_get (current_location, 2) );
  pos_history.insert (pos_history.end(), Pos);
  
  //! float until in destination volume
  while (!inDestination())
  {
	brownianMotion(r, current_location, newPos, timePeriod);
	updateTime(timePeriod);
    gsl_vector_set (current_location, 0, gsl_vector_get (newPos, 0));
	gsl_vector_set (current_location, 1, gsl_vector_get (newPos, 1));
	gsl_vector_set (current_location, 2, gsl_vector_get (newPos, 2));
	//! begin at the motor's current location
    //printf ("motor location: \n");
    //displayPos(current_location);
	
	P1906MOL_Pos Pos;
    Pos.setPos ( gsl_vector_get (current_location, 0),
                 gsl_vector_get (current_location, 1),
			     gsl_vector_get (current_location, 2) );
    pos_history.insert (pos_history.end(), Pos);
  }
}

//! return the elapsed time since the motor time was last initialized
double P1906MOL_Motor::propagationDelay()
{
  return getTime();
}

P1906MOL_Motor::~P1906MOL_Motor ()
{
  //! close down random number generator activity
  gsl_rng_free (r);
  NS_LOG_FUNCTION (this);
}

} // namespace ns3

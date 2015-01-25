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
 *  1906 Component      Motor Instantiation
 * +----------------------+-----------------------+
 * |                      |                       |
 * |    MESSAGE           |  MOTOR CARGO          |
 * |                      |                       |
 * +----------------------------------------------+
 * |                      |                       |
 * |    MESSAGE CARRIER   |  MOLECULAR MOTOR      |
 * |                      |                       |
 * +----------------------------------------------+
 * |                      |                       |
 * |    MOTION            |  BROWNIAN / WALK      |
 * |                      |                       |
 * +----------------------------------------------+
 * |                      |                       |
 * |    FIELD             |  MICROTUBULE          |
 * |                      |                       |
 * +----------------------------------------------+
 * |                      |                       |
 * |    PERTURBATION      |  MOTOR CARGO TYPE     |
 * +----------------------------------------------+
 * |                      |                       |
 * |    SPECIFICITY       |  BINDING TO TARGET    |
 * |                      |                       |
 * +----------------------+-----------------------+
 * 
 *            Motor and Microtubule Network
 *
 *     XXXXX                                             
 *         XXX XX                 XXXX                 
 *                XX X       XXXXX   XXX               
 *            _        XXXXXXXX          XXXXXX         
 *            /|        XXXXX                 XXX       
 *    +-+   /     XXXXXX    XX                 XXX     
 *    +-+      XXX           XX                  XX    
 * MOLECULAR XXX                X                  XXXXX  
 *   MOTOR                      XXXXX                  X
 *                                                     
 *            XXXXX                                    
 *                 X XX X                              
 *                       X X XX                      XX
 *                             XXXX XX XXXXXX XXXXXXX  
 *                                                     
 *                              MICROTUBULES
 * </pre>
 */
  
#include "ns3/log.h"
#include "ns3/object.h"
#include "ns3/nstime.h"
#include "ns3/ptr.h"

#include "ns3/p1906-mol-field-microtubule.h"
#include "ns3/p1906-mol-MathematicaHelper.h"
#include "ns3/p1906-mol-MATLABHelper.h"
#include "ns3/p1906-mol-metrics.h"
#include "ns3/p1906-mol-motor.h"
#include "ns3/p1906-mol-tube.h"

namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("P1906MOL_MicrotubulesField");

TypeId P1906MOL_MicrotubulesField::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::P1906MOL_MicrotubulesField")
    .SetParent<P1906MOL_ExtendedField> ();
  return tid;
}

P1906MOL_MicrotubulesField::P1906MOL_MicrotubulesField ()
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
	  
    Next steps:
	  \todo move this main code into the microtubules-example.cc
	  \todo implement Brownian motion only versus with microtubules as example: need start/end locations and to compute transit time
	  \todo implement metrics or add stubs
	  \todo document code and check into https://code.google.com/p/ieee-p1906-1-reference-code/source/browse/
	  \todo include whether the member should be in a different class, e.g. Motion
	  \todo gradient is the vector field as defined above at any point on a tube: length could be scaled to motor movement rate
	  \todo consider implementing Random Markov Field (RMF) for analysis of motor motion
	  \todo curve fit to tube segments to obtain vector field equation (vector field reconstruction) is currently done in separate Mathematica code
	  \todo implement active network programmability metric: \f$A = \int_S \int_t \Delta f(t) dt dS\f$
	  \todo relate structural entropy to energy, force, and chemical complexity
	  \todo add convection-diffusion using vector field lines
	  \todo could grow tubes from alternating ends so that starting point is in the center
	  \todo return propagation time: need to keep track of time
	  \todo consider many motors operating simultaneously: they would need to update simultaneously
	  
	  \todo plot motor starting in random location versus time to reach destination as function of tube structure
	  \todo plot structural entropy versus delay
	  \todo plot binding time versus delay
	  \todo plot distance travelled versus delay, structural entropy, etc.
  */
  P1906MOL_MathematicaHelper mathematica;
  
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_env_setup();
    
  //! set the microtubule network properties
  setTubeVolume(&ts, 25);
  setTubeLength(&ts, 100);
  setTubeIntraAngle(&ts, 30);
  setTubeInterAngle(&ts, 10);
  setTubeDensity(&ts, 10);
  setTubePersistenceLength(&ts, 50);
  setTubeSegments(&ts, 10);
 
  //! display all the microtubule network properties
  displayTubeChars(&ts);
  
  //! create the microtubules
  genTubes(&ts);
  mathematica.tubes2Mma(tubeMatrix, ts.segPerTube, "tubes.mma");
  printf ("completed tube creation\n");

  //! create the vector field  
  vf = gsl_matrix_alloc (ts.numTubes * ts.segPerTube, 6);
  tubes2VectorField(tubeMatrix, vf);

  //! test the computation of distance
  unitTest_Distance();

  //! test computation of a single overlap
  unitTest_Overlap();

  //! test computation of all segment overlaps
  unitTest_AllOverlaps();

  //! test the computation of a vector field
  unitTest_VectorField();

  //! test persistence length versus entropy plot
  unitTest_PersistenceLengthsVsEntropy();

  vector<P1906MOL_Pos> pts;
  //! test the movement of a motor floating to a tube and walking on the tube
  unitTest_MotorMovement(pts);
  
  //! test plot2mma
  unitTest_Plot2Mma(pts);
  
  //! test motor movement to destination
  //unitTest_MotorMove2Destination(pts);

  NS_LOG_FUNCTION (this);
  NS_LOG_FUNCTION (this << "Created MOL Field Component");
}

//! plot persistence length versus structural entropy
void P1906MOL_MicrotubulesField::persistenceVersusEntropy(struct tubeCharacteristcs_t * ts, gsl_rng * r, gsl_vector * persistenceLengths)
{
  P1906MOL_MathematicaHelper mathematica;
  char plot_filename[256];
  gsl_matrix * tubeMatrix = gsl_matrix_alloc (ts->numTubes * ts->segPerTube, 6);
  gsl_matrix * pve = gsl_matrix_alloc (persistenceLengths->size, 2);
  
  for (size_t i = 0; i < persistenceLengths->size; i++)
  {
    setTubePersistenceLength (ts, gsl_vector_get (persistenceLengths, i));
    genTubes(ts);
	
	sprintf (plot_filename, "tubes_%ld.mma", i);
    mathematica.tubes2Mma(tubeMatrix, ts->segPerTube, plot_filename);
	gsl_matrix_set (pve, i, 0, gsl_vector_get (persistenceLengths,i));
	gsl_matrix_set (pve, i, 1, ts->se);
  }
  
  //! plot the results
  mathematica.plot2Mma(pve, "persistenceVersusEntropy.mma", "persistence length", "structural entropy");
}

//! generate a set of tubes comprised of a total of numSegments in volume with segPerTube segments of segLength 
//! and persistenceLength and return in tubeMatrix
void P1906MOL_MicrotubulesField::genTubes(struct tubeCharacteristcs_t * ts)
{
  //! create a given density of tubes of numSegments in given volume
  //! volume starts at 0, 0, 0 to volume^(1/3) in each dimension  
  tubeMatrix = gsl_matrix_alloc (ts->numTubes * ts->segPerTube, 6);

  //! \todo get actual tube graph properties from biologist
  gsl_vector * startPt = gsl_vector_alloc (3);
  //! hold the values for a tube comprised of many segments: x_start y_start x_start x_end y_end z_end
  gsl_matrix * segMatrix = gsl_matrix_alloc (ts->segPerTube, 6);
  double total_structural_entropy = 0;
  P1906MOL_Tube tube;
  
  // printf ("(genTubes) numTubes: %ld segPerTube: %ld volume: %f\n", ts->numTubes, ts->segPerTube, ts->volume);
  
  //! volume starts at 0, 0, 0 to volume^(1/3) in each dimension
  for(size_t i = 0; i < ts->numTubes; i++)
  {
    //! set the starting location for the tube
    point (startPt, 
	  gsl_ran_gaussian (r, pow(ts->volume, (1/3))),
	  gsl_ran_gaussian (r, pow(ts->volume, (1/3))),
	  gsl_ran_gaussian (r, pow(ts->volume, (1/3))));
	 
    //! create a single tube of many segments
    tube.genTube(ts, r, segMatrix, startPt);
    // displayTube(segMatrix);
	total_structural_entropy += ts->se;
	
	//! copy tube segments to main tube matrix
	for(size_t j = 0; j < ts->segPerTube; j++)
	  for(size_t k = 0; k < 6; k++)
	  {
	    gsl_matrix_set(tubeMatrix, i * ts->segPerTube + j, k, gsl_matrix_get(segMatrix, j, k));
	  }
  }
  
  ts->se = total_structural_entropy;
}

//! test distance calculation
bool P1906MOL_MicrotubulesField::unitTest_Distance()
{
  gsl_vector * startPt = gsl_vector_alloc (3);
  gsl_vector * segment = gsl_vector_alloc (6);
  gsl_vector * pt1 = gsl_vector_alloc (3);
  gsl_vector * pt2 = gsl_vector_alloc (3);
  
  printf ("beginning unitTest_Distance\n");
  point (startPt, 0, 0, 0);
  point (pt1, -1, -1, -1);
  point (pt2, 2, 2, 2);
  line (segment, pt1, pt2);
  double d = distance (startPt, segment);
  printf ("distance: %f\n", d);
  printf ("completed unitTest_Distance\n");
  return true;
}

//! test segment overlap - move to a unit test method
bool P1906MOL_MicrotubulesField::unitTest_Overlap()
{
  gsl_vector * segment3D = gsl_vector_alloc (6);
  gsl_matrix * pts3D = gsl_matrix_alloc (1, 3);
  gsl_matrix * tubeMatrix3D = gsl_matrix_alloc (1, 6);
  gsl_matrix_set_zero (pts3D);
  gsl_vector * pt1 = gsl_vector_alloc (3);
  gsl_vector * pt2 = gsl_vector_alloc (3);
  gsl_vector * pt3 = gsl_vector_alloc (3);
  gsl_vector * pt4 = gsl_vector_alloc (3);
  gsl_vector * tubeSegments = gsl_vector_alloc (1);
  
  printf ("beginning unitTest_Overlap\n");
  point (pt1, 0, 0, 0);
  point (pt2, 5, 5, 0);
  point (pt3, 5, 0, 0);
  point (pt4, 0, 5, 0);
  line (segment3D, pt1, pt2);
  line (tubeMatrix3D, 0, pt3, pt4);
  getOverlap3D(segment3D, tubeMatrix3D, pts3D, tubeSegments);
  if (isPointOverlap(pt1, segment3D))
    printf ("point overlaps\n");
  printf ("completed unitTest_Overlap\n");
  return true;
}

//! test finding all segment overlaps
bool P1906MOL_MicrotubulesField::unitTest_AllOverlaps()
{
  vector<P1906MOL_Pos> pts;
  P1906MOL_MathematicaHelper mathematica;

  printf ("beginning unitTest_AllOverlaps\n");
  getAllOverlaps3D(tubeMatrix, pts);
  mathematica.points2Mma(pts, "pfile.mma");
  printf ("completed unitTest_AllOverlaps\n");
  return true;
}

//! plot persistence length versus entropy
bool P1906MOL_MicrotubulesField::unitTest_PersistenceLengthsVsEntropy()
{
  gsl_vector * persistenceLengths = gsl_vector_alloc (10);
  
  printf ("beginning unitTest_PersistenceLengthsVsEntropy\n");
  for (size_t i = 0; i < 10; i++)
    gsl_vector_set (persistenceLengths, i, i * 100);
  persistenceVersusEntropy(&ts, r, persistenceLengths);
  printf ("completed unitTest_PersistenceLengthsVsEntropy\n");

  return true;
}

//! test plot2mma
bool P1906MOL_MicrotubulesField::unitTest_Plot2Mma(vector<P1906MOL_Pos> pts)
{
  P1906MOL_MathematicaHelper mathematica;
  P1906MOL_Pos Pos;
  double x, y, z;
  
  printf ("beginning unitTest_Plot2Mma\n");
  //! test plotting points - move to a unit test
  gsl_matrix * vals = gsl_matrix_alloc (pts.size(), 2);
  for (size_t i = 0; i < pts.size(); i++)
  {
    Pos = pts.at(i);
	Pos.getPos (&x, &y, &z);
    //for (size_t j = 0; j < 2; j++)
    gsl_matrix_set (vals, i, 0, x);
	gsl_matrix_set (vals, i, 1, y);
  }
  mathematica.plot2Mma(vals, "plottest.mma", "x value", "y value");
  printf ("completed unitTest_Plot2Mma\n");
  return true;
}

//! test tubes2VectorField
bool P1906MOL_MicrotubulesField::unitTest_VectorField()
{
  P1906MOL_MATLABHelper matlab;
  P1906MOL_MathematicaHelper mathematica;

  printf ("beginning unitTest_VectorField\n");
  mathematica.vectorFieldPlotMma(vf, "vectorField.mma");
  matlab.vectorFieldMeshMATLAB(vf, "vectorField.dat");
  printf ("completed plot of vector field\n");
  return true;
}

//! test motor movement to a tube and walking along the tube
bool P1906MOL_MicrotubulesField::unitTest_MotorMovement(vector<P1906MOL_Pos> & pts)
{
  //! \todo modify this to start from transmitter and end when motor reaches receiver returning propagation time
  //! random motion for 100 sec
  //! start near the first tube
  P1906MOL_Motor motor;
  P1906MOL_MathematicaHelper mathematica;
  gsl_vector * startPt = gsl_vector_alloc (3);
  gsl_vector * ll = gsl_vector_alloc (3);
  gsl_vector * ur = gsl_vector_alloc (3);
  double timePeriod = 100;

  printf ("beginning unitTest_MotorMovement\n");
  //! reset the motor's timer
  motor.initTime();
  
  //! start at zero
  point (startPt, 0, 0, 0);
  //! the destination is anywhere in the volume 100...1000 along each dimension
  point (ll, 1000, 1000, 1000);
  point (ur, 2000, 2000, 2000);

  pts.clear(); //! reset the position history
  motor.setStartingPoint(startPt);
  motor.setDestinationVolume(ll, ur);
  motor.float2Destination(timePeriod);
  printf ("(unitTest_MotorMovement) propagation time: %f\n", motor.getTime());
  mathematica.connectedPoints2Mma(motor.pos_history, "motion.mma");
  printf ("completed float2Destination\n");

  //! reset the motor's timer
  motor.initTime();
  
  //! move randomly until overlap with tube
  pts.clear();
  point (startPt, 
    gsl_matrix_get (tubeMatrix, 0, 0) + 10,
	gsl_matrix_get (tubeMatrix, 0, 1),
	gsl_matrix_get (tubeMatrix, 0, 2));
  motor.float2Tube(r, startPt, pts, tubeMatrix, 0.1);
  printf ("(float2Tube) propagation time: %f\n", motor.getTime());
  printf ("completed float to tube\n");	
  mathematica.connectedPoints2Mma(pts, "motion2tube.mma");
  printf ("completed unitTest_MotorMovement\n");
  
  //! start where the motor ended
  P1906MOL_Pos Pos;
  Pos = pts.back();
  double x, y, z;
  Pos.getPos (&x, &y, &z);
  point(startPt, x, y, z);
  displayPoint (startPt);

  //! now walk along the tube
  pts.clear();
  motor.motorWalk(r, startPt, pts, tubeMatrix, ts.segPerTube);
  printf ("(motorWalk) propagation time: %f\n", motor.getTime());
  mathematica.connectedPoints2Mma(pts, "tubewalk.mma");
  printf ("completed unitTest_MotorMovement\n");
  return true;
}

//! test motor movement to destination
bool P1906MOL_MicrotubulesField::unitTest_MotorMove2Destination(vector<P1906MOL_Pos> & pts)
{
  P1906MOL_Motor motor;
  P1906MOL_MathematicaHelper mathematica;
  gsl_vector * startPt = gsl_vector_alloc (3);
  gsl_vector * ll = gsl_vector_alloc (3);
  gsl_vector * ur = gsl_vector_alloc (3);
  double timePeriod = 100;

  printf ("beginning unitTest_MotorMove2Destination\n");
  //! reset the motor's timer
  motor.initTime();
  
  //! start at zero
  point (startPt, 0, 0, 0);
  //! the destination is anywhere in the volume 100...1000 along each dimension
  point (ll, 1000, 1000, 1000);
  point (ur, 2000, 2000, 2000);

  pts.clear(); //! reset the position history
  motor.setStartingPoint(startPt);
  motor.setDestinationVolume(ll, ur);
  motor.move2Destination(tubeMatrix, ts.segPerTube, timePeriod, pts);
  printf ("(unitTest_MotorMove2Destination) propagation time: %f\n", motor.getTime());
  mathematica.connectedPoints2Mma(motor.pos_history, "motion.mma");
  printf ("completed unitTest_MotorMove2Destination\n");
  return true;
}

P1906MOL_MicrotubulesField::~P1906MOL_MicrotubulesField ()
{
  //! close down field activity
  gsl_rng_free (r);

  NS_LOG_FUNCTION (this);
}

} // namespace ns3

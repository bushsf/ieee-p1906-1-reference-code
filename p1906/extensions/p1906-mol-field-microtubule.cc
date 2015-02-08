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
 
//! \details IEEE 1906 Component map to the molecular motor instantiation
//!
//! <pre>
//!  IEEE 1906 Component        Molecular Motor 
//!                              Instantiation
//! +----------------------+-----------------------+
//! |                      |                       |
//! |    MESSAGE           |  MOTOR CARGO          |
//! |                      |                       |
//! +----------------------------------------------+
//! |                      |                       |
//! |    MESSAGE CARRIER   |  MOLECULAR MOTOR      |
//! |                      |                       |
//! +----------------------------------------------+
//! |                      |                       |
//! |    MOTION            |  BROWNIAN / WALK      |
//! |                      |                       |
//! +----------------------------------------------+
//! |                      |                       |
//! |    FIELD             |  MICROTUBULE          |
//! |                      |                       |
//! +----------------------------------------------+
//! |                      |                       |
//! |    PERTURBATION      |  MOTOR CARGO TYPE     |
//! +----------------------------------------------+
//! |                      |                       |
//! |    SPECIFICITY       |  BINDING TO TARGET    |
//! |                      |                       |
//! +----------------------+-----------------------+
//! 
//!            Motor and Microtubule Network
//!
//!     XXXXX                                             
//!         XXX XX                 XXXX                 
//!                XX X       XXXXX   XXX               
//!            _        XXXXXXXX          XXXXXX         
//!            /|        XXXXX                 XXX       
//!    +-+   /     XXXXXX    XX                 XXX     
//!    +-+      XXX           XX                  XX    
//! MOLECULAR XXX                X                  XXXXX  
//!   MOTOR                      XXXXX                  X
//!                                                     
//!            XXXXX                                    
//!                 X XX X                              
//!                       X X XX                      XX
//!                             XXXX XX XXXXXX XXXXXXX  
//!                                                     
//!                              MICROTUBULES
//! </pre>
  
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
#include "ns3/p1906-mol-vol-surface.h"

namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("P1906MOL_MicrotubulesField");

TypeId P1906MOL_MicrotubulesField::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::P1906MOL_MicrotubulesField")
    .SetParent<P1906MOL_ExtendedField> ();
  return tid;
}

/** 
  This class implements persistence length as described in:
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
	  \todo consider implementing Random Markov Field (RMF) for analysis of motor motion
	  \todo curve fit to tube segments to obtain vector field equation (vector field reconstruction) is currently done in separate Mathematica code
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
P1906MOL_MicrotubulesField::P1906MOL_MicrotubulesField ()
{
  P1906MOL_MathematicaHelper mathematica;
  
  //! allocate and start the random number generator
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
  tubeMatrix = gsl_matrix_alloc (ts.numTubes * ts.segPerTube, 6);
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

  //! test the plotting the vector field
  unitTest_VectorField();
  
  //! test the volume surface flux measurement
  unitTest_VolSurface();

  vector<P1906MOL_Pos> pts;
  //! test the movement of a motor floating to a tube and walking along the tube
  unitTest_MotorMovement(pts);
    
  //! test plot2mma to plot the position history stored in pts
  unitTest_Plot2Mma(pts);
  
  //! start at zero and utilize Brownian motion to reach destination
  unitTest_NoTubeMotion();
  
  //! test motor movement to destination using Brownian motion or microtubules if available
  unitTest_MotorMove2Destination(pts);
  
  //! test the reflective barrier volume surface
  unitTest_ReflectiveBarrier();
  
  //! test the FluxMeter
  bool unitTest_FluxMeter();
  
  //! test persistence length versus entropy plot - NB: this test changes the tubeMatrix
  unitTest_PersistenceLengthsVsEntropy();

  NS_LOG_FUNCTION (this);
  NS_LOG_FUNCTION (this << "Created MOL Field Component");
}

//! for each of the persistenceLengths in the vector, generate tubes and plot persistence length versus 
//! structural entropy in persistenceVersusEntropy.mma, also write tubes_<n>.mma for each of the persistenceLengths
void P1906MOL_MicrotubulesField::persistenceVersusEntropy(struct tubeCharacteristcs_t * ts, gsl_vector * persistenceLengths)
{
  P1906MOL_MathematicaHelper mathematica;
  char plot_filename[256];
  //! store the results here
  gsl_matrix * pve = gsl_matrix_alloc (persistenceLengths->size, 2);
  
  //printf ("(persistenceVersusEntropy) persistenceLengths->size: %ld\n", persistenceLengths->size);
  for (size_t i = 0; i < persistenceLengths->size; i++)
  {	
    //printf ("(persistenceVersusEntropy) persistenceLengths(%ld) = %f\n", i, gsl_vector_get (persistenceLengths, i));
    setTubePersistenceLength (ts, gsl_vector_get (persistenceLengths, i));
    genTubes(ts);
	
	//printf ("(persistenceVersusEntropy) i: %ld\n", i);
	//printf ("(persistenceVersusEntropy) tubeMatrix: %ld x %ld\n", tubeMatrix->size1, tubeMatrix->size2);
	//printf ("(persistenceVersusEntropy) segPerTube: %ld\n", ts->segPerTube);
	//gsl_matrix_fprintf (stdout, tubeMatrix, "%f");
	
	//! store the set of tubes
	sprintf (plot_filename, "tubes_%ld.mma", i);
    mathematica.tubes2Mma(tubeMatrix, ts->segPerTube, plot_filename);
	gsl_matrix_set (pve, i, 0, gsl_vector_get (persistenceLengths, i));
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
  //! volume starts at 0, 0, 0 to volume^(1/4) in each dimension
  
  //! don't delete, just refill the tubeMatrix
  //! free the matrix of any previous allocation
  //gsl_matrix_free (tubeMatrix);
  //! re-allocate the tubeMatrix
  //tubeMatrix = gsl_matrix_alloc (ts->numTubes * ts->segPerTube, 6);

  //! \todo get actual tube graph properties from biologist
  gsl_vector * startPt = gsl_vector_alloc (3);
  //! hold the values for a tube comprised of many segments: x_start y_start x_start x_end y_end z_end
  gsl_matrix * segMatrix = gsl_matrix_alloc (ts->segPerTube, 6);
  double total_structural_entropy = 0;
  
  //printf ("(genTubes) tubeMatrix: %ld x %ld\n", tubeMatrix->size1, tubeMatrix->size2);
  //printf ("(genTubes) numTubes: %ld segPerTube: %ld volume: %f\n", ts->numTubes, ts->segPerTube, ts->volume);
  
  //! volume starts at 0, 0, 0 to volume^(1/4) in each dimension
  for(size_t i = 0; i < ts->numTubes; i++)
  {
    //! set the starting location for the tube
    point (startPt, 
	  gsl_ran_gaussian (r, pow(ts->volume, (1/4))),
	  gsl_ran_gaussian (r, pow(ts->volume, (1/4))),
	  gsl_ran_gaussian (r, pow(ts->volume, (1/4))));

	//! create a single tube
    P1906MOL_Tube tube(ts, startPt);
	
    //! create a single tube of many segments
    //tube.genTube(ts, r, segMatrix, startPt);
    //printf ("(genTubes) segMatrix\n");
	//tube.displayTube();
	total_structural_entropy += (ts->se);
	
	tube.getSegmatrix(segMatrix);
	//! copy tube segments to main tube matrix
	for(size_t j = 0; j < ts->segPerTube; j++)
	  for(size_t k = 0; k < 6; k++)
	  {
	    gsl_matrix_set(tubeMatrix, i * ts->segPerTube + j, k, gsl_matrix_get(segMatrix, j, k));
	  }
  }
  
  ts->se = total_structural_entropy;
}

//! \todo test the volume surface as a flux meter and later as a compartmentalization volume
bool P1906MOL_MicrotubulesField::unitTest_VolSurface()
{
  P1906MOL_MathematicaHelper mathematica;
  P1906MOL_Pos c;
  c.setPos (0, 0, 0);
  P1906MOL_VolSurface vs;
  vector<P1906MOL_Pos> ipt;
  gsl_vector * segment = gsl_vector_alloc (6);
  gsl_vector * radius = gsl_vector_alloc (6);
  gsl_matrix * vectors = gsl_matrix_alloc (3, 6);
  double flux;
  
  vs.setVolume (c, 100);
  
  // test tube input/output rates relative to the volume surface
  printf ("beginning unitTest_VolSurface\n");
  
  //! segment start point
  gsl_vector_set (segment, 0, 90);
  gsl_vector_set (segment, 1, 90);
  gsl_vector_set (segment, 2, 90);
  
  //! segment end point
  gsl_vector_set (segment, 3, 110);
  gsl_vector_set (segment, 4, 110);
  gsl_vector_set (segment, 5, 110);
  
  vs.sphereIntersections(segment, ipt);
  
  //! draw the volume surface
  mathematica.volSurfacePlot(vs.center, vs.radius, "volsurface.mma");
  
  for (size_t i = 0; i < ipt.size(); i++)
  {
    printf ("(unitTest_VolSurface) intersecting point(s):\n");
	ipt.at(i).displayPos();
  }
  
  //! plot the intersection points with the volume
  mathematica.points2Mma(ipt, "volsurfintersections.mma");
  
  //! find the vector angle with the sphere surface: 
  //!   compare intersecting segment with sphere radius to determine angle
  
  //! radius start point is the same as the intersecting tube segment start point
  gsl_vector_set (radius, 0, 90);
  gsl_vector_set (radius, 1, 90);
  gsl_vector_set (radius, 2, 90);
  
  //! radius end point is at the center of the sphere
  gsl_vector_set (radius, 3, 0);
  gsl_vector_set (radius, 4, 0);
  gsl_vector_set (radius, 5, 0);
  
  double angle = vs.vectorAngle(segment, radius);
  printf ("(unitTest_VolSurface) angle: %lf\n", angle);
	
  insertVector(vectors, 0, segment);
  insertVector(vectors, 1, radius);
  
  mathematica.vectorPlotMma (vectors, "volsurfvector.mma");
  
  P1906MOL_Pos last_pos;
  P1906MOL_Pos current_pos;
  
  /*
   * reflection 1
   */
  
  last_pos.setPos (40, 40, 40);
  printf ("(unitTest_VolSurface) last_pos\n");
  last_pos.displayPos ();
  current_pos.setPos (110, 110, 110);
  printf ("(unitTest_VolSurface) current_pos\n");
  current_pos.displayPos ();
  
  line (segment, last_pos, current_pos);
  vs.sphereIntersections(segment, ipt);
  vs.getRadiusLine(radius, ipt.front()); //! assuming at least one intersection returned
  line (segment, ipt.front(), current_pos);
  
  insertVector(vectors, 0, segment);
  insertVector(vectors, 1, radius);
  line (segment, last_pos, ipt.front());
  insertVector(vectors, 2, segment);
  
  mathematica.vectorPlotMma (vectors, "volsurfvector1beforereflection.mma");
  
  vs.reflect(last_pos, current_pos);
  printf ("(unitTest_VolSurface) new current_pos\n");
  current_pos.displayPos ();
  
  //! update the reflected segment to display
  line (segment, ipt.front(), current_pos);
  insertVector(vectors, 0, segment);
  insertVector(vectors, 1, radius);
  line (segment, last_pos, ipt.front());
  insertVector(vectors, 2, segment);
  
  mathematica.vectorPlotMma (vectors, "volsurfvector1afterreflection.mma");
  
  /*
   * reflection 2
   */
  
  last_pos.setPos (20, 40, 20);
  printf ("(unitTest_VolSurface) last_pos\n");
  last_pos.displayPos ();
  current_pos.setPos (120, 90, 90);
  printf ("(unitTest_VolSurface) current_pos\n");
  current_pos.displayPos ();

  line (segment, last_pos, current_pos);
  vs.sphereIntersections(segment, ipt);
  vs.getRadiusLine(radius, ipt.front()); //! assuming at least one intersection returned
  line (segment, ipt.front(), current_pos);
  
  insertVector(vectors, 0, segment);
  insertVector(vectors, 1, radius);
  line (segment, last_pos, ipt.front());
  insertVector(vectors, 2, segment);
  
  mathematica.vectorPlotMma (vectors, "volsurfvector2beforereflection.mma");
  
  vs.reflect(last_pos, current_pos);
  printf ("(unitTest_VolSurface) new current_pos\n");
  current_pos.displayPos ();
  
  //! update the reflected segment to display
  line (segment, ipt.front(), current_pos);
  insertVector(vectors, 0, segment);
  
  line (segment, last_pos, current_pos);
  insertVector(vectors, 1, radius);
  line (segment, last_pos, ipt.front());
  insertVector(vectors, 2, segment);
  
  mathematica.vectorPlotMma (vectors, "volsurfvector2afterreflection.mma");
  
  //! measure flow through surface
  flux = vs.fluxMeter(tubeMatrix);
  printf ("(unitTest_VolSurface) flux: %lf\n", flux);
  
  printf ("completed unitTest_VolSurface\n");
  
  return true;
}

//! test the volume surface as a compartmentalization volume
bool P1906MOL_MicrotubulesField::unitTest_ReflectiveBarrier()
{
  P1906MOL_Motor motor;
  P1906MOL_MathematicaHelper mathematica;
  gsl_vector * startPt = gsl_vector_alloc (3);
  double timePeriod = 100;
  P1906MOL_Pos volCenter;

  printf ("beginning unitTest_ReflectiveBarrier\n");
  //! start at zero
  point (startPt, 0, 0, 0);
  volCenter.setPos (0, 0, 0);
  
  //! \todo if the ReflectiveBarrier does include the Receiver destination volume, then motion will continue forever
  motor.addVolumeSurface (volCenter, 1500.0, P1906MOL_VolSurface::ReflectiveBarrier);
  
  //! reset the motor's timer
  motor.initTime();
    
  //! volume surface must overlap with destination in order for the test to end
  motor.displayVolSurfaces();  
  /*
   * move randomly until destination reached
   */
  motor.setStartingPoint(startPt);
  motor.float2Destination(timePeriod);
  //printf ("(unitTest_MotorMovement) propagation time: %f\n", motor.getTime());
  mathematica.connectedPoints2Mma(motor.pos_history, "float2destination.mma");
  printf ("completed unitTest_ReflectiveBarrier\n");
  
  return true;
}

//! test the FluxMeter
bool P1906MOL_MicrotubulesField::unitTest_FluxMeter()
{
  P1906MOL_Motor motor;
  P1906MOL_MathematicaHelper mathematica;
  gsl_vector * startPt = gsl_vector_alloc (3);
  double timePeriod = 100;
  P1906MOL_Pos volCenter;

  printf ("beginning unitTest_ReflectiveBarrier\n");
  //! start at zero
  point (startPt, 0, 0, 0);
  volCenter.setPos (0, 0, 0);
  
  //! add a FluxMeter surface volume
  motor.addVolumeSurface (volCenter, 1500.0, P1906MOL_VolSurface::FluxMeter);
  
  //! reset the motor's timer
  motor.initTime();
    
  //! volume surface must overlap with destination in order for the test to end
  motor.displayVolSurfaces();  
  /*
   * move randomly until destination reached
   */
  motor.setStartingPoint(startPt);
  motor.float2Destination(timePeriod);
  //printf ("(unitTest_MotorMovement) propagation time: %f\n", motor.getTime());
  mathematica.connectedPoints2Mma(motor.pos_history, "float2destination.mma");
  printf ("completed unitTest_ReflectiveBarrier\n");
  
  return true;
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
  double d = P1906MOL_ExtendedField::distance (startPt, segment);
  printf ("distance: %f\n", d);
  printf ("completed unitTest_Distance\n");
  
  return true;
}

//! test finding a segment overlap
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

//! test finding all segment overlaps in a tube network
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
  persistenceVersusEntropy(&ts, persistenceLengths);
  printf ("completed unitTest_PersistenceLengthsVsEntropy\n");

  return true;
}

//! test creating a vector field with tubes2VectorField, creates vectorField.mma and vectorField.dat
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

//! test motor movement using Brownian motion until destination volume reached
bool P1906MOL_MicrotubulesField::unitTest_NoTubeMotion()
{
  P1906MOL_Motor motor;
  P1906MOL_MathematicaHelper mathematica;
  gsl_vector * startPt = gsl_vector_alloc (3);
  double timePeriod = 100;

  printf ("beginning unitTest_NoTubeMotion\n");
  //! reset the motor's timer
  motor.initTime();
  
  //! start at zero
  point (startPt, 0, 0, 0);
  
  //! volume surface must overlap with destination in order for the test to end
  motor.displayVolSurfaces();  
  /*
   * move randomly until destination reached
   */
  motor.setStartingPoint(startPt);
  motor.float2Destination(timePeriod);
  //printf ("(unitTest_MotorMovement) propagation time: %f\n", motor.getTime());
  mathematica.connectedPoints2Mma(motor.pos_history, "float2destination.mma");
  printf ("completed unitTest_NoTubeMotion\n");
  
  return true;
}

//! test motor movement to a tube and then walking along the tube
bool P1906MOL_MicrotubulesField::unitTest_MotorMovement(vector<P1906MOL_Pos> & pts)
{
  P1906MOL_Motor motor;
  P1906MOL_MathematicaHelper mathematica;
  gsl_vector * startPt = gsl_vector_alloc (3);
  //double timePeriod = 100;

  printf ("beginning unitTest_MotorMovement\n");
  //! reset the motor's timer
  motor.initTime();
  
  /*
   * move randomly until overlap with tube
   */
  motor.pos_history.clear();
  point (startPt, 
    gsl_matrix_get (tubeMatrix, 0, 0) + 30, //! start 10 nanometers away from the first tube segment
	gsl_matrix_get (tubeMatrix, 0, 1),
	gsl_matrix_get (tubeMatrix, 0, 2));
  motor.float2Tube(r, startPt, motor.pos_history, tubeMatrix, 0.1, motor.vsl);
  //printf ("completed float2Tube\n");
  //printf ("(unitTest_MotorMovement) float2Tube propagation time: %f\n", motor.getTime());
  //printf ("(unitTest_MotorMovement) float2Tube number of positions: %ld\n", motor.pos_history.size());
  mathematica.connectedPoints2Mma(motor.pos_history, "motion2tube.mma");
  //printf ("completed connectedPoints2Mma\n");
  
  //! start where the motor ended
  P1906MOL_Pos Pos;
  Pos = motor.pos_history.back();
  double x, y, z;
  Pos.getPos (&x, &y, &z);
  point(startPt, x, y, z);
  //printf ("(unitTest_MotorMovement) starting point for motor walk near tube\n");
  //displayPoint (startPt);

  /*
   * now walk along the tube
   */
  motor.pos_history.clear();
  motor.motorWalk(r, startPt, motor.pos_history, tubeMatrix, ts.segPerTube, motor.vsl);
  //printf ("(unitTest_MotorMovement) motorWalk propagation time: %f\n", motor.getTime());
  printf ("(unitTest_MotorMovement) motorWalk number of positions: %ld\n", motor.pos_history.size());
  mathematica.connectedPoints2Mma(motor.pos_history, "motion2end_of_tube.mma");
  //! append the motor history into pts
  pts.insert(pts.end(), motor.pos_history.begin(), motor.pos_history.end());
  printf ("completed unitTest_MotorMovement\n");
  
  return true;
}

//! test plotting to a Mathematica file using plot2mma, creates plottest.mma
bool P1906MOL_MicrotubulesField::unitTest_Plot2Mma(vector<P1906MOL_Pos> & pts)
{
  P1906MOL_MathematicaHelper mathematica;
  P1906MOL_Pos Pos;
  double x, y, z;
  
  printf ("beginning unitTest_Plot2Mma\n");
  //printf ("(unitTest_Plot2Mma) number of points: %ld\n", pts.size());
  if (pts.size() == 0) //! nothing to plot
  {
    printf ("(unitTest_Plot2Mma) nothing to plot\n");
    return false;
  }
  
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

//! test motor movement to destination using Brownian and microtubules if they exist
bool P1906MOL_MicrotubulesField::unitTest_MotorMove2Destination(vector<P1906MOL_Pos> & pts)
{
  P1906MOL_Motor motor;
  P1906MOL_MathematicaHelper mathematica;
  gsl_vector * startPt = gsl_vector_alloc (3);
  double timePeriod = 100;

  printf ("beginning unitTest_MotorMove2Destination\n");
  //! reset the motor's timer
  motor.initTime();
  
  //! start at zero
  point (startPt, 0, 0, 0);

  motor.pos_history.clear(); //! reset the position history
  motor.setStartingPoint(startPt);

  motor.move2Destination(tubeMatrix, ts.segPerTube, timePeriod, motor.pos_history);
  //printf ("(unitTest_MotorMove2Destination) propagation time: %f\n", motor.getTime());
  mathematica.connectedPoints2Mma(motor.pos_history, "motion2destination.mma");
  //! append the motor history into pts
  pts.insert(pts.end(), motor.pos_history.begin(), motor.pos_history.end());
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

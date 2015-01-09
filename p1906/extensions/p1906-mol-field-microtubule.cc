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
 * Author: Steve Bush - GE Global Research
 *                      bushsf@research.ge.com
 *                      http://www.amazon.com/author/stephenbush
 */

#include "ns3/log.h"

#include "p1906-mol-field-microtubule.h"
#include "ns3/p1906-field.h"

namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("P1906MOL_MicrotubulesField");

TypeId P1906MOL_MicrotubulesField::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::P1906MOL_MicrotubulesField")
    .SetParent<P1906Field> ();
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

  //! properties of the microtubule network
  tubeCharacteristcs_t ts;
  
  //! set the microtubule network properties
  setTubeVolume(&ts, 25);
  setTubeLength(&ts, 100);
  setTubeIntraAngle(&ts, 30);
  setTubeInterAngle(&ts, 10);
  setTubeDensity(&ts, 10);
  setTubePersistenceLength(&ts, 50);
  setTubeSegments(&ts, 10);
  
  //! simulated time
  struct simtime_t t;
  
  //! reset the simulated time
  initTime(&t);
  
  //! random number generation structures and initialization
  const gsl_rng_type * T;
  gsl_rng * r;
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_env_setup();
  
  //! display all the microtubule network properties
  displayTubeChars(&ts);
  
  //! create a given density of tubes of numSegments in given volume
  //! volume starts at 0, 0, 0 to volume^(1/3) in each dimension  
  gsl_matrix * tubeMatrix = gsl_matrix_alloc (ts.numTubes * ts.segPerTube, 6);
  genTubes(&ts, r, tubeMatrix);
  // printf ("tube matrix\n");
  // displayTube(tubeMatrix);
  tubes2Mma(tubeMatrix, ts.segPerTube, "tubes.mma");
  printf ("completed tube creation\n");
  
  //! find tube overlap points
  //! the maximum number of overlaps for n lines is n^2
  gsl_matrix * pts = gsl_matrix_alloc (10000, 3);
  size_t numPts = getAllOverlaps3D(tubeMatrix, pts);
  points2Mma(pts, numPts, "pfile.mma");
  printf ("completed finding tube overlaps\n");
  
  //! holds the vector field
  gsl_matrix * vf = gsl_matrix_alloc (ts.numTubes * ts.segPerTube, 6);
  
  //! print the vector field
  tubes2VectorField(tubeMatrix, vf);
  vectorFieldPlotMma(vf, "vectorField.mma");
  vectorFieldMeshMATLAB(vf, "vectorField.dat");
  printf ("completed plot of vector field\n");
  
  //! random motion for 100 sec
  //! start near the first tube
  gsl_vector * startPt = gsl_vector_alloc (3);
    
  numPts = freeFloat(r, startPt, pts, 100, 100, &t);
  connectedPoints2Mma(pts, numPts, "motion.mma");
  printf ("completed random motion\n");

  //! move randomly until overlap with tube
  point (startPt, 
    gsl_matrix_get (tubeMatrix, 0, 0) + 10,
	gsl_matrix_get (tubeMatrix, 0, 1),
	gsl_matrix_get (tubeMatrix, 0, 2));
  float2Tube(r, startPt, pts, &numPts, tubeMatrix, 0.1, &t);
  printf ("completed float to tube\n");
	
  connectedPoints2Mma(pts, numPts, "motion2tube.mma");
  printf ("completed connectedPoints2Mma\n");
  printf ("numPts: %ld\n", numPts);
  printf ("pts->size1: %ld\n", pts->size1);
  point(startPt,
    gsl_matrix_get (pts, numPts - 1, 0),
	gsl_matrix_get (pts, numPts - 1, 1),
	gsl_matrix_get (pts, numPts - 1, 2));
  printf ("point after random walk\n");
  displayPoint (startPt);
  printf ("completed displayPoint\n");

  //! now walk along the tube
  motorWalk(r, startPt, pts, &numPts, tubeMatrix, ts.segPerTube, &t);
  connectedPoints2Mma(pts, numPts, "tubewalk.mma");
  printf ("completed walk along tube\n");
  
  gsl_vector * segment = gsl_vector_alloc (6);
  gsl_vector * pt1 = gsl_vector_alloc (3);
  gsl_vector * pt2 = gsl_vector_alloc (3);
  point (startPt, 0, 0, 0);
  point (pt1, -1, -1, -1);
  point (pt2, 2, 2, 2);
  line (segment, pt1, pt2);
  double d = distance (startPt, segment);
  printf ("distance: %f\n", d);
  printf ("completed distance test\n");
  
  //! see http://www.gnu.org/software/gsl/manual/html_node/Level-3-GSL-BLAS-Interface.html#Level-3-GSL-BLAS-Interface
  /* gsl_blas_dgemm (CblasNoTrans, CblasTrans,
                  1.0, segStartX, segAngle,
                  0.0, segEndX); */
  /* result of element-wise multiplication is in first arg */				  
  /* gsl_matrix_mul_elements (segEndX, segAngle); */
  
  //! test segment overlap
  gsl_vector * segment3D = gsl_vector_alloc (6);
  gsl_matrix * pts3D = gsl_matrix_alloc (1, 3);
  gsl_matrix * tubeMatrix3D = gsl_matrix_alloc (1, 6);
  gsl_matrix_set_zero (pts3D);
  gsl_vector * pt3 = gsl_vector_alloc (3);
  gsl_vector * pt4 = gsl_vector_alloc (3);
  gsl_vector * tubeSegments = gsl_vector_alloc (1);
  
  point (pt1, 0, 0, 0);
  point (pt2, 5, 5, 0);
  point (pt3, 5, 0, 0);
  point (pt4, 0, 5, 0);
  line (segment3D, pt1, pt2);
  line (tubeMatrix3D, 0, pt3, pt4);
  getOverlap3D(segment3D, tubeMatrix3D, pts3D, tubeSegments);
  if (isPointOverlap(pt1, segment3D))
    printf ("point overlaps\n");
  printf ("completed overlap test\n");
	
  //! test plotting points
  gsl_matrix * vals = gsl_matrix_alloc (numPts, 2);
  for (size_t i = 0; i < numPts; i++)
    for (size_t j = 0; j < 2; j++)
      gsl_matrix_set (vals, i, j, gsl_matrix_get (pts, i, j));
  plot2Mma(vals, "plottest.mma", "x value", "y value");
  printf ("completed plot test\n");
  
  gsl_vector * persistenceLengths = gsl_vector_alloc (10);
  for (size_t i = 0; i < 10; i++)
    gsl_vector_set (persistenceLengths, i, i * 100);
  perstenceVersusEntropy(&ts, r, persistenceLengths);
  
  //! close down field activity
  gsl_rng_free (r);

  NS_LOG_FUNCTION (this);
  NS_LOG_FUNCTION (this << "Created MOL Field Component");
}

//! compute and return the active network programmability metric at pt
//! \todo use vector field and divergence to compute this value
void P1906MOL_MicrotubulesField::activeNetworkProgrammability(gsl_matrix * vf, gsl_vector * pt)
{
  //! sum the vector input and output to point pt
  //! may add a radius argument to define a sphere for determining flux through the surface of the sphere
  //! activeNetworkProgrammability is fundamentally a change to the vector field
}

//! display the vector field in Mathematica format for VectorPlot3D in file fname
void P1906MOL_MicrotubulesField::vectorFieldPlotMma(gsl_matrix * vf, const char *fname)
{
  FILE * pFile;

  pFile = fopen (fname,"w");
  
  fprintf (pFile, "ListVectorPlot3D[{");
  for (size_t i = 0; i < vf->size1; i++)
  {
    fprintf (pFile, "{{%f, %f, %f}, {%f, %f, %f}}", 
	  gsl_matrix_get (vf, i, 0),
	  gsl_matrix_get (vf, i, 1),
	  gsl_matrix_get (vf, i, 2),
	  gsl_matrix_get (vf, i, 3),
	  gsl_matrix_get (vf, i, 4),
	  gsl_matrix_get (vf, i, 5));
    if (i < (vf->size1 - 1)) fprintf (pFile, ", ");
  }
  //! fprintf (pFile, "}, ");
  //! fprintf (pFile, ", PlotStyle -> {Dashed, Thick, Red}");
  fprintf (pFile, "}]\n");
  //! option to print vertex labels
  //! fprintf (pFile, "}, VertexLabeling -> True]\n");
  
  fclose(pFile);
}

//! write the vector field in Mathematica format using regular spacing between samples in file fname
void P1906MOL_MicrotubulesField::vectorFieldMeshMma(gsl_matrix * vf, const char *fname)
{
  FILE * pFile;

  pFile = fopen (fname,"w");
  
  fprintf (pFile, "ListVectorPlot3D[{");
  for (size_t i = 0; i < vf->size1; i++)
  {
    fprintf (pFile, "{{%f, %f, %f}, {%f, %f, %f}}", 
	  gsl_matrix_get (vf, i, 0),
	  gsl_matrix_get (vf, i, 1),
	  gsl_matrix_get (vf, i, 2),
	  gsl_matrix_get (vf, i, 3),
	  gsl_matrix_get (vf, i, 4),
	  gsl_matrix_get (vf, i, 5));
    if (i < (vf->size1 - 1)) fprintf (pFile, ", ");
  }
  //! fprintf (pFile, "}, ");
  //! fprintf (pFile, ", PlotStyle -> {Dashed, Thick, Red}");
  fprintf (pFile, "}]\n");
  //! option to print vertex labels
  //! fprintf (pFile, "}, VertexLabeling -> True]\n");
  
  fclose(pFile);
}

//! return vector field vf based upon tubes in tubeMatrix
void P1906MOL_MicrotubulesField::tubes2VectorField(gsl_matrix * tubeMatrix, gsl_matrix * vf)
{
  gsl_matrix * v = gsl_matrix_alloc (tubeMatrix->size1, 3);
  gsl_matrix * pt = gsl_matrix_alloc (tubeMatrix->size1, 3);
  
  for (size_t i = 0; i < tubeMatrix->size1; i++)
    for (size_t j = 0; j < 3; j++)
    {
	  gsl_matrix_set (pt, i, j, gsl_matrix_get (tubeMatrix, i, j));
      gsl_matrix_set (v, i, j, gsl_matrix_get (tubeMatrix, i, j + 3));
    }
	
  // printf ("(tubes2VectorField) set pt and v\n");
  gsl_matrix_sub (v, pt);
  // printf ("(tubes2VectorField) subtracted pt from v\n");
  
  for (size_t i = 0; i < tubeMatrix->size1; i++)
    for (size_t j = 0; j < 3; j++)
    {
	  gsl_matrix_set (vf, i, j, gsl_matrix_get (pt, i, j));
      gsl_matrix_set (vf, i, j + 3, gsl_matrix_get (v, i, j));
    }  
  // printf ("(tubes2VectorField) set vf\n");
}

//! return simulation time
double P1906MOL_MicrotubulesField::propagationDelay(struct simtime_t * t)
{
  return t->time;
}

//! initialize simulation time
void P1906MOL_MicrotubulesField::initTime(struct simtime_t * t)
{
  t->time = 0;
}

//! update simulation time
void P1906MOL_MicrotubulesField::updateTime(struct simtime_t * t, double event_time)
{
  t->time += event_time;
}

//! print tube characteristics to standard output
void P1906MOL_MicrotubulesField::displayTubeChars(struct tubeCharacteristcs_t * ts)
{
  printf ("volume = %f\n", ts->volume);
  printf ("mean_tube_length = %f\n", ts->mean_tube_length);
  printf ("mean_intra_tube_angle = %f\n", ts->mean_intra_tube_angle);
  printf ("mean_inter_tube_angle = %f\n", ts->mean_inter_tube_angle);
  printf ("mean_tube_density = %f\n", ts->mean_tube_density);
  printf ("segLength = %f\n", ts->segLength);
  printf ("numSegments = %ld\n", ts->numSegments);
}

//! set the volume in which tube centers exist
void P1906MOL_MicrotubulesField::setTubeVolume(struct tubeCharacteristcs_t * ts, double volume)
{
  ts->volume = volume;
}

//! set the mean tube length
void P1906MOL_MicrotubulesField::setTubeLength(struct tubeCharacteristcs_t * ts, double mean_tube_length)
{
  ts->mean_tube_length = mean_tube_length;
  ts->segLength = ts->mean_tube_length / 5;
}

//! set the mean angle between segments of the tube
void P1906MOL_MicrotubulesField::setTubeIntraAngle(struct tubeCharacteristcs_t * ts, double mean_intra_tube_angle)
{
  ts->mean_intra_tube_angle = mean_intra_tube_angle;
}

//! set the mean angle between tubes
void P1906MOL_MicrotubulesField::setTubeInterAngle(struct tubeCharacteristcs_t * ts, double mean_inter_tube_angle)
{
  ts->mean_inter_tube_angle = mean_inter_tube_angle;
}

//! this is really the segment density
void P1906MOL_MicrotubulesField::setTubeDensity(struct tubeCharacteristcs_t * ts, double mean_tube_density)
{
  ts->mean_tube_density = mean_tube_density;
  ts->numSegments = ts->mean_tube_density * ts->volume;
}

//! set the persistence length of each tube
void P1906MOL_MicrotubulesField::setTubePersistenceLength(struct tubeCharacteristcs_t * ts, double persistenceLength)
{
  ts->persistenceLength = persistenceLength;
}

//! set the number of segments per tube
void P1906MOL_MicrotubulesField::setTubeSegments(struct tubeCharacteristcs_t * ts, size_t segPerTube)
{
  ts->segPerTube = segPerTube;
  ts->numTubes = floor(ts->numSegments / ts->segPerTube);
}

//! plot persistence length versus structural entropy
void P1906MOL_MicrotubulesField::perstenceVersusEntropy(struct tubeCharacteristcs_t * ts, gsl_rng * r, gsl_vector * persistenceLengths)
{
  char plot_filename[256];
  gsl_matrix * tubeMatrix = gsl_matrix_alloc (ts->numTubes * ts->segPerTube, 6);
  gsl_matrix * pve = gsl_matrix_alloc (persistenceLengths->size, 2);
  
  for (size_t i = 0; i < persistenceLengths->size; i++)
  {
    setTubePersistenceLength (ts, gsl_vector_get (persistenceLengths, i));
    genTubes(ts, r, tubeMatrix);
	
	sprintf (plot_filename, "tubes_%ld.mma", i);
    tubes2Mma(tubeMatrix, ts->segPerTube, plot_filename);
	gsl_matrix_set (pve, i, 0, gsl_vector_get (persistenceLengths,i));
	gsl_matrix_set (pve, i, 1, ts->se);
  }
  
  //! plot the results
  plot2Mma(pve, "persistenceVersusEntropy.mma", "persistence length", "structural entropy");
}

//! return the index of the nearest tube in tubeMatrix within a given radius from pt otherwise return -1
int P1906MOL_MicrotubulesField::findNearestTube(gsl_vector *pt, gsl_matrix *tubeMatrix, double radius)
{
  double shortestDistance = GSL_POSINF;
  double d = 0;
  size_t closestSegment = -1;
  gsl_vector *segment = gsl_vector_alloc (6);
  
  for (size_t i = 0; i < tubeMatrix->size1; i++)
  {
    line(segment, tubeMatrix, i);
	d = distance(pt, segment);
	if ((d < shortestDistance) && (d <= radius))
	{
	  shortestDistance = d;
	  closestSegment = i;
	}	
  }
  
  // printf ("shortestDistance: %f\n", shortestDistance);
  return closestSegment;
}

//! return the distance between two point pt1 and point pt2
double P1906MOL_MicrotubulesField::distanceP(gsl_vector *pt1, gsl_vector *pt2)
{
  gsl_vector * d = gsl_vector_alloc (3);
  
  //! copy pt1 to d
  gsl_vector_memcpy (d, pt1);
  //! subtract pt2 from pt1 in d
  gsl_vector_sub (d, pt2);
  //! return the Euclidean distance
  return gsl_blas_dnrm2 (d);
}

//! return the shortest distance between point pt and the line or point segment_or_point
//! if segment_or_point is a vector of length 3, then it is a point
//! if segment_or_point is a vector of length 6, then it is a line segment described by two end points
double P1906MOL_MicrotubulesField::distance(gsl_vector *pt, gsl_vector *segment_or_point)
{
  gsl_vector * td = gsl_vector_alloc (3);
  gsl_vector * segment = gsl_vector_alloc (6);
  gsl_vector * pt1 = gsl_vector_alloc (3);
  gsl_vector * pt2 = gsl_vector_alloc (3);
  gsl_vector * res = gsl_vector_alloc (3);
  gsl_vector * res1 = gsl_vector_alloc (3);
  gsl_vector * res2 = gsl_vector_alloc (3);
  double d;
	  
  switch (segment_or_point->size)
  {
    case 3:
      gsl_vector_memcpy (pt1, pt);
      gsl_vector_memcpy (pt2, segment_or_point);
      gsl_vector_memcpy (td, pt1);
      gsl_vector_sub (td, pt2);
      return gsl_blas_dnrm2 (td);
    case 6:
   	  gsl_vector_memcpy (segment, segment_or_point);
	  
      //! distance = |pt - pt1| x |pt - pt2| / |pt2 - pt1| where x is the cross product
      //! \f$|\vec{pt} - \vec{pt1}| \times |\vec{pt} - \vec{pt2}| / |\vec{pt2} - \vec{pt1}|\f$
      point (pt1,
        gsl_vector_get (segment, 0),
	    gsl_vector_get (segment, 1),
	    gsl_vector_get (segment, 2));
      point (pt2,
        gsl_vector_get (segment, 3),
	    gsl_vector_get (segment, 4),
	    gsl_vector_get (segment, 5));
	
      gsl_vector_memcpy (res1, pt);
      gsl_vector_sub (res1, pt1);
  
      // printf ("res1\n");
      // displayPoint (res1);
  
      gsl_vector_memcpy (res2, pt);
      gsl_vector_sub (res2, pt2);
  
      // printf ("res2\n");
      // displayPoint (res2);
  
      cross_product (res1, res2, res);
      // printf ("res\n");
      // displayPoint (res);
  
      gsl_vector_memcpy (res1, pt2);
      gsl_vector_sub (res1, pt1);
      // printf ("res1\n");
      // displayPoint (res1);
  
      d = gsl_blas_dasum (res) / gsl_blas_dasum (res1);
      // printf ("distance: %g\n", d);
  
      return d;
      /* gsl_blas_sdsdot (0, x, y, result); */
    default:
      printf ("(distance) invalid argument to distance\n");
	  return -1;
  }
}

//! the vector cross product is normal to the input vectors and has magnitude equivalent to the volume of a parallelogram formed by the input vectors
//! u and v are the input vectors and product is the output vector
void P1906MOL_MicrotubulesField::cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *product)
{
        double p1 = gsl_vector_get(u, 1)*gsl_vector_get(v, 2)
                - gsl_vector_get(u, 2)*gsl_vector_get(v, 1);
 
        double p2 = gsl_vector_get(u, 2)*gsl_vector_get(v, 0)
                - gsl_vector_get(u, 0)*gsl_vector_get(v, 2);
 
        double p3 = gsl_vector_get(u, 0)*gsl_vector_get(v, 1)
                - gsl_vector_get(u, 1)*gsl_vector_get(v, 0);
 
        gsl_vector_set(product, 0, p1);
        gsl_vector_set(product, 1, p2);
        gsl_vector_set(product, 2, p3);
}

//! \todo verify that motorWalk is working properly
void P1906MOL_MicrotubulesField::motorWalk(gsl_rng * r, gsl_vector * startPt, gsl_matrix *pts, size_t * numPts, gsl_matrix * tubeMatrix, size_t segPerTube, struct simtime_t * t)
{
  /** See "Movements of Molecular Motors," Reinhard Lipowsky
	movement speed: ~1 um / sec
	bound time ~2 sec
	assumes startPt is on a tube in tubeMatrix
  */
  gsl_vector * segment = gsl_vector_alloc(6);
  gsl_vector * pt1 = gsl_vector_alloc(3);
  gsl_vector * pt2 = gsl_vector_alloc(3);
  int p = 0;
  double radius = 15;
  double movementRate = 1000; //! nm / sec
  //! double bindingTime = 2; sec

  //! find the tube the motor is starting on
  int seg = findNearestTube(startPt, tubeMatrix, radius); //! \todo set tube radius (thickness) globally
  
  //! store the starting point
  for (int j = 0; j < 3; j++)
	  gsl_matrix_set(pts, p, j, gsl_vector_get(startPt, j));
  p++;
  
  //! walk along tube for distance determined by bound time
  //! segments are sequential in tubeMatrix of length segPerTube
  int segOfTube = seg % segPerTube; //! the segment within the tube
  int segToGo = segPerTube - segOfTube; //! segments until end of tube
  
  // printf ("seg: %d segOfTube: %d segToGo: %d\n", seg, segOfTube, segToGo);
  for (int i = seg; i < (seg + segToGo); i++)
  {
    //! walk to end of segment
    line(segment, tubeMatrix, i);
	for (int j = 0; j < 3; j++)
	  gsl_matrix_set(pts, p, j, gsl_vector_get(segment, j + 3));
	point (pt1, 
	  gsl_matrix_get (pts, p - 1, 0),
	  gsl_matrix_get (pts, p - 1, 1),
	  gsl_matrix_get (pts, p - 1, 2));
    point (pt2, 
	  gsl_matrix_get (pts, p, 0),
	  gsl_matrix_get (pts, p, 1),
	  gsl_matrix_get (pts, p, 2));
    updateTime(t, distanceP(pt1, pt2) / movementRate);
	p++;
  }
  
  *numPts = p;
}

//! return true if point pt lies on segment, false otherwise
bool P1906MOL_MicrotubulesField::isPointOverlap(gsl_vector * pt, gsl_vector * segment)
{
	bool overlap = false;
	/**
	  Determine whether the points A(2,4,2) B(3,7,-3) and C(1,3,3) lie on a straight line.
	  Start by determining a direction vector AB = [a, b, c] for the line connecting points A and B. 
	  This vector gives the slopes of the three parametric equations describing the line between the points. 
	  AB = [\Delta x, \Delta y, \Delta z] = [3 - 2, 7 - 4, -3 - 2] = [1, 3, -5]. 
	  Three parametric equations for a line are: 
        x = x1 + at 
        y = y1 + bt 
        z = z1 + ct, 
      So if we let A = (x1, y1, z1) 
        x = 2 + t 
        y = 4 + 3t 
        z = 2 - 5t 
      These three parametric equations describe the line connecting A and B. 
	  To see if C is on this line, we can substitute its x coordinate into the x parametric equation to find t: 
        1 = 2 + t, 
      so t = -1. Now substitute this t into the other two equations: 
        y = 4 + 3(-1) = 1; 
        z = 2 - 5(-1) = 7. 
      So (1, 1, 7) is a point on the line, but C = (1, 3, 3) is not.
    */
	double px = gsl_vector_get (pt, 0);
	double py = gsl_vector_get (pt, 1);
	double pz = gsl_vector_get (pt, 2);
	
	double sx1 = gsl_vector_get (segment, 0);
	double sy1 = gsl_vector_get (segment, 1);
	double sz1 = gsl_vector_get (segment, 2);
	
	double sx2 = gsl_vector_get (segment, 0);
	double sy2 = gsl_vector_get (segment, 1);
	double sz2 = gsl_vector_get (segment, 2);
	
	//! direction vector
	gsl_vector * dV = gsl_vector_alloc (3);
	
	//! slope of each dimension
	gsl_vector_set (dV, 0, sx2 - sx1);
	gsl_vector_set (dV, 1, sy2 - sy1);
	gsl_vector_set (dV, 2, sz2 - sz1);
	
	if ( ( ((px - sx1) / gsl_vector_get (dV, 0)) == ((py - sy1) / gsl_vector_get (dV, 1)) ) &&  
		 ( ((py - sy1) / gsl_vector_get (dV, 1)) == ((pz - sz1) / gsl_vector_get (dV, 2)) ) )
	  overlap = true;
	
	return overlap;
}

//! motor is unbound from tube and floats via Brownian motion until it makes contact with a tube
//!   pts - locations of the motor during its random walk comprised of npts
//!   startPt - where the motor began its random walk
//!   timePeriod - length of each step of the walk
//!   returns the index of the contact segment in tubeMatrix
int P1906MOL_MicrotubulesField::float2Tube(gsl_rng * r, gsl_vector * startPt, gsl_matrix * pts, size_t * npts, gsl_matrix * tubeMatrix, double timePeriod, struct simtime_t * t)
{
  gsl_vector * currentPos = gsl_vector_alloc (3);
  gsl_vector * newPos = gsl_vector_alloc (3);
  int numPts = 0; //! total number of points traversed
  double timeout = 200; //! end after this many steps
  int ts; //! nearest tube segment
  double radius = 15;
  
  //! begin at the starting point
  point (currentPos, 
    gsl_vector_get (startPt, 0), 
	gsl_vector_get (startPt, 1), 
	gsl_vector_get (startPt, 2));

  //! float to the nearest tube within a given radius
  for (int i = 0; i < timeout; i++)
  {
    gsl_matrix_set (pts, i, 0, gsl_vector_get (currentPos, 0));
	gsl_matrix_set (pts, i, 1, gsl_vector_get (currentPos, 1));
	gsl_matrix_set (pts, i, 2, gsl_vector_get (currentPos, 2));
	ts = findNearestTube(currentPos, tubeMatrix, radius);
	numPts++; //! consider starting position the first point
	if (ts != -1)
	{
	  printf ("motor contact with segment: %d\n", ts);
	  *npts = numPts;
	  return ts;
	}
	brownianMotion(r, currentPos, newPos, timePeriod);
	updateTime(t, timePeriod);
    gsl_vector_set (currentPos, 0, gsl_vector_get (newPos, 0));
	gsl_vector_set (currentPos, 1, gsl_vector_get (newPos, 1));
	gsl_vector_set (currentPos, 2, gsl_vector_get (newPos, 2));
  }
  
  *npts = numPts;
  return ts;
}

//! implements a motor floating via Brownian motion for time steps with step lengths of timePeriod
int P1906MOL_MicrotubulesField::freeFloat(gsl_rng * r, gsl_vector * startPt, gsl_matrix * pts, int time, double timePeriod, struct simtime_t * t)
{
  gsl_vector * currentPos = gsl_vector_alloc (3);
  gsl_vector * newPos = gsl_vector_alloc (3);
  point (currentPos, gsl_vector_get (currentPos, 0), gsl_vector_get (currentPos, 1), gsl_vector_get (currentPos, 2));
  int numPts = 0;
  
  for (int i = 0; i < time; i++)
  {
    gsl_matrix_set (pts, i, 0, gsl_vector_get (currentPos, 0));
	gsl_matrix_set (pts, i, 1, gsl_vector_get (currentPos, 1));
	gsl_matrix_set (pts, i, 2, gsl_vector_get (currentPos, 2));
	brownianMotion(r, currentPos, newPos, timePeriod);
	updateTime(t, timePeriod);
    gsl_vector_set (currentPos, 0, gsl_vector_get (newPos, 0));
	gsl_vector_set (currentPos, 1, gsl_vector_get (newPos, 1));
	gsl_vector_set (currentPos, 2, gsl_vector_get (newPos, 2));
    numPts++;
  }
  
  return numPts;
}

//! generate a set of tubes comprised of a total of numSegments in volume with segPerTube segments of segLength 
//! and persistenceLength and return in tubeMatrix
void P1906MOL_MicrotubulesField::genTubes(struct tubeCharacteristcs_t * ts, gsl_rng * r, gsl_matrix * tubeMatrix)
{
  //! \todo get actual tube graph properties from biologist
  gsl_vector * startPt = gsl_vector_alloc (3);
  //! hold the values for a tube comprised of many segments: x_start y_start x_start x_end y_end z_end
  gsl_matrix * segMatrix = gsl_matrix_alloc (ts->segPerTube, 6);
  double total_structural_entropy = 0;
  
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
    genTube(ts, r, segMatrix, startPt);
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

//! set point pt with x, y, and z coordinates
void P1906MOL_MicrotubulesField::point(gsl_vector * pt, double x, double y, double z)
{
  gsl_vector_set (pt, 0, x);
  gsl_vector_set (pt, 1, y);
  gsl_vector_set (pt, 2, z);
}

//! set line with end points pt1 and pt2
void P1906MOL_MicrotubulesField::line(gsl_vector * line, gsl_vector * pt1, gsl_vector * pt2)
{
  for (int i = 0; i < 3; i++)
    gsl_vector_set (line, i, gsl_vector_get (pt1, i));

  for (int i = 0; i < 3; i++)
    gsl_vector_set (line, i + 3, gsl_vector_get (pt2, i));
}

//! return the segment in location mp from tubeMatrix
void P1906MOL_MicrotubulesField::line(gsl_vector * segment, gsl_matrix * tubeMatrix, int mp)
{
  for (size_t i = 0; i < 6; i++)
    gsl_vector_set (segment, i, gsl_matrix_get(tubeMatrix, mp, i));
}

//! set the segment at location mp in the line matrix with pt1 and pt2
void P1906MOL_MicrotubulesField::line(gsl_matrix * line, int mp, gsl_vector * pt1, gsl_vector * pt2)
{
  for (int i = 0; i < 3; i++)
    gsl_matrix_set (line, mp, i, gsl_vector_get (pt1, i));

  for (int i = 0; i < 3; i++)
    gsl_matrix_set (line, mp, i + 3, gsl_vector_get (pt2, i));
}

//! print the tube segments in segMatrix
void P1906MOL_MicrotubulesField::displayTube(gsl_matrix *segMatrix)
{
  for (size_t i = 0; i < segMatrix->size1; i++)
  {
    for (size_t j = 0; j < segMatrix->size2; j++)
	{
	  printf ("segMatrix(%ld,%ld) = %g\t", i, j, gsl_matrix_get (segMatrix, i, j));
    }
	printf ("\n");
  }
}

//! print the points in pts
void P1906MOL_MicrotubulesField::displayPoints(gsl_matrix *pts)
{
  size_t numPts = pts->size1;
  
  for (size_t i = 0; i < numPts; i++)
    printf ("Point: %g %g %g\n", 
	  gsl_matrix_get(pts, i, 0), 
	  gsl_matrix_get(pts, i, 1),
	  gsl_matrix_get(pts, i, 2));
}

//! print the points (vertices) in pts in Mathematica format into file fname and include edges between the vertices
void P1906MOL_MicrotubulesField::connectedPoints2Mma(gsl_matrix *pts, size_t numPts, const char* fname)
{
  FILE * pFile;

  pFile = fopen (fname,"w");
  size_t pt = 1;
  
  fprintf (pFile, "GraphPlot3D[{");
  for (size_t i = 0; i < numPts; i++)
  {
    fprintf (pFile, "%ld -> ", pt);
	pt++;
	fprintf (pFile, "%ld", pt);
    if (i < (numPts - 1)) fprintf (pFile, ", ");
  }
  fprintf (pFile, "}, ");
  pt = 1;
  fprintf (pFile, "VertexCoordinateRules ->{");
  for (size_t i = 0; i < numPts; i++)
  {
    fprintf (pFile, "%ld -> {%f, %f, %f}",
      pt,
	  gsl_matrix_get(pts, i, 0),
      gsl_matrix_get(pts, i, 1),
	  gsl_matrix_get(pts, i, 2));
    pt++;
	if (i < (numPts - 1)) fprintf (pFile, ", ");
  }
  fprintf (pFile, "}");
  fprintf (pFile, ", PlotStyle -> {Dashed, Thick, Red}");
  fprintf (pFile, "]\n");
  //! option to print vertex labels
  //! fprintf (pFile, "}, VertexLabeling -> True]\n");
  
  fclose(pFile);
}

//! print the first numPts points pts in Mathematica format in file fname
void P1906MOL_MicrotubulesField::points2Mma(gsl_matrix *pts, size_t numPts, const char* fname)
{
  FILE * pFile;

  pFile = fopen (fname,"w");
  
  fprintf (pFile, "Graphics3D[{PointSize[Large], Blue, ");
  for (size_t i = 0; i < numPts; i++)
  {
    fprintf (pFile, "Point[{%f, %f, %f}]", 
	  gsl_matrix_get(pts, i, 0), 
	  gsl_matrix_get(pts, i, 1),
	  gsl_matrix_get(pts, i, 2));
	if (i < numPts - 1) fprintf (pFile, ", ");
  }
  fprintf (pFile, "}]\n");
  
  fclose(pFile);
}

//! print a plot of x,y values in vals in Mathematica format into file fname
void P1906MOL_MicrotubulesField::plot2Mma(gsl_matrix *vals, const char* fname, const char* xlabel, const char* ylabel)
{
  FILE * pFile;
  size_t numVals = vals->size1;
  
  //! ListLinePlot[Quantity[{0, 3, 6, 8, 10, 11, 11, 16, 20, 22}, "Centimeters"], AxesLabel -> Automatic]

  pFile = fopen (fname,"w");
  
  fprintf (pFile, "ListLinePlot[{");
  for (size_t i = 0; i < numVals; i++)
  {
    fprintf (pFile, "{%f, %f}", 
	  gsl_matrix_get(vals, i, 0), 
	  gsl_matrix_get(vals, i, 1));
    if (i < (numVals - 1)) fprintf (pFile, ", ");
  }
  fprintf (pFile, "}");
  fprintf (pFile, ", AxesLabel -> {\"%s\", \"%s\"}, GridLines -> Automatic", xlabel, ylabel);
  fprintf (pFile, "]\n");
  
  fclose(pFile);
}

//! print the point pt
void P1906MOL_MicrotubulesField::displayPoint(gsl_vector *pt)
{
  printf ("Point: %g %g %g\n", 
	gsl_vector_get(pt, 0), 
	gsl_vector_get(pt, 1),
	gsl_vector_get(pt, 2));
}

//! print the first numPts points pts
void P1906MOL_MicrotubulesField::displayPoints(gsl_matrix *pts, size_t numPts)
{
  for (size_t i = 0; i < numPts; i++)
    printf ("Point: %g %g %g\n", 
	  gsl_matrix_get(pts, i, 0), 
	  gsl_matrix_get(pts, i, 1),
	  gsl_matrix_get(pts, i, 2));
}

//! print the position pts
void P1906MOL_MicrotubulesField::displayPos(gsl_vector *pt)
{
  printf ("Position: %g %g %g\n", 
    gsl_vector_get(pt, 0), 
	gsl_vector_get(pt, 1),
	gsl_vector_get(pt, 2));
}

//! a unit test of the findClosestPoint function
bool P1906MOL_MicrotubulesField::unitTestfindClosestPoint()
{
  bool pass = false;
  gsl_vector * pt = gsl_vector_alloc (3);
  gsl_vector * result = gsl_vector_alloc (6);
  gsl_matrix * vf = gsl_matrix_alloc (3, 6);
    
  gsl_matrix_set (vf, 0, 0, 0);
  gsl_matrix_set (vf, 0, 1, 0);
  gsl_matrix_set (vf, 0, 2, 0);
  gsl_matrix_set (vf, 0, 3, 0);
  gsl_matrix_set (vf, 0, 4, 0);
  gsl_matrix_set (vf, 0, 5, 0);

  gsl_matrix_set (vf, 1, 0, 1);
  gsl_matrix_set (vf, 1, 1, 1);
  gsl_matrix_set (vf, 1, 2, 1);
  gsl_matrix_set (vf, 1, 3, 1);
  gsl_matrix_set (vf, 1, 4, 1);
  gsl_matrix_set (vf, 1, 5, 1);

  gsl_matrix_set (vf, 2, 0, 5);
  gsl_matrix_set (vf, 2, 1, 5);
  gsl_matrix_set (vf, 2, 2, 5);
  gsl_matrix_set (vf, 2, 3, 5);
  gsl_matrix_set (vf, 2, 4, 5);
  gsl_matrix_set (vf, 2, 5, 5);
  
  point (pt, 1, 1, 1);

  findClosestPoint(pt, vf, result);
  
  printf ("result %f %f %f %f %f %f\n", 
    gsl_vector_get (result, 0),
	gsl_vector_get (result, 1),
	gsl_vector_get (result, 2),
	gsl_vector_get (result, 3),
	gsl_vector_get (result, 4),
	gsl_vector_get (result, 5));
  
  if (gsl_vector_get (result, 3) == 0 &&
      gsl_vector_get (result, 4) == 0 &&
      gsl_vector_get (result, 5) == 0)
    pass = true;
  else
    pass = false;
	
  return pass;
}

//! a unit test for getOverlap
bool P1906MOL_MicrotubulesField::unitTestgetOverlap()
{
  //! \todo create unit tests
  bool passTests = false;
  gsl_vector * segment = gsl_vector_alloc (6);
  gsl_matrix * tubeMatrix = gsl_matrix_alloc (1, 6);
  gsl_matrix * pts = gsl_matrix_alloc (1, 6);
  gsl_matrix_set_zero (pts);
  gsl_vector * pt1 = gsl_vector_alloc (3);
  gsl_vector * pt2 = gsl_vector_alloc (3);
  gsl_vector * pt3 = gsl_vector_alloc (3);
  gsl_vector * pt4 = gsl_vector_alloc (3);
  gsl_vector *tubeSegments = gsl_vector_alloc (1);
  
  //! Test 1: Simple cross
  
  point (pt1, 25, 0, 0);
  point (pt2, 0, 25, 0);
  point (pt3, 0, 25, 0);
  point (pt4, 25, 0, 0);
  
  line (segment, pt1, pt2);
  line (tubeMatrix, 0, pt3, pt4);
 
  int numPts = getOverlap3D(segment, tubeMatrix, pts, tubeSegments);
  printf ("numPts: %d\n", numPts);
  displayPoints (pts);
  
  //! Test 2: Moving cross
  for (int theta = 0; theta < 2 * M_PI; theta += M_PI / 6);
  //! convert theta to x, y values
  
  //! Test 2: Rotating segment
  
  return passTests;  
}

//! return all overlapping points in tubeMatrix in the list of points pts
size_t P1906MOL_MicrotubulesField::getAllOverlaps3D(gsl_matrix *tubeMatrix, gsl_matrix *pts)
{
  size_t numSegments = tubeMatrix->size1;
  gsl_vector *segment = gsl_vector_alloc (6);
  gsl_matrix *tmpPts = gsl_matrix_alloc (numSegments * numSegments, 3);
  size_t numPts = 0, totPts = 0;
  size_t pp = 0;
  gsl_vector *tubeSegments = gsl_vector_alloc (numSegments * numSegments);
  
  for (size_t i = 0; i < numSegments; i++)
  {
    //! move segment from tubes to segment
	for (size_t k = 0; k < 6; k++)
	  gsl_vector_set (segment, k, gsl_matrix_get (tubeMatrix, i, k));
	//! check for overlap
    numPts = getOverlap3D(segment, tubeMatrix, tmpPts, tubeSegments);
	//! store overlapping points
	for (size_t k = 0; k < numPts; k++)
	{
	  for (size_t j = 0; j < 3; j++)
	  {
	    // printf ("tmpPts(%ld, %ld) = %g\n", k, j, gsl_matrix_get (tmpPts, k, j));
	    gsl_matrix_set (pts, pp, j, gsl_matrix_get (tmpPts, k, j));
        // printf ("pts(%ld, %ld) = %g\n", pp, j, gsl_matrix_get (pts, pp, j));
	  }
	  pp++;
	}
    totPts += numPts;
  }
  
  return totPts;
}

//! return the number of overlapping points in pts and the index of the tubeMatrix segments overlapped in tubeSegments
int P1906MOL_MicrotubulesField::getOverlap3D(gsl_vector *segment, gsl_matrix *tubeMatrix, gsl_matrix *pts, gsl_vector *tubeSegments)
{
  /** 
    all points defined by x, y, z values
	line A -> B: (a1, a2, a3) -> (b1, b2, b3)
	line C -> D: (c1, c2, c3) -> (d1, d2, d3)
  
    solve: 
    a1 + t * (b1 - a1) == c1 + s * (d1 - c1)
	a2 + t * (b2 - a2) == c2 + s * (d2 - c2)
	a3 + t * (b3 - a3) == c3 + s * (d3 - c3)

    A x =  b where x(0) = t and x(1) = s
	(b1 - a1) * t - (d1 - c1) * s == c1 - a1
	(b2 - a2) * t - (d2 - c2) * s == c2 - a2
	(b3 - a3) * t - (d3 - c3) * s == c3 - a3
  */
  size_t numSegments = tubeMatrix->size1;
  double a1 = gsl_vector_get (segment, 0), a2 = gsl_vector_get (segment, 1), a3 = gsl_vector_get (segment, 2);
  double b1 = gsl_vector_get (segment, 3), b2 = gsl_vector_get (segment, 4), b3 = gsl_vector_get (segment, 5);
  double c1, c2, c3;
  double d1, d2, d3;
  
  gsl_matrix * A = gsl_matrix_alloc (3, 2);
  gsl_matrix * V = gsl_matrix_alloc (2, 2);
  gsl_vector * b = gsl_vector_alloc (3);
  gsl_vector * x = gsl_vector_alloc (2);
  gsl_vector * S = gsl_vector_alloc (2);
  gsl_vector * work = gsl_vector_alloc (2);
  gsl_vector * pt = gsl_vector_alloc (3);
  size_t numPts = 0;
  size_t ts = 0;
  
  for(size_t i = 0; i < numSegments; i++)
  {
    c1 = gsl_matrix_get (tubeMatrix, i, 0), c2 = gsl_matrix_get (tubeMatrix, i, 1), c3 = gsl_matrix_get (tubeMatrix, i, 2);
    d1 = gsl_matrix_get (tubeMatrix, i, 3), d2 = gsl_matrix_get (tubeMatrix, i, 4), d3 = gsl_matrix_get (tubeMatrix, i, 5);
    
	// printf ("a1 = %g, a2 = %g, a3 = %g\n", a1, a2, a3);
	// printf ("b1 = %g, b2 = %g, b3 = %g\n", b1, b2, b3);
	// printf ("c1 = %g, c2 = %g, c3 = %g\n", c1, c2, c3);
	// printf ("d1 = %g, d2 = %g, d3 = %g\n", d1, d2, d3);
  
    gsl_matrix_set(A, 0, 0, b1 - a1);
    gsl_matrix_set(A, 0, 1, - (d1 - c1));
    gsl_matrix_set(A, 1, 0, b2 - a2);
    gsl_matrix_set(A, 1, 1, - (d2 - c2));
	gsl_matrix_set(A, 2, 0, b3 - a3);
    gsl_matrix_set(A, 2, 1, - (d3 - c3));
    
    gsl_vector_set(b, 0, c1 - a1);
    gsl_vector_set(b, 1, c2 - a2);
	gsl_vector_set(b, 2, c3 - a3);
	
	// printf ("A = %g %g\n    %g %g\n    %g %g\n", 
	//  gsl_matrix_get(A, 0, 0), gsl_matrix_get(A, 0, 1), 
	//  gsl_matrix_get(A, 1, 0), gsl_matrix_get(A, 1, 1),
	//  gsl_matrix_get(A, 2, 0), gsl_matrix_get(A, 2, 1));
	//printf ("b = %g\n    %g\n    %g\n", 
	//  gsl_vector_get(b, 0), 
	//  gsl_vector_get(b, 1),
	//  gsl_vector_get(b, 2));
  
    //! A is M x N work is tmp storage of length N, A is replaced with U, V is N x N, S is M x N
    gsl_linalg_SV_decomp (A, V, S, work);
    gsl_linalg_SV_solve (A, V, S, b, x);

    // printf ("x = \n");
    // gsl_vector_fprintf(stdout, x, "%g");
	
	//! if x(0) = t, x(1) = s exist, then substitute and solve for overlapping point
	
	gsl_vector_set (pt, 0, a1 + gsl_vector_get(x, 0) * (b1 - a1));
	gsl_vector_set (pt, 1, a2 + gsl_vector_get(x, 0) * (b2 - a2));
	gsl_vector_set (pt, 2, a3 + gsl_vector_get(x, 0) * (b3 - a3));
    
	// printf ("pt = \n");
    // gsl_vector_fprintf(stdout, pt, "%g");
	
    if(!gsl_isnan (gsl_vector_get(x, 0)) && !gsl_isnan (gsl_vector_get(x, 1)))
	  //! the segments only overlap if their lengths are long enough
	  if(
	    //! A <= x <= B and C <= x <= D
	    min(a1, b1) <= gsl_vector_get(pt, 0) && 
		min(a2, b2) <= gsl_vector_get(pt, 1) && 
		min(a3, b3) <= gsl_vector_get(pt, 2) && 
		
		max(a1, b1) >= gsl_vector_get(pt, 0) && 
		max(a2, b2) >= gsl_vector_get(pt, 1) && 
		max(a3, b3) >= gsl_vector_get(pt, 2) && 
		
		min(c1, d1) <= gsl_vector_get(pt, 0) && 
		min(c2, d2) <= gsl_vector_get(pt, 1) && 
		min(c3, d3) <= gsl_vector_get(pt, 2) && 
		
		max(c1, d1) >= gsl_vector_get(pt, 0) && 
		max(c2, d2) >= gsl_vector_get(pt, 1) && 
		max(c3, d3) >= gsl_vector_get(pt, 2)
	  )
	  {
	    gsl_matrix_set(pts, numPts, 0, gsl_vector_get(pt, 0));
	    gsl_matrix_set(pts, numPts, 1, gsl_vector_get(pt, 1));
		gsl_matrix_set(pts, numPts, 2, gsl_vector_get(pt, 2));
	    // printf ("numPts(%ld) = %g %g %g\n", numPts, 
		//  gsl_matrix_get(pts, numPts, 0), 
		//  gsl_matrix_get(pts, numPts, 1),
		//  gsl_matrix_get(pts, numPts, 2));
		gsl_vector_set (tubeSegments, ts++, i);
	    numPts++;
	  }
	}
  
  gsl_vector_free (x);
  
  return numPts;
}

//! return newPos based upon Brownian motion from currentPos over timePeriod
int P1906MOL_MicrotubulesField::brownianMotion(gsl_rng * r, gsl_vector * currentPos, gsl_vector * newPos, double timePeriod)
{
  //! the new position is Gaussian with variance proportional to time taken: W_t - W_s ~ N(0, t - s)
  double sigma = timePeriod; /* sigma should be proportional to time */
  
  point (newPos, 
    gsl_vector_get (currentPos, 0) + gsl_ran_gaussian (r, sigma), /* x distance */
	gsl_vector_get (currentPos, 1) + gsl_ran_gaussian (r, sigma), /* y distance */
	gsl_vector_get (currentPos, 2) + gsl_ran_gaussian (r, sigma)  /* z distance */
  );
  
  return 0;
}

//! write the vector field in MATLAB format using regular spacing between samples in the file fname
//! The result of this file is loaded into Mathematica and the vector field is reconstructed from the samples via interpolation
//! Then the vector field operators are applied - see bushsf@research.ge.com for results
void P1906MOL_MicrotubulesField::vectorFieldMeshMATLAB(gsl_matrix * vf, const char *fname)
{
  FILE * pFile;

  pFile = fopen (fname,"w");
  
  // find the mesh volume limits
  double xMin, xMax;
  double yMin, yMax;
  double zMin, zMax;
  double uMin, uMax;
  double vMin, vMax;
  double wMin, wMax;
  
  gsl_vector * x = gsl_vector_alloc(vf->size1);
  gsl_vector * y = gsl_vector_alloc(vf->size1);
  gsl_vector * z = gsl_vector_alloc(vf->size1);
  gsl_vector * u = gsl_vector_alloc(vf->size1);
  gsl_vector * v = gsl_vector_alloc(vf->size1);
  gsl_vector * w = gsl_vector_alloc(vf->size1);
  
  //! extract the location and vector components from the field
  gsl_matrix_get_col (x, vf, 0);
  gsl_matrix_get_col (y, vf, 1);
  gsl_matrix_get_col (z, vf, 2);
  gsl_matrix_get_col (u, vf, 3);
  gsl_matrix_get_col (v, vf, 4);
  gsl_matrix_get_col (w, vf, 5);
  
  //! find the min and max values for each component
  gsl_vector_minmax (x, &xMin, &xMax);
  gsl_vector_minmax (y, &yMin, &yMax);
  gsl_vector_minmax (z, &zMin, &zMax);
  gsl_vector_minmax (u, &uMin, &uMax);
  gsl_vector_minmax (v, &vMin, &vMax);
  gsl_vector_minmax (w, &wMin, &wMax);
  
  printf ("xMin: %f yMin: %f zMin: %f uMin: %f vMin: %f wMin: %f\n", xMin, yMin, zMin, vMin, uMin, wMin);
  printf ("xMax: %f yMax: %f zMax: %f uMax: %f vMax: %f wMax: %f\n", xMax, yMax, zMax, vMax, uMax, wMax);
	
  gsl_vector * pt1 = gsl_vector_alloc (3);
  gsl_vector * pt2 = gsl_vector_alloc (3);
  gsl_vector * vec = gsl_vector_alloc (3);
  gsl_vector * closest = gsl_vector_alloc (6);
  
  double xStepsize = (xMax - xMin) / 10.0;
  double yStepsize = (yMax - yMin) / 10.0;
  double zStepsize = (zMax - zMin) / 10.0;
  
  printf ("xStepsize: %f yStepsize: %f zStepsize: %f\n", xStepsize, yStepsize, zStepsize);
  
  //! step through equidistant points in a volume and store the vector values at each point
  for (double i = xMin; i < xMax; i += xStepsize)
    for (double j = yMin; j < yMax; j += yStepsize)
	  for (double k = zMin; k < zMax; k += zStepsize)
      {
	    //! find the closest point to the current location
	    point (pt1, i, j, k);
		// printf ("pt1\n");
		// displayPoint (pt1);
	    findClosestPoint (pt1, vf, closest);
		//printf ("(findClosestPoint) closest vector: %f %f %f %f %f %f\n",
        //  gsl_vector_get (closest, 0),
	    //  gsl_vector_get (closest, 1),
	    //  gsl_vector_get (closest, 2),
	    //  gsl_vector_get (closest, 3),
	    //  gsl_vector_get (closest, 4),
	    //  gsl_vector_get (closest, 5));
		point (pt2, 
		  gsl_vector_get (closest, 0),
		  gsl_vector_get (closest, 1),
		  gsl_vector_get (closest, 2));
		  
		// printf ("(vectorFieldMeshMATLAB) distance between pt and vector location: %f\n", distance (pt1, pt2));
		//! check if distance within range
		if (distance (pt1, pt2) > 2.0 * xStepsize)
		{
		  //! if not, store the null vector
		  // printf ("(vectorFieldMeshMATLAB) using null vector\n");
          point (vec, 
		    0.0, 
			0.0, 
			0.0);
		  // displayPoint (vec);
	    } 
		else 
		{
		  // printf ("(vectorFieldMeshMATLAB) using vector\n");
		  //! otherwise, store the vector value
		  point (vec,
	        gsl_vector_get (closest, 3),
	        gsl_vector_get (closest, 4),
	        gsl_vector_get (closest, 5));
		  // displayPoint (vec);
		}
		//! print current location and stored vector value
        fprintf (pFile, "%f %f %f %f %f %f\n",
	      i, 
		  j,
		  k, 
		  gsl_vector_get (vec, 0),
          gsl_vector_get (vec, 1),
          gsl_vector_get (vec, 2));
      }
	  
  fclose(pFile);
}

//! return the location of the vector from vf that is closest to the point pt in result
void P1906MOL_MicrotubulesField::findClosestPoint(gsl_vector * pt, gsl_matrix * vf, gsl_vector *result)
{
  //! the current closest point
  gsl_vector * cpt = gsl_vector_alloc (3);
  //! the current closest vector 
  gsl_vector * cv = gsl_vector_alloc (3);
  //! the test point
  gsl_vector * tpt = gsl_vector_alloc (3);
  
  //! display the vector field
  // for (size_t i = 0; i < vf->size1; i++)
  // {
  //  printf ("vf: %f %f %f %f %f %f\n",
  //	  gsl_matrix_get (vf, i, 0),
  //	  gsl_matrix_get (vf, i, 1),
  //	  gsl_matrix_get (vf, i, 2),
  //	  gsl_matrix_get (vf, i, 3),
  //	  gsl_matrix_get (vf, i, 4),
  //	  gsl_matrix_get (vf, i, 5));
  // }
  
  point (cpt, GSL_POSINF, GSL_POSINF, GSL_POSINF);
  
  for (size_t i = 0; i < vf->size1; i++)
  {
    //! test the next vector field location 
    point (tpt, 
	  gsl_matrix_get (vf, i, 0),
	  gsl_matrix_get (vf, i, 1),
	  gsl_matrix_get (vf, i, 2));
	
	// printf ("pt\n");
	// displayPoint (pt);
	// printf ("tpt\n");
	// displayPoint (tpt);
	// printf ("cpt\n");
	// displayPoint (cpt);
	
	//! tpt is point under test from vector field
	//! cpt is the current shortest distance point
	// printf ("(findClosestPoint) distance to tpt: %f\n", distance (tpt, pt));
	// printf ("(findClosestPoint) distance to cpt: %f\n", distance (cpt, pt));
	
	//! update the location and vector if closer
	if ((distance (tpt, pt) < distance (cpt, pt)) || gsl_isnan(distance (cpt, pt)))
	{
	  // printf ("storing point\n");
	  point (cpt, 
	    gsl_vector_get (tpt, 0),
	    gsl_vector_get (tpt, 1),
	    gsl_vector_get (tpt, 2));
	  // printf ("cpt\n");
	  // displayPoint (cpt);
      point (cv, 
	    gsl_matrix_get (vf, i, 3),
	    gsl_matrix_get (vf, i, 4),
	    gsl_matrix_get (vf, i, 5));
	  // printf ("cv\n");
	  // displayPoint (cv);
	}
  }
  
  // printf ("(findClosestPoint) final cv\n");
  // displayPoint (cv);
  
  //! return the result
  gsl_vector_set (result, 0, gsl_vector_get (cpt, 0));
  gsl_vector_set (result, 1, gsl_vector_get (cpt, 1));
  gsl_vector_set (result, 2, gsl_vector_get (cpt, 2));
  gsl_vector_set (result, 3, gsl_vector_get (cv, 0));
  gsl_vector_set (result, 4, gsl_vector_get (cv, 1));
  gsl_vector_set (result, 5, gsl_vector_get (cv, 2));
  
  // printf ("(findClosestPoint) result: %f %f %f %f %f %f\n",
  //  gsl_vector_get (result, 0),
  //	gsl_vector_get (result, 1),
  //	gsl_vector_get (result, 2),
  //	gsl_vector_get (result, 3),
  //	gsl_vector_get (result, 4),
  //	gsl_vector_get (result, 5));
}

//! write a list of vectors into file fname in MATLAB loadable format
void P1906MOL_MicrotubulesField::vectorFieldPlotMATLAB(gsl_matrix * vf, const char *fname)
{
  FILE * pFile;

  pFile = fopen (fname,"w");
  
  for (size_t i = 0; i < vf->size1; i++)
  {
    fprintf (pFile, "%f %f %f %f %f %f\n", 
	  gsl_matrix_get (vf, i, 0),
	  gsl_matrix_get (vf, i, 1),
	  gsl_matrix_get (vf, i, 2),
	  gsl_matrix_get (vf, i, 3),
	  gsl_matrix_get (vf, i, 4),
	  gsl_matrix_get (vf, i, 5));
  }
  
  fclose(pFile);
}

//! print all tubes in tubeMatrix into a Mathematica file fname with segments per tube of segPerTube
void P1906MOL_MicrotubulesField::tubes2Mma(gsl_matrix *tubeMatrix, size_t segPerTube, const char* fname)
{
  /** save tubes to file in the form of 
    GraphPlot3D[{1 -> 2, 1 -> 4, 1 -> 5, 2 -> 3, 2 -> 6, 3 -> 4, 3 -> 7, 4 -> 8, 5 -> 6, 5 -> 8, 6 -> 7, 7 -> 8}, 
	VertexCoordinateRules -> {1 -> {0, 1, 2}, 2 -> {-1, 0, 2}, 3 -> {0, -1, 2}, 4 -> {1, 0, 2}, 5 -> {0, 2, 0}, 6 -> {-2, 0, 0}, 7 -> {0, -2, 0}, 8 -> {2, 0, 0}}]
  */
  FILE * pFile;

  pFile = fopen (fname,"w");
  
  size_t numSegments = tubeMatrix->size1;
  size_t pt = 1;
  size_t numTubes = numSegments / segPerTube;
  
  //! fprintf (pFile, "numTubes: %ld segPerTube: %ld\n", numTubes, segPerTube);
  
  fprintf (pFile, "GraphPlot3D[{");
  for (size_t i = 0; i < numTubes; i++)
  {
    for (size_t j = 0; j < segPerTube; j++)
	{
      fprintf (pFile, "%ld -> ", pt);
	  pt++;
	  fprintf (pFile, "%ld", pt);
      if ((i < (numTubes - 1)) || (j < (segPerTube - 1))) fprintf (pFile, ", ");
	}
    pt++;
  }
  fprintf (pFile, "}, ");
  pt = 1;
  fprintf (pFile, "VertexCoordinateRules ->{");
    for (size_t i = 0; i < numTubes; i++)
	{
      for (size_t j = 0; j < segPerTube; j++)
	  {
	    if (j == 0) //! only print the ends after the first one
	    {
          fprintf (pFile, "%ld -> {%g, %g, %g}, ",
		    pt,
	        gsl_matrix_get(tubeMatrix, i * segPerTube + j, 0),
		    gsl_matrix_get(tubeMatrix, i * segPerTube + j, 1),
		    gsl_matrix_get(tubeMatrix, i * segPerTube + j, 2));
		  pt++;
        }
		fprintf (pFile, "%ld -> {%g, %g, %g}", 
		  pt,
		  gsl_matrix_get(tubeMatrix, i * segPerTube + j, 3),
		  gsl_matrix_get(tubeMatrix, i * segPerTube + j, 4),
		  gsl_matrix_get(tubeMatrix, i * segPerTube + j, 5));
		pt++;
	    if ((i < (numTubes - 1)) || (j < (segPerTube - 1))) fprintf (pFile, ", ");
	  }
	}
  fprintf (pFile, "}]\n");
  //! option to print vertex labels
  //! fprintf (pFile, "}, VertexLabeling -> True]\n");
  
  fclose (pFile);
}

//! create and return a microtubule of given persistence length in segMatrix and its structural entropy in se
int P1906MOL_MicrotubulesField::genTube(struct tubeCharacteristcs_t * ts, gsl_rng *r, gsl_matrix *segMatrix, gsl_vector *startPt)
{
  /**
    for 3D consider generating two angles per tube: \theta and \psi in spherical coordinates, length is already specified 
	x = r sin \theta cos \psi
	y = r sin \theta sin \psi
	z = r cos \theta
	
	\todo randomize initial tube orientations
  */
   
  gsl_matrix * segAngleTheta = gsl_matrix_alloc (ts->numSegments, 1);
  gsl_matrix * segAnglePsi = gsl_matrix_alloc (ts->numSegments, 1);

  // printf ("(genTube) startX = %g startY = %g startZ = %g numSegments = %ld segLength = %g persistenceLength = %g\n", 
  //  gsl_vector_get (startPt, 0), 
  //	gsl_vector_get (startPt, 1),  
  //	gsl_vector_get (startPt, 2), 
  //	ts->numSegments, 
  //	ts->segLength, 
  //	ts->persistenceLength);
  
  genPersistenceLength(r, segAngleTheta, ts->segLength, ts->persistenceLength);
  genPersistenceLength(r, segAnglePsi, ts->segLength, ts->persistenceLength);
  
  ts->se = sEntropy (segAngleTheta) + sEntropy (segAnglePsi);
  
  // for (size_t i = 0; i < ts->numSegments; i++)
  //  for (size_t j = 0; j < 1; j++)
  //	{
  //	  printf ("segAngleTheta(%ld,%ld) = %g\n", i, j, gsl_matrix_get (segAngleTheta, i, j));
  //	  printf ("segAnglePsi(%ld,%ld) = %g\n", i, j, gsl_matrix_get (segAnglePsi, i, j));
  //  }
	
  // printf ("(genTube) segments allocated %ld %ld numSegments %ld\n", segMatrix->size1, segMatrix->size2, ts->numSegments);
  
  for (size_t i = 0; i < ts->segPerTube; i++)
  {
    //! set the start points of the segment
    if(i == 0)
    { //! use starting coordinates
      gsl_matrix_set(segMatrix, 0, 0, gsl_vector_get (startPt, 0));
      gsl_matrix_set(segMatrix, 0, 1, gsl_vector_get (startPt, 1));
	  gsl_matrix_set(segMatrix, 0, 2, gsl_vector_get (startPt, 2));
    }
    else
    { //! use end of last segment
      gsl_matrix_set(segMatrix, i, 0, gsl_matrix_get(segMatrix, i - 1, 3));
      gsl_matrix_set(segMatrix, i, 1, gsl_matrix_get(segMatrix, i - 1, 4));
	  gsl_matrix_set(segMatrix, i, 2, gsl_matrix_get(segMatrix, i - 1, 5));
    }
	//! set the end points of the segment
	double x = ts->segLength * sin(gsl_matrix_get(segAngleTheta, i, 0)) * cos(gsl_matrix_get(segAnglePsi, i, 0));
	double y = ts->segLength * sin(gsl_matrix_get(segAngleTheta, i, 0)) * sin(gsl_matrix_get(segAnglePsi, i, 0));
	double z = ts->segLength * cos(gsl_matrix_get(segAngleTheta, i, 0));
	gsl_matrix_set(segMatrix, i, 3, x + gsl_matrix_get(segMatrix, i, 0));
    gsl_matrix_set(segMatrix, i, 4, y + gsl_matrix_get(segMatrix, i, 1));
	gsl_matrix_set(segMatrix, i, 5, z + gsl_matrix_get(segMatrix, i, 2));
  }
  
  return 0;
}

//! generate angles for a structure with the given persistence length and segment length and return the angles in setAngle
double P1906MOL_MicrotubulesField::genPersistenceLength(gsl_rng *r, gsl_matrix *segAngle, double segLength, double persistenceLength)
{
  /** 
    the angle distribution is Gaussian with zero mean and variance \f$\sigma^2 = \sqrt(2 \Delta s / l_p)\f$,
    where \f$l_p\f$ is persistence length and \f$\Delta\f$ is the segment length; angles are in radians
	
	See http://www.uic.edu/classes/phys/phys450/MARKO/N015.html.
	
	Fortunately there is a relation that we can take advantage of to figure out the stiffness of the torsion spring 
	between each segment. This relation is \f$K = L_p k_BT\f$. From this we can see that the persistence length is a characteristic 
	length scale that relates the bending rigidity of the DNA chain to the thermal energy. Larger persistence lengths indicate 
	a stiffer chain that is less susceptible to thermal fluctuations.

    What we're interested in, of course, is the energy cost for bending the DNA into the nice round coil it must be in to 
	fit in the capsid. Like a linear spring where the energy is of the form \f$E = 1/2 k x^2\f$, the energy to bend the DNA chain 
	is of the form \f$E = 1/2 K L (dt/ds)^2\f$, where \f$dt/ds\f$ is the derivative of the tangent vector to each Kuhn segment of the 
	chain. Since the chain is being bent into a circle of radius \f$R\f$, \f$dt/ds = 1/R\f$, and thus the bending energy is 
	\f$E = 1/2 K L/R^2 = 1/2 L_p L/R^2 k_BT\f$. If we assume that all of the DNA is coiled around the outer radius of the 
	capsid (R = 42 nm), then the energy cost due to bending is \f$E_{bend} = 1/2 * (50 nm)*(6800 nm)/(42 nm)2 k_BT ~ 10^2 k_BT\f$. 
	Note that this is approximately the same energy cost as we found for the loss of entropy. From the model shown on 
	the first page of the website, you may have noticed that some of the DNA is in fact coiled at smaller radii than the 
	outermost "shell". The energy cost for bending that part of the DNA will be much higher since E is proportional to \f$1/R^2\f$, 
	but even so, it is unlikely that this energy cost greatly exceeds \f$~10^3 k_BT\f$.
	
	Ref: Multiscale Modeling in Biomechanics and Mechanobiology edited by Suvranu De, Wonmuk Hwang, Ellen Kuh pp. 68-69
  */
  double angle; 

  for(size_t i = 0; i < segAngle->size1; i++)
  {
    angle = gsl_ran_gaussian (r, sqrt(2.0 * segLength / persistenceLength));
	//! a radian is 180/Pi degrees
	// printf ("angle(%ld) = %g\n", i, angle * 180.00 / M_PI);
    gsl_matrix_set(segAngle, i, 0, angle);
  }
    
  sEntropy(segAngle);
  
  return 0;
}

//! compute the persistence length of a set of segments
//! \todo never implemented; may not be needed
double P1906MOL_MicrotubulesField::getPersistenceLength(gsl_matrix *segMatrix)
{
  //! segMatrix rows are x_start x_end y_start y_end (may not need numSegments)
  
  //! may not need this method
  
  return 0;
}

//! return the information entropy of a tube segment defined by its list of angles in segAngle
double P1906MOL_MicrotubulesField::sEntropy(gsl_matrix *segAngle)
{
  //! \f$ H(x) = - sum( P(x) \log P(x) ) \f$
  //! bin the values in order to find P(x)
  //! see gsl_histogram_pdf * gsl_histogram_pdf_alloc (size_t n)
  //! see https://www.gnu.org/software/gsl/manual/html_node/The-histogram-probability-distribution-struct.html
  
  gsl_histogram * h = gsl_histogram_alloc (100); /* need to determine optimal number of bins */
  //! get the max segAngle
  double minAngle = gsl_matrix_min(segAngle);
  double maxAngle = gsl_matrix_max(segAngle);
  gsl_histogram_set_ranges_uniform (h, minAngle, maxAngle);

  for (size_t i = 0; i < segAngle->size1; i++)
    for (size_t j = 0; j < 1; j++)
	{
	  //! printf ("segAngle(%ld,%ld) = %g\t", i, j, gsl_matrix_get (segAngle, i, j)); */
	  gsl_histogram_increment (h, gsl_matrix_get (segAngle, i, j));
    }

  double H = 0;
  int numBins = gsl_histogram_bins (h);
  for(int i = 0; i < numBins; i++)
  {
	double p = gsl_histogram_get (h, i) / gsl_histogram_sum(h);
	//! printf("p = %g\n", p);
	if (p > 0)
	{
	  H = H - p * gsl_sf_log(p);
	  //! printf("H = %g\n", H);
	}
  }
	
  // printf("H = %g\n", H);
  //! gsl_histogram_fprintf (stdout, h, "%g", "%g");
  gsl_histogram_free (h);
  
  return H;
}

P1906MOL_MicrotubulesField::~P1906MOL_MicrotubulesField ()
{
  NS_LOG_FUNCTION (this);
}

} // namespace ns3

/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 *  Copyright � 2014 by IEEE.
 *
 *  This source file is an essential part of IEEE P1906.1,
 *  Recommended Practice for Nanoscale and Molecular
 *  Communication Framework.
 *  Verbatim copies of this source file may be used and
 *  distributed without restriction. Modifications to this source
 *  file as permitted in IEEE P1906.1 may also be made and
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

/*
 * Description:
 * this file models the molecular communication based on microtubules structures in class P1906MOL_MicrotubulesField
 */
 
//! \todo (DONE) set microtubule field parameters from within microtubules-example.cc (expose the parameters to the command line)
//! \todo (WRITTEN) pass the microtubule structures to both instances of field for the Nodes?
//! \todo (DONE) remove the ts structure from the set methods
//! \todo access diffusion (MOL) parameters from extended model (currently in Specificity; may need to move to motion)
//! \todo ensure that brownianMotion and motorWalk are using realistic movement rates
//! \todo plot results
//! \todo change all printfs to ns-3 LOG output
//! \todo review documentation of each method
//! \todo (OK AS-IS) modify the unique float2destination.mma file names to be shorter
//! \todo modify _RUN_MOTOR_CHANNEL_CAPACITY_.sh to generate plot of distance versus Brownian motion
//! \todo extra-credit: modify _RUN_MOTOR_CHANNEL_CAPACITY_.sh to generate plot of distance versus Brownian motion WITH TUBES @ given orientations
//! \todo extra-credit: modify _RUN_MOTOR_CHANNEL_CAPACITY_.sh to generate plot of structural entropy WITH TUBES @ given orientations
//! \todo update/complete the extensions/README.txt file
//! \todo rerun Doxygen for final SVN update
 
//! \details 1906 Component map to the molecular motor instantiation
//! <pre>
//!  1906 Component             Molecular Motor 
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

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "ns3/p1906-helper.h"
#include "ns3/p1906-net-device.h"
#include "ns3/p1906-medium.h"

#include "ns3/p1906-mol-specificity.h"

#include "ns3/p1906-mol-extended-perturbation.h"
#include "ns3/p1906-mol-extended-field.h"
#include "ns3/p1906-mol-extended-motion.h"
#include "ns3/p1906-mol-field-microtubule.h"
#include "ns3/p1906-mol-extended-communication-interface.h"
#include "ns3/p1906-mol-extended-transmitter-communication-interface.h"
#include "ns3/p1906-mol-extended-receiver-communication-interface.h"

using namespace ns3;


int main (int argc, char *argv[])
{	

  //set of parameters
  double nodeDistance = 0.001; 								//  [m]
  double nbOfMoleculas = 50000; 							//  [pJ]
  double pulseInterval = 1.;								//  [ms]
  double diffusionCoefficient = 1;							//  [nm^2/ns]
  double tube_volume = 25; 									// [nm^3]
  double mean_tube_length = 100; 							// [nm]
  double mean_intra_tube_angle = 30; 						// [degrees]
  double mean_inter_tube_angle = 10; 						// [degrees]
  double mean_tube_density = 10; 							// [tube segments/nm^3]
  double tube_persistenceLength = 50; 						// [nm]
  size_t segPerTube = 10; 									// [segments/microtubule]

  CommandLine cmd;
  cmd.AddValue("nodeDistance", "nodeDistance", nodeDistance);
  cmd.AddValue("nbOfMoleculas", "nbOfMoleculas", nbOfMoleculas);
  cmd.AddValue("diffusionCoefficient", "diffusionCoefficient", diffusionCoefficient);
  cmd.AddValue("tube_volume", "tube_volume", tube_volume);
  cmd.AddValue("mean_tube_length", "mean_tube_length", mean_tube_length);
  cmd.AddValue("mean_intra_tube_angle", "mean_intra_tube_angle", mean_intra_tube_angle);
  cmd.AddValue("mean_inter_tube_angle", "mean_inter_tube_angle", mean_inter_tube_angle);
  cmd.AddValue("mean_tube_density", "mean_tube_density", mean_tube_density);
  cmd.AddValue("tube_persistenceLength", "tube_persistenceLength", tube_persistenceLength);
  cmd.AddValue("segPerTube", "tube_persistenceLength", segPerTube);
  cmd.Parse(argc, argv);

  Time::SetResolution(Time::NS);

  // Create P1906 Helper
  P1906Helper helper;
  helper.EnableLogComponents ();

  // Create nodes (typical operation of ns-3)
  NodeContainer n;
  NetDeviceContainer d;
  n.Create (2);

  // Create a medium and the Motion component
  Ptr<P1906Medium> medium = CreateObject<P1906Medium> ();
  //! may need to set other motion properties...
  Ptr<P1906MOL_ExtendedMotion> motion = CreateObject<P1906MOL_ExtendedMotion> ();
  
  motion->SetDiffusionCoefficient (diffusionCoefficient);
  medium->SetP1906Motion (motion);

  // Create Device 1 and related components/entities
  Ptr<P1906NetDevice> dev1 = CreateObject<P1906NetDevice> ();
  Ptr<P1906MOL_ExtendedCommunicationInterface> c1 = CreateObject<P1906MOL_ExtendedCommunicationInterface> ();
  Ptr<P1906MOLSpecificity> s1 = CreateObject<P1906MOLSpecificity> ();
  Ptr<P1906MOL_MicrotubulesField> fi1 = CreateObject<P1906MOL_MicrotubulesField> ();
  fi1->setTubeVolume (volume);
  fi1->setTubeLength (mean_tube_length);
  fi1->setTubeIntraAngle (mean_intra_tube_angle);
  fi1->setTubeInterAngle (mean_inter_tube_angle);
  fi1->setTubeDensity (mean_tube_density);
  fi1->setTubePersistenceLength (persistenceLength);
  fi1->setTubeSegments (segPerTube);
  
  //! this class creates the motor (message carrier)
  Ptr<P1906MOL_ExtendedPerturbation> p1 = CreateObject<P1906MOL_ExtendedPerturbation> ();
  //! don't need these properties...
  p1->SetPulseInterval (MilliSeconds(pulseInterval));
  p1->SetMolecules (nbOfMoleculas);
  s1->SetDiffusionCoefficient (diffusionCoefficient);
  
  printf ("(microtubules-example) Device 1 created\n");

  // Create Device 2 and related components/entities
  Ptr<P1906NetDevice> dev2 = CreateObject<P1906NetDevice> ();
  Ptr<P1906MOL_ExtendedCommunicationInterface> c2 = CreateObject<P1906MOL_ExtendedCommunicationInterface> ();
  Ptr<P1906MOLSpecificity> s2 = CreateObject<P1906MOLSpecificity> ();
  //! SFB: this class creates the microtubules
  Ptr<P1906MOL_MicrotubulesField> fi2 = CreateObject<P1906MOL_MicrotubulesField> ();
  //! this class creates the motor (message carrier)
  Ptr<P1906MOL_ExtendedPerturbation> p2 = CreateObject<P1906MOL_ExtendedPerturbation> ();
  //! don't need these properties...
  p2->SetPulseInterval (MilliSeconds(pulseInterval));
  p2->SetMolecules (nbOfMoleculas);
  s2->SetDiffusionCoefficient (diffusionCoefficient);
  
  printf ("(microtubules-example) Device 2 created\n");

  //set devices positions
  Ptr<ListPositionAllocator> positionAlloc =
		  CreateObject<ListPositionAllocator> ();
  positionAlloc->Add (Vector(0, 0, 0));
  positionAlloc->Add (Vector(nodeDistance, 0, 0));
  MobilityHelper mobility;
  mobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
  mobility.SetPositionAllocator(positionAlloc);
  mobility.Install(n);

  // Connect devices, nodes, medium, components and entities
  d.Add (dev1);
  d.Add (dev2);
  helper.Connect(n.Get (0), dev1, medium, c1, fi1, p1, s1);
  helper.Connect(n.Get (1), dev2, medium, c2, fi2, p2, s2);
  
  printf ("(microtubules-example) Connected devices, nodes, medium, components and entities\n");

  // Create a message to send into the network
  int pktSize = 1; //bytes
  uint8_t *buffer  = new uint8_t[pktSize];
  for (int i = 0; i < pktSize; i++)
    {
	  buffer[i] = 0; //empty information
    }
  
  Ptr<Packet> message = Create<Packet>(buffer, pktSize);
  
  printf ("(microtubules-example) Packet created\n");

  //! c1 is the P1906MOLCommunicationInterface for Node 1
  c1->HandleTransmission (message);
  
  printf ("(microtubules-example) c1->HandleTransmission (message)\n");

  Simulator::Stop (Seconds (0.01));
  Simulator::Run ();

  Simulator::Destroy ();
  return 0;
}

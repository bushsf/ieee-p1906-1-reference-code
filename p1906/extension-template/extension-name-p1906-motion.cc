/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 *  Copyright � 2015 by IEEE.
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
 * Author: Giuseppe Piro - Telematics Lab Research Group
 *                         Politecnico di Bari
 *                         giuseppe.piro@poliba.it
 *                         telematics.poliba.it/piro
 */


#include "ns3/log.h"

#include "extension-name-p1906-motion.h"
#include "extension-name-p1906-communication-interface.h"
#include "extension-name-p1906-message-carrier.h"
#include "extension-name-p1906-field.h"

namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("ExtensionNameP1906Motion");

TypeId ExtensionNameP1906Motion::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::ExtensionNameP1906Motion")
    .SetParent<P1906Motion> ();
  return tid;
}

ExtensionNameP1906Motion::ExtensionNameP1906Motion ()
{
  NS_LOG_FUNCTION (this);
}

ExtensionNameP1906Motion::~ExtensionNameP1906Motion ()
{
  NS_LOG_FUNCTION (this);
}


double
ExtensionNameP1906Motion::ComputePropagationDelay (Ptr<P1906CommunicationInterface> src,
		                             Ptr<P1906CommunicationInterface> dst,
		                             Ptr<P1906MessageCarrier> message,
		                             Ptr<P1906Field> field)
{
  NS_LOG_FUNCTION (this << "Return the defaul value: 0s");
  return 0.;
}


Ptr<P1906MessageCarrier>
ExtensionNameP1906Motion::CalculateReceivedMessageCarrier(Ptr<P1906CommunicationInterface> src,
		                                    Ptr<P1906CommunicationInterface> dst,
		                                    Ptr<P1906MessageCarrier> message,
		                                    Ptr<P1906Field> field)
{
  NS_LOG_FUNCTION (this << "apply default behavior: do nothing");
  return message;
}


} // namespace ns3

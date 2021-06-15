//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in ECCE_DEMP.h.
//
// ECCE_DEMP(const std::string &name = "ECCE_DEMP")
// everything is keyed to ECCE_DEMP, duplicate names do work but it makes
// e.g. finding culprits in logs difficult or getting a pointer to the module
// from the command line
//
// ECCE_DEMP::~ECCE_DEMP()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the node tree
//
// int ECCE_DEMP::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer. You
// can create historgrams here or put objects on the node tree but be aware that
// modules which haven't been registered yet did not put antyhing on the node tree
//
// int ECCE_DEMP::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's action
// depends on what else is around. Last chance to put nodes under the DST Node
// We mix events during readback if branches are added after the first event
//
// int ECCE_DEMP::process_event(PHCompositeNode *topNode)
// called for every event. Return codes trigger actions, you find them in
// $OFFLINE_MAIN/include/fun4all/Fun4AllReturnCodes.h
//   everything is good:
//     return Fun4AllReturnCodes::EVENT_OK
//   abort event reconstruction, clear everything and process next event:
//     return Fun4AllReturnCodes::ABORT_EVENT; 
//   proceed but do not save this event in output (needs output manager setting):
//     return Fun4AllReturnCodes::DISCARD_EVENT; 
//   abort processing:
//     return Fun4AllReturnCodes::ABORT_RUN
// all other integers will lead to an error and abort of processing
//
// int ECCE_DEMP::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs clearing
// after each event, this is the place to do that. The nodes under the DST node are cleared
// by the framework
//
// int ECCE_DEMP::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int ECCE_DEMP::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int ECCE_DEMP::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void ECCE_DEMP::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
//____________________________________________________________________________..

#include "ECCE_DEMP.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

#include <stdio.h>

#include <fun4all/Fun4AllHistoManager.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// G4Hits includes
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include <TFile.h>
#include <TNtuple.h>

#include <cassert>
#include <sstream>
#include <string>
#include <iostream>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

using namespace std;


//____________________________________________________________________________..
//ECCE_DEMP::ECCE_DEMP(const std::string &name):
// SubsysReco(name)
//{
//  std::cout << "ECCE_DEMP::ECCE_DEMP(const std::string &name) Calling ctor" << std::endl;
//}


ECCE_DEMP::ECCE_DEMP(const std::string &name, const std::string& filename):
 SubsysReco(name)
 , outfilename(filename)
{
  std::cout << "ECCE_DEMP_example::Diff_Tagg_example(const std::string &name) Calling ctor" << std::endl;


  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(m_RandomGenerator, seed);

}




//____________________________________________________________________________..
ECCE_DEMP::~ECCE_DEMP()
{

  gsl_rng_free(m_RandomGenerator);

  std::cout << "ECCE_DEMP::~ECCE_DEMP() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int ECCE_DEMP::Init(PHCompositeNode *topNode)
{

  hm = new Fun4AllHistoManager(Name());
  // create and register your histos (all types) here
  // TH1 *h1 = new TH1F("h1",....)
  // hm->registerHisto(h1);
  outfile = new TFile(outfilename.c_str(), "RECREATE");
  g4hitntuple = new TNtuple("hitntup", "G4Hits", "x0:y0:z0:x1:y1:z1:edep");

  std::cout << "ECCE_DEMP::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  event_itt = 0;


  h2_ZDC_XY = new TH2F("ZDC_XY", "ZDC XY", 200, -50, 50, 200, -50, 50);

  h2_ZDC_XY_double = new TH2F("ZDC_XY_double", "ZDC XY Double gamma", 200, -50, 50, 200, -50, 50);

  h1_E_dep = new TH1F("E_dep", "E Dependence", 120, 0.0, 60.0);

  h1_E_dep_smeared = new TH1F("E_dep_smeared", "E Dependence Smeared", 120, 0.0, 60.0);



  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP::InitRun(PHCompositeNode *topNode)
{
  std::cout << "ECCE_DEMP::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP::process_event(PHCompositeNode *topNode)
{
//  std::cout << "ECCE_DEMP::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;

  
  ZDC_hit = 0;

  event_itt++; 
 
  if(event_itt%100 == 0)
     std::cout << "Event Processing Counter: " << event_itt << endl;

  process_g4hits_ZDC(topNode);

  process_g4hits_RomanPots(topNode);


  ////-------------------------
  ////Example for Getting the Hadron end cap hits and clusters
  //// uncommenting the following line:
  ///// for inner Hadron end cap  
  // process_g4hits(topNode, "HCALIN");
  // process_g4clusters(topNode, "HCALIN");
  //
  ///// for outer Hadron end cap  
  // process_g4hits(topNode, "HCALOUT");
  // process_g4clusters(topNode, "HCALOUT");
  //
  //// for Electron end cap  
  // process_g4hits(topNode, "EEMC");
  // process_g4clusters(topNode, "EEMC");
  //

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP::ResetEvent(PHCompositeNode *topNode)
{
//  std::cout << "ECCE_DEMP::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
//
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP::EndRun(const int runnumber)
{
  std::cout << "ECCE_DEMP::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP::End(PHCompositeNode *topNode)
{
  std::cout << "ECCE_DEMP::End(PHCompositeNode *topNode) This is the End..." << std::endl;


//  h2_ZDC_XY->Write();
//  h2_ZDC_XY_double->Write();
//  
//  h1_E_dep->Write();
//  h1_E_dep_smeared->Write();

  outfile->cd();
  g4hitntuple->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP::Reset(PHCompositeNode *topNode)
{
 std::cout << "ECCE_DEMP::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void ECCE_DEMP::Print(const std::string &what) const
{
  std::cout << "ECCE_DEMP::Print(const std::string &what) const Printing info for " << what << std::endl;
}


//***************************************************
//


int ECCE_DEMP::process_g4hits_ZDC(PHCompositeNode* topNode)
{
  ostringstream nodename;

  // loop over the G4Hits
  nodename.str("");
//  nodename << "G4HIT_" << detector;
//  nodename << "G4HIT_" << "ZDC";
  nodename << "G4HIT_" << "ZDC";
//  nodename << "G4HIT_" << "EEMC";

  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());

//  cout << "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa " << endl;


  float smeared_E;


  if (hits) {
//    // this returns an iterator to the beginning and the end of our G4Hits
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)

    {
//	HIT_IN_ZDC=true;
	ZDC_hit++;
    }

    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++) {


//	ZDC_hit++;

//      cout << "AAA" << endl;
      // the pointer to the G4Hit is hit_iter->second
      g4hitntuple->Fill(hit_iter->second->get_x(0),
                        hit_iter->second->get_y(0),
                        hit_iter->second->get_z(0),
                        hit_iter->second->get_x(1),
                        hit_iter->second->get_y(1),
                        hit_iter->second->get_z(1),
                        hit_iter->second->get_edep());



//      cout << hit_iter->second->get_x(0)-90 << "   " << hit_iter->second->get_y(0) << endl;


      h2_ZDC_XY->Fill(hit_iter->second->get_x(0)-90, hit_iter->second->get_y(0)); 
//
      smeared_E = EMCAL_Smear(hit_iter->second->get_edep());
//
      if (ZDC_hit == 2 ) {

//      cout << hit_iter->second->get_x(0)-90 << "   " << hit_iter->second->get_y(0) << endl;

        h2_ZDC_XY_double->Fill(hit_iter->second->get_x(0)-90, hit_iter->second->get_y(0)); 
//      h1_E_dep->Fill(hit_iter->second->get_edep()); 
//
        h1_E_dep->Fill(hit_iter->second->get_edep()); 
        h1_E_dep_smeared->Fill(smeared_E); 
//
      }
//
//
//
////	hit_iter->get_avg_t();

    }
  }

//  cout << "BB" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}


//***************************************************
// Getting the RomanPots hits

int ECCE_DEMP::process_g4hits_RomanPots(PHCompositeNode* topNode)
{
  ostringstream nodename;


  cout << "Entering Romanpot?" << endl;

  // loop over the G4Hits
  nodename.str("");
//  nodename << "G4HIT_" << detector;
//  nodename << "G4HIT_" << "ZDC";
  nodename << "G4HIT_" << "RomanPots_0";
//  nodename << "G4HIT_" << "EEMC";

  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());


  if (hits) {
//    // this returns an iterator to the beginning and the end of our G4Hits
    PHG4HitContainer::ConstRange hit_range = hits->getHits();

    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++) {


	cout << "Roman pot hits? " << endl;
	cout << "This is where you can fill your loop " << endl;

      }
    }

  return Fun4AllReturnCodes::EVENT_OK;
}


//*****************************************************

float ECCE_DEMP::EMCAL_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.45*.45/E + 0.075*0.075);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}


//*****************************************************

float ECCE_DEMP::HCAL_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.50*.50/E + 0.1*0.1);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}

//*****************************************************

float ECCE_DEMP::PbWO4_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

//*****************************************************

float ECCE_DEMP::Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.1;         /// Position resolution 0.1 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}




int ECCE_DEMP::process_g4hits(PHCompositeNode* topNode, const std::string& detector)
{
  ostringstream nodename;

  // loop over the G4Hits
  nodename.str("");
  nodename << "G4HIT_" << detector;
  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
  if (hits)
  {
    // this returns an iterator to the beginning and the end of our G4Hits
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)

    {
      // the pointer to the G4Hit is hit_iter->second
      g4hitntuple->Fill(hit_iter->second->get_x(0),
                        hit_iter->second->get_y(0),
                        hit_iter->second->get_z(0),
                        hit_iter->second->get_x(1),
                        hit_iter->second->get_y(1),
                        hit_iter->second->get_z(1),
                        hit_iter->second->get_edep());
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}





int ECCE_DEMP::process_g4clusters(PHCompositeNode* topNode, const string& detector)
{
  ostringstream nodename;

  // loop over the G4Hits
  nodename.str("");
  nodename << "CLUSTER_" << detector;
  RawClusterContainer* clusters = findNode::getClass<RawClusterContainer>(topNode, nodename.str().c_str());
  if (clusters)
  {
    RawClusterContainer::ConstRange cluster_range = clusters->getClusters();
    for (RawClusterContainer::ConstIterator cluster_iter = cluster_range.first; cluster_iter != cluster_range.second; cluster_iter++)
    {
      clusterntuple->Fill(cluster_iter->second->get_phi(),
                          cluster_iter->second->get_z(),
                          cluster_iter->second->get_energy(),
                          cluster_iter->second->getNTowers());
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}


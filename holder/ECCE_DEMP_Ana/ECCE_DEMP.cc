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
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

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

// Track includes
#include <trackbase_historic/SvtxTrackMap.h>

// Jet includes
#include <g4eval/JetEvalStack.h>
#include <g4jets/JetMap.h>

// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

/// HEPMC truth includes
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

/// Fun4All includes
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TH2F.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cassert>
#include <sstream>
#include <string>
#include <iostream>
#include <stdexcept>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

using namespace std;

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

  std::cout << "ECCE_DEMP::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  event_itt = 0;

  gDirectory->mkdir("Pion_Info");
  gDirectory->cd("Pion_Info");
  h1_pi_px = new TH1F("pi_px", "#pi p_{x} Distribution", 200, -20, 20);
  h1_pi_py = new TH1F("pi_py", "#pi p_{y} Distribution", 200, -20, 20);
  h1_pi_pz = new TH1F("pi_pz", "#pi p_{z} Distribution", 200, -50, 50); 
  h1_pi_p = new TH1F("pi_p", "#pi p Distribution", 200, 0, 50);
  h1_pi_E = new TH1F("pi_E", "#pi E Distribution", 200, 0, 50);
  h1_pi_Theta = new TH1F("pi_Theta", "#pi #theta Distribution; #theta [deg]", 200, 0, 50);
  h1_pi_Phi = new TH1F("pi_Phi", "#pi #phi Distribution; #phi [deg]", 360, -180, 180);
  h2_piTrack_ThetaPhi = new TH2F("piTrack_ThetaPhi", "#pi Track #theta vs #phi; #theta [deg]; #phi [deg]", 120, 0, 60, 720, -180, 180);
  h2_piTrack_pTheta = new TH2F("piTrack_pTheta", "#pi Track #theta vs P; #theta [deg]; P [GeV/c]", 120, 0, 60, 500, 0, 50);
  h2_piTrack_ThetaPhi_Smeared = new TH2F("piTrack_ThetaPhi_Smeared", "#pi Track #theta vs #phi; #theta [deg]; #phi [deg]", 120, 0, 60, 720, -180, 180);
  h2_piTrack_pTheta_Smeared = new TH2F("piTrack_pTheta_Smeared", "#pi Track #theta vs P; #theta [deg]; P [GeV/c]", 120, 0, 60, 500, 0, 50);
  gDirectory->cd("../");

  gDirectory->mkdir("Pion_Truth_Info");
  gDirectory->cd("Pion_Truth_Info");
  h1_piTruth_p = new TH1F("piTruth_p", "#pi #frac{#Delta p}{Truth p} Distribution (%); %", 100, -50, 50);
  h1_piTruth_px = new TH1F("piTruth_px", "#pi #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_piTruth_py = new TH1F("piTruth_py", "#pi #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_piTruth_pz = new TH1F("piTruth_pz", "#pi #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_piTruth_E = new TH1F("piTruth_E", "#pi #frac{#Delta E}{Truth E} Distribution (%); %", 100, -50, 50);
  h1_piTruth_p_Smeared = new TH1F("piTruth_p_Smeared", "#pi #frac{#Delta p}{Truth p} Distribution (%); %", 100, -50, 50);
  h1_piTruth_px_Smeared = new TH1F("piTruth_px_Smeared", "#pi #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_piTruth_py_Smeared = new TH1F("piTruth_py_Smeared", "#pi #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_piTruth_pz_Smeared = new TH1F("piTruth_pz_Smeared", "#pi #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_piTruth_E_Smeared = new TH1F("piTruth_E_Smeared", "#pi #frac{#Delta E}{Truth E} Distribution (%); %", 100, -50, 50);
  h2_piTruth_pxpy = new TH2F("piTruth_pxpy", "#pi #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  h2_piTruth_pxpy_Smeared = new TH2F("piTruth_pxpy_Smeared", "#pi #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  gDirectory->cd("../");
  
  gDirectory->mkdir("Scattered_Electron_Info");
  gDirectory->cd("Scattered_Electron_Info");
  h1_e_px = new TH1F("e_px", "e' p_{x} Distribution", 200, -10, 10);
  h1_e_py = new TH1F("e_py", "e' p_{y} Distribution", 200, -10, 10);
  h1_e_pz = new TH1F("e_pz", "e' p_{z} Distribution", 200, -10, 0); 
  h1_e_p = new TH1F("e_p", "e' p Distribution", 200, 0, 10);
  h1_e_E = new TH1F("e_E", "e' E Distribution", 200, 0, 10);
  h1_e_Theta = new TH1F("e_Theta", "e' #theta Distribution; #theta [deg]", 200, 110, 160);
  h1_e_Phi = new TH1F("e_Phi", "e' #phi Distribution; #phi [deg]", 360, -180, 180);
  h2_eTrack_ThetaPhi = new TH2F("eTrack_ThetaPhi", "e' Track #theta vs #phi; #theta [deg]; #phi [deg]", 140, 110, 180, 720, -180, 180);
  h2_eTrack_pTheta = new TH2F("eTrack_pTheta", "e' Track #theta vs P; #theta [deg]; P [GeV/c]", 140, 110, 180, 100, 0, 10);
  h2_eTrack_ThetaPhi_Smeared = new TH2F("eTrack_ThetaPhi_Smeared", "e' Track #theta vs #phi; #theta [deg]; #phi [deg]", 140, 110, 180, 720, -180, 180);
  h2_eTrack_pTheta_Smeared = new TH2F("eTrack_pTheta_Smeared", "e' Track #theta vs P; #theta [deg]; P [GeV/c]", 140, 110, 180, 100, 0, 10);
  gDirectory->cd("../");

  gDirectory->mkdir("Scattered_Electron_Truth_Info");
  gDirectory->cd("Scattered_Electron_Truth_Info");
  h1_eTruth_p = new TH1F("eTruth_p", "e' #frac{#Delta p}{Truth p} Distribution (%) ; %", 100, -50, 50);
  h1_eTruth_px = new TH1F("eTruth_px", "#e' #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_eTruth_py = new TH1F("eTruth_py", "#e' #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_eTruth_pz = new TH1F("eTruth_pz", "e' #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_eTruth_E = new TH1F("eTruth_E", "e' #frac{#Delta E}{Truth E} Distribution (%) ; %", 100, -50, 50);
  h1_eTruth_p_Smeared = new TH1F("eTruth_p_Smeared", "e' #frac{#Delta p}{Truth p} Distribution (%) ; %", 100, -50, 50);
  h1_eTruth_px_Smeared = new TH1F("eTruth_px_Smeared", "#e' #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_eTruth_py_Smeared = new TH1F("eTruth_py_Smeared", "#e' #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_eTruth_pz_Smeared = new TH1F("eTruth_pz_Smeared", "e' #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_eTruth_E_Smeared = new TH1F("eTruth_E_Smeared", "e' #frac{#Delta E}{Truth E} Distribution (%) ; %", 100, -50, 50);
  h2_eTruth_pxpy = new TH2F("eTruth_pxpy", "e' #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);  
  h2_eTruth_pxpy_Smeared = new TH2F("eTruth_pxpy_Smeared", "e' #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  gDirectory->cd("../");

  gDirectory->mkdir("Neutron_Info");
  gDirectory->cd("Neutron_Info");
  h1_n_px = new TH1F("n_px", "n p_{x} Distribution", 320, -4, 4);
  h1_n_py = new TH1F("n_py", "n p_{y} Distribution", 200, -2.5, 2.5);
  h1_n_pz = new TH1F("n_pz", "n p_{z} Distribution", 240, 0, 120); 
  h1_n_p = new TH1F("n_p", "n p Distribution", 240, 0, 120);
  h1_n_E = new TH1F("n_E", "n E Distribution", 240, 0, 120);
  h1_n_Theta = new TH1F("n_Theta", "n #theta Distribution; #theta [deg]", 300, 0, 3);
  h1_n_Phi = new TH1F("n_Phi", "n #phi Distribution; #phi [deg]", 400, -20, 20);
  h2_nTrack_ThetaPhi = new TH2F("nTrack_ThetaPhi", "n Track #theta vs #phi; #theta [deg]; #phi [deg]", 100, 0, 5, 100, -50, 50);
  h2_nTrack_pTheta = new TH2F("nTrack_pTheta", "n Track #theta vs P; #theta [deg]; P [GeV/c]", 100, 0, 5, 1000, 0, 100);
  h2_nTrack_ThetaPhi_Smeared = new TH2F("nTrack_ThetaPhi_Smeared", "n Track #theta vs #phi; #theta [deg]; #phi [deg]", 100, 0, 5, 100, -50, 50);
  h2_nTrack_pTheta_Smeared = new TH2F("nTrack_pTheta_Smeared", "n Track #theta vs P; #theta [deg]; P [GeV/c]", 100, 0, 5, 1000, 0, 100);
  gDirectory->cd("../");

  gDirectory->mkdir("Neutron_Truth_Info");
  gDirectory->cd("Neutron_Truth_Info");
  h1_nTruth_p = new TH1F("nTruth_p", "n #frac{#Delta p}{Truth p} Distribution (%) ; %", 100, -50, 50);
  h1_nTruth_px = new TH1F("nTruth_px", "#n #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_nTruth_py = new TH1F("nTruth_py", "#n #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_nTruth_pz = new TH1F("nTruth_pz", "n #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_nTruth_E = new TH1F("nTruth_E", "n #frac{#Delta E}{Truth E} Distribution (%) ; %", 100, -50, 50);
  h1_nTruth_p_Smeared = new TH1F("nTruth_p_Smeared", "n #frac{#Delta p}{Truth p} Distribution (%) ; %", 100, -50, 50);
  h1_nTruth_px_Smeared = new TH1F("nTruth_px_Smeared", "#n #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_nTruth_py_Smeared = new TH1F("nTruth_py_Smeared", "#n #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_nTruth_pz_Smeared = new TH1F("nTruth_pz_Smeared", "n #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_nTruth_E_Smeared = new TH1F("nTruth_E_Smeared", "n #frac{#Delta E}{Truth E} Distribution (%) ; %", 100, -50, 50);
  h2_nTruth_pxpy = new TH2F("nTruth_pxpy", "n #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  h2_nTruth_pxpy_Smeared = new TH2F("nTruth_pxpy_Smeared", "n #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  gDirectory->cd("../");

  /*goDirectory->mkdir("Neutron_lowE_Info"); // added by Maggie Kerr; July 28, 2021
  goDirectory->cd("Neutron_lowE_Info");
  h1_nlowETruth_E = new TH1F("nlowETruth_E","n #frac{#Delta E}{Truth E} Distribution (%) for < 70 GeV ; %", 100, -50, 50);
  h1_nlowETruth_p = new TH1F("nlowETruth_P","n #frac{#Delta P}{Truth P} Distribution (%) for < 70 GeV ; %", 100, -50, 50);
  h1_nlowETruth_px = new TH1F("nlowETruth_px","n #frac{#Delta px}{Truth py} Distribution (%) for < 70 GeV ; %", 100, -50, 50);
  h1_nlowETruth_py = new TH1F("nlowETruth_py_","n #frac{#Delta py}{Truth py} Distribution (%) for < 70 GeV ; %", 100, -50, 50);
  h1_nlowETruth_pz = new TH1F("nlowETruth_pz_","n #frac{#Delta pz}{Truth pz} Distribution (%) for < 70 GeV ; %", 100, -50, 50);
  h1_nlowETruth_theta = new TH1F("nlowETruth_theta","n #frac{#Delta #theta}{Truth #theta} Distribution (%) for < 70 GeV ; %", 100, -50, 50);
  h1_nlowETruth_phi = new TH1F("nlowETruth_phi","n #frac{#Delta #phi}{Truth #phi} Distribution (%) for < 70 GeV ; %", 100, -50, 50);
  h2_nlowETruth_pxpy = new TH2F("nlowETruth_pxpy","n #frac{#Delta px}{Truth px} vs #frac{#Delta py}{Truth py} Distribution (%) for < 70 GeV ; %", 100, -50, 50, 100, -50, 50);
  h1_nlowETruth_E_Smeared = new TH1F("nlowETruth_E_Smeared","n #frac{#Delta E}{Truth E} Distribution (%) for < 70 GeV (Smeared); %", 100, -50, 50);
  h1_nlowETruth_p_Smeared = new TH1F("nlowETruth_P_Smeared","n #frac{#Delta P}{Truth P} Distribution (%) for < 70 GeV (Smeared); %", 100, -50, 50);
  h1_nlowETruth_px_Smeared = new TH1F("nlowETruth_px_Smeared","n #frac{#Delta px}{Truth px} Distribution (%) for < 70 GeV (Smeared); %", 100, -50, 50);
  h1_nlowETruth_py_Smeared = new TH1F("nlowETruth_py_Smeared","n #frac{#Delta py}{Truth py} Distribution (%) for < 70 GeV (Smeared); %", 100, -50, 50);
  h1_nlowETruth_pz_Smeared = new TH1F("nlowETruth_pz_Smeared","n #frac{#Delta pz}{Truth pz} Distribution (%) for < 70 GeV (Smeared); %", 100, -50, 50);
  h1_nlowETruth_theta_Smeared = new TH1F("nlowETruth_theta_Smeared","n #frac{#Delta #theta}{Truth #theta} Distribution (%) for < 70 GeV (Smeared); %", 100, -50, 50);
  h1_nlowETruth_phi_Smeared = new TH1F("nlowETruth_phi_Smeared","n #frac{#Delta #phi}{Truth #phi} Distribution (%) for < 70 GeV (Smeared); %", 100, -50, 50);
  h2_nlowETruth_pxpy_Smeared = new TH2F("nlowETruth_pxpy_Smeared","n #frac{#Delta px}{Truth px} vs #frac{#Delta py}{Truth py} Distribution (%) for < 70 GeV (Smeared); %", 100, -50, 50, 100, -50, 50);
  goDirectory->cd("../");

  //goDirectory->mkdir("Pion_NoTrack_Info");
  //goDirectory->cd("Pion_NoTrack_Info");
  //goDirectory->cd("../");

  //goDirectory->mkdir("Electron_NoTrack_Info");
  //goDirectory->cd("Electron_NoTrack_Info");
  //goDirectory->cd("../");*/

  gDirectory->mkdir("PMiss_Info");
  gDirectory->cd("PMiss_Info");
  h1_pmiss_px = new TH1F("pmiss_px", "p_{miss} p_{x} Distribution", 800, -10, 10);
  h1_pmiss_py = new TH1F("pmiss_py", "p_{miss} p_{y} Distribution", 200, -2.5, 2.5);
  h1_pmiss_pz = new TH1F("pmiss_pz", "p_{miss} p_{z} Distribution", 240, 0, 120); 
  h1_pmiss_p = new TH1F("pmiss_p", "p_{miss} p Distribution", 240, 0, 120);
  h1_pmiss_E = new TH1F("pmiss_E", "p_{miss} E Distribution", 240, 0, 120);
  h1_pmiss_Theta = new TH1F("pmiss_Theta", "p_{miss} #theta Distribution; #theta [deg]", 1000, 0, 10);
  h1_pmiss_Phi = new TH1F("pmiss_Phi", "p_{miss} #phi Distribution; #phi [deg]", 720, -180, 180);
  gDirectory->cd("../");
  
  gDirectory->mkdir("Virtual_Photon_Info");
  gDirectory->cd("Virtual_Photon_Info");
  h1_gamma_px = new TH1F("gamma_px", "#gamma p_{x} Distribution", 200, -10, 10);
  h1_gamma_py = new TH1F("gamma_py", "#gamma p_{y} Distribution", 200, -10, 10);
  h1_gamma_pz = new TH1F("gamma_pz", "#gamma p_{z} Distribution", 200, -10, 0); 
  h1_gamma_p = new TH1F("gamma_p", "#gamma p Distribution", 200, 0, 10);
  h1_gamma_E = new TH1F("gamma_E", "#gamma E Distribution", 200, 0, 10);
  h1_gamma_Theta = new TH1F("gamma_Theta", "#gamma #theta Distribution; #theta [deg]", 360, -180, 180);
  h1_gamma_Phi = new TH1F("gamma_Phi", "#gamma #phi Distribution; #phi [deg]", 360, -180, 180);
  gDirectory->cd("../");
 
  gDirectory->mkdir("Kinematics_Info");
  gDirectory->cd("Kinematics_Info");
  h1_Q2_Dist = new TH1F("Q2_Dist", "Q^{2} Distribution", 200, 0, 50);
  h1_W_Dist = new TH1F("W_Dist", "W Distribution", 500, 0, 50);
  h1_t_Dist = new TH1F("t_Dist", "t Distribution", 100, 0, 10);
  h1_xb_Dist = new TH1F("xb_Dist", "x_{b} Distribution", 100, 0, 1);
  h1_xi_Dist = new TH1F("xi_Dist", "#xi Distribution", 100, 0, 1);
  gDirectory->cd("../");

  h2_ZDC_XY = new TH2F("ZDC_XY", "ZDC XY", 200, -50, 50, 200, -50, 50);
  h2_ZDC_XY_Smeared = new TH2F("ZDC_XY_Smeared", "ZDC XY", 200, -50, 50, 200, -50, 50);

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
  
  ZDC_hit = 0;
  EEMC_hit = 0;
  Double_t Pi = TMath::ACos(-1);
  eBeam4Vect.SetPxPyPzE(0,0,-5,5);
  pBeam4Vect.SetPxPyPzE(-100*TMath::Sin(0.05),100*TMath::Sin(0.05)*TMath::Sin(Pi),100*TMath::Cos(0.05),100);

  event_itt++; 
 
  if(event_itt%100 == 0)
     std::cout << "Event Processing Counter: " << event_itt << endl;

  //process_g4hits_ZDC(topNode);

  if (Check_n(topNode) == true && Check_ePi(topNode) == true){ // For event, check if it look like we have an e/pi/n in the event
    // Get track map for e'/pi info
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    // Get ZDC hits for neutron info
    PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_ZDC");
    // Get MC truth info
    PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    // Get the primary particle range
    PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();

    if (!trackmap)
      {
	trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
	if (!trackmap)
	  {
	    cout
	      << "ECCE_DEMP::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
	    exit(-1);
	  }
      }

    // Loop over our tracks, assign info to e'/pi
    for (SvtxTrackMap::Iter iter = trackmap->begin();
	 iter != trackmap->end();
	 ++iter)
      {
	SvtxTrack* track = iter->second;
	if ( track->get_pz() > 0 && track->get_charge() == 1){ // +ve z direction -> pions, crappy way of selecting them for now w/o truth info
	  piVect.SetXYZ(track->get_px(), track->get_py(), track->get_pz());
	  piVectSmeared.SetXYZ(Position_Smear(track->get_px()), Position_Smear(track->get_py()), Position_Smear(track->get_pz()));
	  pi4Vect.SetPxPyPzE(track->get_px(), track->get_py(), track->get_pz(), sqrt(pow(piVect.Mag(), 2)+pow(mPi,2)));
	  pi4VectSmeared.SetPxPyPzE(piVectSmeared.x(), piVectSmeared.y(), piVectSmeared.z(), sqrt(pow(piVect.Mag(), 2)+pow(mPi,2)));
	}

	else if (track->get_pz() < 0  && track->get_charge() == -1 ){ // -ve z direction -> electrons, crappy way of selecting them for now w/o truth info
	  eVect.SetXYZ(track->get_px(), track->get_py(), track->get_pz());
	  eVectSmeared.SetXYZ(Position_Smear(track->get_px()), Position_Smear(track->get_py()), Position_Smear(track->get_pz()));
	  e4Vect.SetPxPyPzE(track->get_px(), track->get_py(), track->get_pz(), sqrt(pow(eVect.Mag(), 2)+pow(mElec,2)));
	  e4VectSmeared.SetPxPyPzE(eVectSmeared.x(), eVectSmeared.y(), eVectSmeared.z(), sqrt(pow(eVect.Mag(), 2)+pow(mElec,2)));
	}
      }
    
    // Loop over the hts in the zdc
    if (hits) {
      // this returns an iterator to the beginning and the end of our G4Hits
      PHG4HitContainer::ConstRange hit_range = hits->getHits();
      for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
	{
	  nZDCPos.SetXYZ(hit_iter->second->get_x(0), hit_iter->second->get_y(0), hit_iter->second->get_z(0));
	  nZDCPosSmeared.SetXYZ(Position_Smear(hit_iter->second->get_x(0)), Position_Smear(hit_iter->second->get_y(0)), Position_Smear(hit_iter->second->get_z(0)));
	  nEDep = hit_iter->second->get_edep();
	  nEDepSmeared = EMCAL_Smear(hit_iter->second->get_edep());
	  nTheta = nZDCPos.Theta();
	  nThetaSmeared = nZDCPosSmeared.Theta();
	  nPhi = nZDCPos.Phi();
	  nPhiSmeared = nZDCPosSmeared.Phi();
	  nPMag = sqrt((pow(nEDep,2)) - (pow(mNeut,2)));
	  nPMagSmeared = sqrt((pow(nEDepSmeared,2)) - (pow(mNeut,2)));
	  n4Vect.SetPxPyPzE(nPMag*sin(nTheta)*cos(nPhi), nPMag*sin(nTheta)*sin(nPhi), nPMag*cos(nTheta), nEDep);	  
	  n4VectSmeared.SetPxPyPzE(nPMagSmeared*sin(nThetaSmeared)*cos(nPhiSmeared), nPMagSmeared*sin(nThetaSmeared)*sin(nPhiSmeared), nPMagSmeared*cos(nThetaSmeared), nEDepSmeared);
	}
    }

    if (!truthinfo)
      {
	cout << PHWHERE
	     << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"
	     << endl;
	return Fun4AllReturnCodes::EVENT_OK;
      }

    /// Loop over the G4 truth (stable) particles
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
	 iter != range.second;
	 ++iter)
      {
	/// Get this truth particle
	const PHG4Particle *truth = iter->second;
	if ( truth->get_pid() == 11){ // PDG 11 -> Scattered electron
	  e4VectTruth.SetPxPyPzE(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
	}
	else if (truth->get_pid() == 211){ // PDG 211 -> Pion 
	  pi4VectTruth.SetPxPyPzE(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
	}
	else if (truth->get_pid() == 2112){ // PDG 2112 -> Neutron
	  n4VectTruth.SetPxPyPzE(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
	}
      }

    
    // Now have relevant information from this event, fill some histograms and calculate some stuff

    // Calculate kinematic quantities
    virtphoton4Vect = eBeam4Vect - e4Vect;
    t4Vect = virtphoton4Vect - pi4Vect;
    pmiss4Vect = (eBeam4Vect + pBeam4Vect) - (e4Vect+pi4Vect);
    Q2 = -1*(virtphoton4Vect.Mag2());
    W = (virtphoton4Vect+pBeam4Vect).Mag();
    t = -(t4Vect.Mag2());
    xb =  Q2/(2*(pBeam4Vect.Dot(virtphoton4Vect)));
    xi = xb/(2-xb);

    // Truth versions of kinematic quantities
    virtphoton4VectTruth = eBeam4Vect - e4VectTruth;
    t4VectTruth = virtphoton4VectTruth - pi4VectTruth;
    pmiss4VectTruth = (eBeam4Vect + pBeam4Vect) - (e4VectTruth+pi4VectTruth);
    Q2_truth = -1*(virtphoton4VectTruth.Mag2());
    W_truth = (virtphoton4VectTruth+pBeam4Vect).Mag();
    t_truth = -(t4VectTruth.Mag2());
    xb_truth =  Q2_truth/(2*(pBeam4Vect.Dot(virtphoton4VectTruth)));
    xi_truth = xb_truth/(2-xb_truth);
      
    // Fill histograms
    h1_pi_px->Fill(pi4Vect.Px());
    h1_pi_py->Fill(pi4Vect.Py());
    h1_pi_pz->Fill(pi4Vect.Pz());
    h1_pi_p->Fill(pi4Vect.P());
    h1_pi_E->Fill(pi4Vect.E());
    h1_pi_Theta->Fill(pi4Vect.Theta()*TMath::RadToDeg());
    h1_pi_Phi->Fill(pi4Vect.Phi()*TMath::RadToDeg());
    h1_e_px->Fill(e4Vect.Px());
    h1_e_py->Fill(e4Vect.Py());
    h1_e_pz->Fill(e4Vect.Pz());
    h1_e_p->Fill(e4Vect.P());
    h1_e_E->Fill(e4Vect.E());
    h1_e_Theta->Fill(e4Vect.Theta()*TMath::RadToDeg());
    h1_e_Phi->Fill(e4Vect.Phi()*TMath::RadToDeg());
    h1_n_px->Fill(n4Vect.Px());
    h1_n_py->Fill(n4Vect.Py());
    h1_n_pz->Fill(n4Vect.Pz());
    h1_n_p->Fill(n4Vect.P());
    h1_n_E->Fill(n4Vect.E());
    h1_n_Theta->Fill(n4Vect.Theta()*TMath::RadToDeg());
    h1_n_Phi->Fill(n4Vect.Phi()*TMath::RadToDeg());
    h1_pmiss_px->Fill(pmiss4Vect.Px());
    h1_pmiss_py->Fill(pmiss4Vect.Py());
    h1_pmiss_pz->Fill(pmiss4Vect.Pz());
    h1_pmiss_p->Fill(pmiss4Vect.P());
    h1_pmiss_E->Fill(pmiss4Vect.E());
    h1_pmiss_Theta->Fill(pmiss4Vect.Theta()*TMath::RadToDeg());
    h1_pmiss_Phi->Fill(pmiss4Vect.Phi()*TMath::RadToDeg());
    h1_gamma_px->Fill(virtphoton4Vect.Px());
    h1_gamma_py->Fill(virtphoton4Vect.Py());
    h1_gamma_pz->Fill(virtphoton4Vect.Pz());
    h1_gamma_p->Fill(virtphoton4Vect.P());
    h1_gamma_E->Fill(virtphoton4Vect.E());
    h1_gamma_Theta->Fill(virtphoton4Vect.Theta()*TMath::RadToDeg());
    h1_gamma_Phi->Fill(virtphoton4Vect.Phi()*TMath::RadToDeg());

    h1_Q2_Dist->Fill(Q2);
    h1_W_Dist->Fill(W);
    h1_t_Dist->Fill(t);
    h1_xb_Dist->Fill(xb);
    h1_xi_Dist->Fill(xi);

    h1_piTruth_p->Fill((pi4Vect.P()-pi4VectTruth.P())/(pi4VectTruth.P())*100);
    h1_piTruth_px->Fill((pi4Vect.Px()-pi4VectTruth.Px())/(pi4VectTruth.Px())*100);
    h1_piTruth_py->Fill((pi4Vect.Py()-pi4VectTruth.Py())/(pi4VectTruth.Py())*100);
    h1_piTruth_pz->Fill((pi4Vect.Pz()-pi4VectTruth.Pz())/(pi4VectTruth.Pz())*100);
    h1_piTruth_E->Fill((pi4Vect.E()-pi4VectTruth.E())/(pi4VectTruth.E())*100);
    h1_eTruth_p->Fill((e4Vect.P()-e4VectTruth.P())/(e4VectTruth.P())*100);
    h1_eTruth_px->Fill((e4Vect.Px()-e4VectTruth.Px())/(e4VectTruth.Px())*100);
    h1_eTruth_py->Fill((e4Vect.Py()-e4VectTruth.Py())/(e4VectTruth.Py())*100);
    h1_eTruth_pz->Fill((e4Vect.Pz()-e4VectTruth.Pz())/(e4VectTruth.Pz())*100);
    h1_eTruth_E->Fill((e4Vect.E()-e4VectTruth.E())/(e4VectTruth.E())*100);
    h1_nTruth_p->Fill((n4Vect.P()-n4VectTruth.P())/(n4VectTruth.P())*100);
    h1_nTruth_px->Fill((n4Vect.Px()-n4VectTruth.Px())/(n4VectTruth.Px())*100);
    h1_nTruth_py->Fill((n4Vect.Py()-n4VectTruth.Py())/(n4VectTruth.Py())*100);
    h1_nTruth_pz->Fill((n4Vect.Pz()-n4VectTruth.Pz())/(n4VectTruth.Pz())*100);
    h1_nTruth_E->Fill((n4Vect.E()-n4VectTruth.E())/(n4VectTruth.E())*100);

    h1_piTruth_p_Smeared->Fill((pi4VectSmeared.P()-pi4VectTruth.P())/(pi4VectTruth.P())*100);
    h1_piTruth_px_Smeared->Fill((pi4VectSmeared.Px()-pi4VectTruth.Px())/(pi4VectTruth.Px())*100);
    h1_piTruth_py_Smeared->Fill((pi4VectSmeared.Py()-pi4VectTruth.Py())/(pi4VectTruth.Py())*100);
    h1_piTruth_pz_Smeared->Fill((pi4VectSmeared.Pz()-pi4VectTruth.Pz())/(pi4VectTruth.Pz())*100);
    h1_piTruth_E_Smeared->Fill((pi4VectSmeared.E()-pi4VectTruth.E())/(pi4VectTruth.E())*100);
    h1_eTruth_p_Smeared->Fill((e4VectSmeared.P()-e4VectTruth.P())/(e4VectTruth.P())*100);
    h1_eTruth_px_Smeared->Fill((e4VectSmeared.Px()-e4VectTruth.Px())/(e4VectTruth.Px())*100);
    h1_eTruth_py_Smeared->Fill((e4VectSmeared.Py()-e4VectTruth.Py())/(e4VectTruth.Py())*100);
    h1_eTruth_pz_Smeared->Fill((e4VectSmeared.Pz()-e4VectTruth.Pz())/(e4VectTruth.Pz())*100);
    h1_eTruth_E_Smeared->Fill((e4VectSmeared.E()-e4VectTruth.E())/(e4VectTruth.E())*100);
    h1_nTruth_p_Smeared->Fill((n4VectSmeared.P()-n4VectTruth.P())/(n4VectTruth.P())*100);
    h1_nTruth_px_Smeared->Fill((n4VectSmeared.Px()-n4VectTruth.Px())/(n4VectTruth.Px())*100);
    h1_nTruth_py_Smeared->Fill((n4VectSmeared.Py()-n4VectTruth.Py())/(n4VectTruth.Py())*100);
    h1_nTruth_pz_Smeared->Fill((n4VectSmeared.Pz()-n4VectTruth.Pz())/(n4VectTruth.Pz())*100);
    h1_nTruth_E_Smeared->Fill((n4VectSmeared.E()-n4VectTruth.E())/(n4VectTruth.E())*100);
    
    h2_ZDC_XY->Fill(nZDCPos.x()-90, nZDCPos.y());
    h2_ZDC_XY_Smeared->Fill(nZDCPosSmeared.x()-90, nZDCPosSmeared.y());

    h2_piTrack_ThetaPhi->Fill((pi4Vect.Theta()*TMath::RadToDeg()), (pi4Vect.Phi()*TMath::RadToDeg()));
    h2_piTrack_pTheta->Fill((pi4Vect.Theta()*TMath::RadToDeg()), pi4Vect.P());
    h2_eTrack_ThetaPhi->Fill((e4Vect.Theta()*TMath::RadToDeg()), (e4Vect.Phi()*TMath::RadToDeg()));
    h2_eTrack_pTheta->Fill((e4Vect.Theta()*TMath::RadToDeg()), e4Vect.P());
    h2_nTrack_ThetaPhi->Fill((n4Vect.Theta()*TMath::RadToDeg()), (n4Vect.Phi()*TMath::RadToDeg()));
    h2_nTrack_pTheta->Fill((n4Vect.Theta()*TMath::RadToDeg()), n4Vect.P());
    h2_piTrack_ThetaPhi_Smeared->Fill((pi4VectSmeared.Theta()*TMath::RadToDeg()), (pi4VectSmeared.Phi()*TMath::RadToDeg()));
    h2_piTrack_pTheta_Smeared->Fill((pi4VectSmeared.Theta()*TMath::RadToDeg()), pi4VectSmeared.P());
    h2_eTrack_ThetaPhi_Smeared->Fill((e4VectSmeared.Theta()*TMath::RadToDeg()), (e4VectSmeared.Phi()*TMath::RadToDeg()));
    h2_eTrack_pTheta_Smeared->Fill((e4VectSmeared.Theta()*TMath::RadToDeg()), e4VectSmeared.P());
    h2_nTrack_ThetaPhi_Smeared->Fill((n4VectSmeared.Theta()*TMath::RadToDeg()), (n4VectSmeared.Phi()*TMath::RadToDeg()));
    h2_nTrack_pTheta_Smeared->Fill((n4VectSmeared.Theta()*TMath::RadToDeg()), n4VectSmeared.P());
    
    h2_piTruth_pxpy->Fill((pi4Vect.Px()-pi4VectTruth.Px())/(pi4VectTruth.Px())*100, (pi4Vect.Py()-pi4VectTruth.Py())/(pi4VectTruth.Py())*100);
    h2_eTruth_pxpy->Fill((e4Vect.Px()-e4VectTruth.Px())/(e4VectTruth.Px())*100, (e4Vect.Py()-e4VectTruth.Py())/(e4VectTruth.Py())*100);
    h2_nTruth_pxpy->Fill((n4Vect.Px()-n4VectTruth.Px())/(n4VectTruth.Px())*100, (n4Vect.Py()-n4VectTruth.Py())/(n4VectTruth.Py())*100);
    h2_piTruth_pxpy_Smeared->Fill((pi4VectSmeared.Px()-pi4VectTruth.Px())/(pi4VectTruth.Px())*100, (pi4VectSmeared.Py()-pi4VectTruth.Py())/(pi4VectTruth.Py())*100);
    h2_eTruth_pxpy_Smeared->Fill((e4VectSmeared.Px()-e4VectTruth.Px())/(e4VectTruth.Px())*100, (e4VectSmeared.Py()-e4VectTruth.Py())/(e4VectTruth.Py())*100);
    h2_nTruth_pxpy_Smeared->Fill((n4VectSmeared.Px()-n4VectTruth.Px())/(n4VectTruth.Px())*100, (n4VectSmeared.Py()-n4VectTruth.Py())/(n4VectTruth.Py())*100);

    /*if (nEDep < 70) {
    h1_nlowETruth_E->Fill(((n4Vect.E()-n4VectTruth.E())/(n4VectTruth.E()))*100);
    h1_nlowETruth_p->Fill(((n4Vect.P()-n4VectTruth.P())/(n4VectTruth.P()))*100);
    h1_nlowETruth_px->Fill(((n4Vect.Px()-n4VectTruth.Px())/(n4VectTruth.Px()))*100);
    h1_nlowETruth_py->Fill(((n4Vect.Py()-n4VectTruth.Py())/(n4VectTruth.Py()))*100);
    h1_nlowETruth_pz->Fill(((n4Vect.Pz()-n4VectTruth.Pz())/(n4VectTruth.Pz()))*100);
    h1_nlowETruth_theta->Fill((((n4Vect.Theta()*TMath::RadToDeg())-(n4VectTruth.Theta()*TMath::RadToDeg()))/((n4VectTruth.Theta()*TMath::RadToDeg())))*100);
    h1_nlowETruth_phi->Fill((((n4Vect.Phi()*TMath::RadToDeg())-(n4VectTruth.Phi()*TMath::RadToDeg()))/((n4VectTruth.Phi()*TMath::RadToDeg())))*100);
    h2_nlowETruth_pxpy->Fill(((n4Vect.Px()-n4VectTruth.Px())/(n4VectTruth.Px()))*100,((n4Vect.Py()-n4VectTruth.Py())/(n4VectTruth.Py()))*100);
    h1_nlowETruth_E_Smeared->Fill(((n4VectSmeared.E()-n4VectTruth.E())/(n4VectTruth.E()))*100);
    h1_nlowETruth_p_Smeared->Fill(((n4VectSmeared.P()-n4VectTruth.P())/(n4VectTruth.P()))*100);
    h1_nlowETruth_px_Smeared->Fill(((n4VectSmeared.Px()-n4VectTruth.Px())/(n4VectTruth.Px()))*100);
    h1_nlowETruth_py_Smeared->Fill(((n4VectSmeared.Py()-n4VectTruth.Py())/(n4VectTruth.Py()))*100);
    h1_nlowETruth_pz_Smeared->Fill(((n4VectSmeared.Pz()-n4VectTruth.Pz())/(n4VectTruth.Pz()))*100);
    h1_nlowETruth_theta_Smeared->Fill((((n4VectSmeared.Theta()*TMath::RadToDeg())-(n4VectTruth.Theta()*TMath::RadToDeg()))/((n4VectTruth.Theta()*TMath::RadToDeg())))*100);
    h1_nlowETruth_phi_Smeared->Fill((((n4VectSmeared.Phi()*TMath::RadToDeg())-(n4VectTruth.Phi()*TMath::RadToDeg()))/((n4VectTruth.Phi()*TMath::RadToDeg())))*100);
    h2_nlowETruth_pxpy_Smeared->Fill(((n4VectSmeared.Px()-n4VectTruth.Px())/(n4VectTruth.Px()))*100,((n4VectSmeared.Py()-n4VectTruth.Py())/(n4VectTruth.Py()))*100);
    }*/

  }

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

  outfile->cd();
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

//***************************************************

bool ECCE_DEMP::Check_ePi(PHCompositeNode* topNode)
{
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!trackmap)
    {
      trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
      if (!trackmap)
    	{
    	  cout
    	    << "ECCE_DEMP::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
    	  exit(-1);
    	}
    }
  int nTracks = 0;
  Bool_t ElecTrack = kFALSE;
  Bool_t PionTrack = kFALSE;
  // Iterate over tracks
  for (SvtxTrackMap::Iter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter)
    {
      SvtxTrack* track = iter->second;
      nTracks++;
      if ( track->get_pz() > 0 && track->get_charge() == 1){ // +ve z direction -> pions, crappy way of selecting them for now w/o truth info
	PionTrack = kTRUE;
      }
      else if (track->get_pz() < 0  && track->get_charge() == -1 ){ // -ve z direction -> electrons, crappy way of selecting them for now w/o truth info
	ElecTrack = kTRUE;
      }
    }
  
  if( PionTrack == kTRUE && ElecTrack == kTRUE && nTracks == 2){ // Both a pion and an electron track, only 2 tracks
    return true;
  }
  else{
    return false;
  }
}


//***************************************************

bool ECCE_DEMP::Check_n(PHCompositeNode* topNode)
{
  // loop over the G4Hits

  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_ZDC");

  int ZDCHits = 0;
  Bool_t nZDChit = kFALSE;

  if (hits) {
    // this returns an iterator to the beginning and the end of our G4Hits
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
      {
	ZDCHits++;
	if (hit_iter->second->get_edep() > 40){ // Hit in ZDC with roughly correct energy for neutron
	  nZDChit = kTRUE;
	}
      }
  }

  if( nZDChit == kTRUE ){ // Hit in ZDC with correct energy
    return true;
  }
  else if (nZDChit == kFALSE){
    return false;
  }
  else return false;
}

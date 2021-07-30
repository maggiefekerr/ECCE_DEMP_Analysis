// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ECCE_DEMP_ANA_H
#define ECCE_DEMP_ANA_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <memory>
#include <string>
#include <utility>  // std::pair, std::make_pair                                                                                                                                                            
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"


class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class JetEvalStack;

class ECCE_DEMP : public SubsysReco
{
 public:

  ECCE_DEMP(const std::string &name = "ECCE_DEMP", const std::string &fname = "MyNtuple.root");

  virtual ~ECCE_DEMP();

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  int process_g4hits_ZDC(PHCompositeNode *);

  int process_g4hits(PHCompositeNode *, const std::string&);
  int process_g4clusters(PHCompositeNode *, const std::string&);

  int process_g4tracks(PHCompositeNode *);

  bool Check_ePi(PHCompositeNode *);
  bool Check_n(PHCompositeNode *);

  void use_initial_vertex(const bool b = true) {initial_vertex = b;}

  //private:

 protected:

  //! flag to use initial vertex in track evaluator 
  bool initial_vertex = false;

  std::string detector;
  std::string outfilename;
  Fun4AllHistoManager *hm;

  TFile *outfile;
  TNtuple *g4hitntuple;
  TNtuple *g4trackntuple;
  TNtuple *ZDChitntuple;
  TNtuple *EEMChitntuple;
  TNtuple *EEMCalclusterntuple;
  TNtuple *truthtree; // added by Maggie Kerr
  unsigned long long int event_itt;
  gsl_rng* m_RandomGenerator;

  //*********************************
  // Energy and Position smearing

  float EMCAL_Smear(float E);
  float HCAL_Smear(float E);
  float PbWO4_Smear(float E);
  float Position_Smear(float E);

  //---------------------
  // From ejana

  double true_q2;
  double true_x;
  double true_s_e;
  double true_xpi;
  double true_ypi;
  double true_tpi;

  double have_true_dis_info = false;
  
  bool  HIT_IN_ZDC;
  bool  HIT_IN_EMCAL;
  bool  HIT_IN_HCAL;
  bool  CLUS_IN_EMCAL;
  bool  CLUS_IN_HCAL;
  bool  HIT_IN_HEC;	

  double e_beam_energy;
  double ion_beam_energy;

  double crossing_angle;

  TLorentzVector r_lelectron;
//  TLorentzVector r_lproton;

  TLorentzVector r_lscatelec;
  TLorentzVector r_l_scat_nucleon;

  TLorentzVector lproton;

  // Particle Masses
  Double_t mPi = 0.139570;
  Double_t mElec = 0.000510998950;
  Double_t mNeut = 0.93965420;

  // Quantities we want to determine
  TVector3 eVect;
  TVector3 piVect;
  TVector3 nZDCPos;
  TLorentzVector e4Vect;
  TLorentzVector pi4Vect;
  TLorentzVector n4Vect;
  TLorentzVector e4VectTruth;
  TLorentzVector pi4VectTruth;
  TLorentzVector n4VectTruth;
  Double_t nEDep;
  Double_t nTheta;
  Double_t nPhi;
  Double_t nPMag;

  Int_t ZDC_hit;
  Int_t EEMC_hit;

  TH2F* h2_ZDC_XY; 
  TH2F* h2_ZDC_XY_nEnergy;
  TH2F* h2_ZDC_XY_nEnergy_Smeared;

  TH1F* h1_ZDC_E_dep;
  TH1F* h1_ZDC_E_dep_smeared;

  TH1F* h1_EMCal_E_dep;
  TH1F* h1_EMCal_E_dep_smeared;

  TH1F* h1_HCalIn_E_dep;
  TH1F* h1_HCalIn_E_dep_smeared;

  TH1F* h1_HCalOut_E_dep;
  TH1F* h1_HCalOut_E_dep_smeared;

  TH1F* h1_eTrack_px;
  TH1F* h1_eTrack_py;
  TH1F* h1_eTrack_pz;
  TH1F* h1_eTrack_p;
  TH1F* h1_eTrack_theta;
  TH1F* h1_eTrack_phi;

  TH1F* h1_piTrack_px;
  TH1F* h1_piTrack_py;
  TH1F* h1_piTrack_pz;
  TH1F* h1_piTrack_p;
  TH1F* h1_piTrack_theta;
  TH1F* h1_piTrack_phi;

  TH1F* h1_nTracksDist;

  TH2F* h2_ePiTrackDist;

  TH2F* h2_eTrack_ThetaPhi;
  TH2F* h2_eTrack_pTheta;

  TH2F* h2_piTrack_ThetaPhi;
  TH2F* h2_piTrack_pTheta;

  // Histogram for coincidence analysis routine

  // TH1F* h1_piTruth_p;
  // TH1F* h1_piTruth_px;
  // TH1F* h1_piTruth_py;
  // TH1F* h1_piTruth_pz;
  // TH1F* h1_piTruth_E;
  // TH1F* h1_eTruth_p;
  // TH1F* h1_eTruth_px;
  // TH1F* h1_eTruth_py;
  // TH1F* h1_eTruth_pz;
  // TH1F* h1_eTruth_E;
  // TH1F* h1_nTruth_p;
  // TH1F* h1_nTruth_px;
  // TH1F* h1_nTruth_py;
  // TH1F* h1_nTruth_pz;
  // TH1F* h1_nTruth_E;

// added by Maggie Kerr

// for hits in ZDC with < 70 GeV of energy, _per converts to a percent comparison
  TH2F* h2_ZDC_TvDep_E;
  // TH1F* h1_ZDC_TvDep_E_per;
  TH2F* h2_ZDC_TvDep_p;
  // TH1F* h1_ZDC_TvDep_p_per;
  TH2F* h2_ZDC_TvDep_px;
  // TH1F* h1_ZDC_TvDep_px_per;
  TH2F* h2_ZDC_TvDep_py;
  // TH1F* h1_ZDC_TvDep_py_per;
  TH2F* h2_ZDC_TvDep_pz;
  // TH1F* h1_ZDC_TvDep_pz_per;
  TH2F* h2_ZDC_TvDep_theta;
  // TH1F* h1_ZDC_TvDep_theta_per;
  TH2F* h2_ZDC_TvDep_phi;
  // TH1F* h1_ZDC_TvDep_phi_per;
  TH2F* h2_ZDC_xy_lt70; // records position of hit

// for e/pi events without track information, _per converts to a percent comparison
  // TH1F* h1_truth_nohits_eE;
  // TH1F* h1_truth_nohits_epx;
  // TH1F* h1_truth_nohits_epy;
  // TH1F* h1_truth_nohits_epz;
  // TH1F* h1_truth_nohits_etheta;
  // TH1F* h1_truth_nohits_ephi;

  // TH1F* h1_truth_nohits_eE_per;
  // TH1F* h1_truth_nohits_ep_per;
  // TH1F* h1_truth_nohits_epx_per;
  // TH1F* h1_truth_nohits_epy_per;
  // TH1F* h1_truth_nohits_epz_per;
  // TH1F* h1_truth_nohits_etheta_per;
  // TH1F* h1_truth_nohits_ephi_per;

  // TH1F* h1_truth_nohits_piE;
  // TH1F* h1_truth_nohits_pipx;
  // TH1F* h1_truth_nohits_pipy;
  // TH1F* h1_truth_nohits_pipz;
  // TH1F* h1_truth_nohits_pitheta;
  // TH1F* h1_truth_nohits_piphi;

  // TH1F* h1_truth_nohits_piE_per;
  // TH1F* h1_truth_nohits_pip_per;
  // TH1F* h1_truth_nohits_pipx_per;
  // TH1F* h1_truth_nohits_pipy_per;
  // TH1F* h1_truth_nohits_pipz_per;
  // TH1F* h1_truth_nohits_pitheta_per;
  // TH1F* h1_truth_nohits_piphi_per;

// more comparisons
  TH2F* h2_ZDC_ThetavMom;
};

#endif // ECCE_DEMP_ANA_H

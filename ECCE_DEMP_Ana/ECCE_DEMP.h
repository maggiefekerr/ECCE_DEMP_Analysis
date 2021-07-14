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
  TVector3 eVectSmeared;
  TVector3 piVect;
  TVector3 piVectSmeared;
  TVector3 nZDCPos;
  TVector3 nZDCPosSmeared;
  TLorentzVector e4Vect;
  TLorentzVector e4VectSmeared;
  TLorentzVector pi4Vect;
  TLorentzVector pi4VectSmeared;
  TLorentzVector n4Vect;
  TLorentzVector n4VectSmeared;
  TLorentzVector e4VectTruth;
  TLorentzVector pi4VectTruth;
  TLorentzVector n4VectTruth;
  Double_t nEDep;
  Double_t nEDepSmeared;
  Double_t nTheta;
  Double_t nThetaSmeared;
  Double_t nPhi;
  Double_t nPhiSmeared;
  Double_t nPMag;
  Double_t nPMagSmeared;

  Int_t ZDC_hit;
  Int_t EEMC_hit;

  // Histogram for coincidence analysis routine

  // Resolution test plots for unsmeared vectors
  TH1F* h1_piTruth_p;
  TH1F* h1_piTruth_px;
  TH1F* h1_piTruth_py;
  TH1F* h1_piTruth_pz;
  TH1F* h1_piTruth_E;
  TH1F* h1_eTruth_p;
  TH1F* h1_eTruth_px;
  TH1F* h1_eTruth_py;
  TH1F* h1_eTruth_pz;
  TH1F* h1_eTruth_E;
  TH1F* h1_nTruth_p;
  TH1F* h1_nTruth_px;
  TH1F* h1_nTruth_py;
  TH1F* h1_nTruth_pz;
  TH1F* h1_nTruth_E;
  // Resolution test plots with smeared vectors
  TH1F* h1_piTruth_p_Smeared;
  TH1F* h1_piTruth_px_Smeared;
  TH1F* h1_piTruth_py_Smeared;
  TH1F* h1_piTruth_pz_Smeared;
  TH1F* h1_piTruth_E_Smeared;
  TH1F* h1_eTruth_p_Smeared;
  TH1F* h1_eTruth_px_Smeared;
  TH1F* h1_eTruth_py_Smeared;
  TH1F* h1_eTruth_pz_Smeared;
  TH1F* h1_eTruth_E_Smeared;
  TH1F* h1_nTruth_p_Smeared;
  TH1F* h1_nTruth_px_Smeared;
  TH1F* h1_nTruth_py_Smeared;
  TH1F* h1_nTruth_pz_Smeared;
  TH1F* h1_nTruth_E_Smeared;

  // 2D distributions 
  TH2F* h2_ZDC_XY;
  TH2F* h2_ZDC_XY_Smeared;
  // Particle Theta/Phi and Theta/p distributions
  TH2F* h2_eTrack_ThetaPhi;
  TH2F* h2_eTrack_pTheta;
  TH2F* h2_piTrack_ThetaPhi;
  TH2F* h2_piTrack_pTheta;
  TH2F* h2_nTrack_ThetaPhi;
  TH2F* h2_nTrack_pTheta;
  TH2F* h2_eTrack_ThetaPhi_Smeared;
  TH2F* h2_eTrack_pTheta_Smeared;
  TH2F* h2_piTrack_ThetaPhi_Smeared;
  TH2F* h2_piTrack_pTheta_Smeared;
  TH2F* h2_nTrack_ThetaPhi_Smeared;
  TH2F* h2_nTrack_pTheta_Smeared;
  // 2D resolution test plots
  TH2F* h2_eTruth_pxpy;
  TH2F* h2_piTruth_pxpy;
  TH2F* h2_nTruth_pxpy;
  TH2F* h2_eTruth_pxpy_Smeared;
  TH2F* h2_piTruth_pxpy_Smeared;
  TH2F* h2_nTruth_pxpy_Smeared;
 
};

#endif // ECCE_DEMP_ANA_H

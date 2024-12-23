#include "Framework/EventProcessor.h"
#include "Ecal/Event/EcalVetoResult.h"
#include "Recon/Event/TriggerResult.h"
#include "Hcal/Event/HcalVetoResult.h"
#include "DetDescr/HcalID.h"
#include "SimCore/Event/SimParticle.h"
#include "DetDescr/SimSpecialID.h"
#include "SimCore/Event/SimTrackerHit.h"
#include "Tracking/include/Tracking/Event/Track.h"
#include <math.h>

// v0: Just a few cuts, establish minimal scenario
// v1: More strict cuts
// v2: Fix for N-1 plots, fix NReadoutHits binning
// v3: Max cell dep is not depth but deposition, non-fid option
// v4: Adding HCAL PE variables, change to truth level recoil ele momentum info
// v5: Adding Target SP recoil ele variables, trigger eff, Target SP angles
// v6: Add alternative cutflow where fiducial is before trigger, add BDT score (Gabrielle), BDT vs PE
// v7: Update to latest LDMX-SW, redo trigger consitently
// v8: Moving to v1.4.1
// v9: Add the BDT cutflow too
// v10: Add all events bin to BDT cutflow, change MIP < 3
// v11: Fiducial with getFiducial()
// v12: Have BDT cutflow start with trigger
// v13: Add new cutflow with lin-reg, together with HCAL as start
// v14: Have 34 layers, add tracking plots

class CutBasedDM : public framework::Analyzer {
public:
  CutBasedDM(const std::string& name, framework::Process& p)
  : framework::Analyzer(name, p) {}
  ~CutBasedDM() = default;
  void configure(framework::config::Parameters &ps);
  void onProcessStart();
  void analyze(const framework::Event& event) final;
  template <typename T, size_t n>
  bool passPreselection(T (&passedCutsArray)[n], bool verbose);
  std::tuple<int, const ldmx::SimParticle *> getRecoilEle(
    const std::map<int, ldmx::SimParticle> &particleMap);
  std::string trigger_collName_;
  std::string trigger_passName_;
  std::string tagger_track_collection_;
  std::string recoil_track_collection_;
  bool fiducial_;
  bool acts_tracking_;
};


void CutBasedDM::configure(framework::config::Parameters &ps) {
  trigger_collName_ = ps.getParameter<std::string>("trigger_name");
  trigger_passName_ = ps.getParameter<std::string>("trigger_pass");
  fiducial_ = ps.getParameter<bool>("fiducial");
  acts_tracking_ = ps.getParameter<bool>("acts_tracking", false);
  tagger_track_collection_ = ps.getParameter<std::string>("tagger_track_collection", "TaggerTracks");
  recoil_track_collection_ = ps.getParameter<std::string>("recoil_track_collection", "RecoilTracks");
  
  return;
}

void CutBasedDM::onProcessStart(){
  getHistoDirectory();

  histograms_.create("TrigEffVsMissingE", "Triggered?", 2, -0.5, 1.5, "Missing ECAL energy [MeV]", 100, 0.0, 10000.0);
  histograms_.create("TrigEffVsRecoilPTAtTarget", "Triggered?", 2, -0.5, 1.5, "Recoil p_{T} @Target [MeV]", 800, 0.0, 10000.0);
  histograms_.create("RecoilX", "", 20, -0.5, 19.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
  histograms_.create("RecoilXAtTarget", "", 20, -0.5, 19.5, "RecoilX @Target [mm]", 90, -450.0, 450.0);
  histograms_.create("AvgLayerHit", "", 20, -0.5, 19.5, "Avg hit layer", 34, -0.5, 33.5);
  histograms_.create("DeepestLayerHit", "", 20, -0.5, 19.5, "Deepest hit layer", 34, -0.5, 33.5);
  histograms_.create("EcalBackEnergy", "", 20, -0.5, 19.5, "Ecal back energy [MeV]", 100, 0.0, 3000.0);
  histograms_.create("EpAng", "", 20, -0.5, 19.5, "EpAng", 100, 0.0, 90.0);
  histograms_.create("EpSep", "", 20, -0.5, 19.5, "EpSep", 100, 0.0, 1000.0);
  histograms_.create("FirstNearPhLayer", "", 20, -0.5, 19.5, "First near PhLayer", 34, -0.5, 33.5);
  histograms_.create("MaxCellDep", "", 20, -0.5, 19.5, "Max cell deposition [MeV]", 100, 0.0, 800.0);
  histograms_.create("NReadoutHits", "", 20, -0.5, 19.5, "#Readout hits", 150, -0.5, 149.5);
  histograms_.create("StdLayerHit", "", 20, -0.5, 19.5, "Std of hit layers", 70, -0.5, 34.5);
  histograms_.create("Straight", "", 20, -0.5, 19.5, "Straight tracks", 15, -0.5, 14.5);
  histograms_.create("LinRegNew", "", 20, -0.5, 19.5, "Linear regression tracks", 15, -0.5, 14.5);
  histograms_.create("SummedDet", "", 20, -0.5, 19.5, "Summed ECAL energy [MeV]", 100, 0.0, 10000.0);
  histograms_.create("SummedTightIso", "", 20, -0.5, 19.5, "Summed ECAL energy with tight iso [MeV]", 100, 0.0, 10000.0);
  histograms_.create("ShowerRMS", "", 20, -0.5, 19.5, "Shower RMS [mm]", 100, 0.0, 250.0);
  histograms_.create("XStd", "", 20, -0.5, 19.5, "Shower RMS_{X} [mm]", 100, 0.0, 250.0);
  histograms_.create("YStd", "", 20, -0.5, 19.5, "Shower RMS_{Y} [mm]", 100, 0.0, 250.0);
  histograms_.create("BDTDiscr", "", 20, -0.5, 19.5, "BDT discriminating score", 100, 0.0, 1.0);
  histograms_.create("BDTDiscrLog", "", 20, -0.5, 19.5, "-log(1-BDT discriminating score)", 100, 0.0, 5.0);
  
  histograms_.create("RecoilPT", "", 20, -0.5, 19.5, "Recoil p_{T} [MeV]", 800, 0.0, 10000.0);
  histograms_.create("RecoilPZ", "", 20, -0.5, 19.5, "Recoil p_{Z} [MeV]", 800, -10.0, 8010.0);
  histograms_.create("RecoilP", "", 20, -0.5, 19.5, "Recoil p [MeV]", 800, 0.0, 10000.0);
  histograms_.create("RecoilPTAtTarget", "", 20, -0.5, 19.5, "Recoil p_{T} @Target [MeV]", 800, 0.0, 10000.0);
  histograms_.create("RecoilPZAtTarget", "", 20, -0.5, 19.5, "Recoil p_{Z} @Target [MeV]", 800, -10.0, 8010.0);
  histograms_.create("RecoilPAtTarget", "", 20, -0.5, 19.5, "Recoil p @Target [MeV]", 800, 0.0, 10000.0);
  histograms_.create("RecoilTheta", "", 20, -0.5, 19.5, "Recoil theta @Target", 90, 0.0, 90.0);
  histograms_.create("RecoilPhi", "", 20, -0.5, 19.5, "Recoil phi @Target", 360, -180.0, 180.0);


  histograms_.create("Hcal_MaxPE", "", 20, -0.5, 19.5, "HCAL max photo-electron hits", 65, -0.5, 64.5);
  histograms_.create("Hcal_TotalPE", "", 20, -0.5, 19.5, "HCAL total photo-electron hits", 100, -0.5, 200.5);
  histograms_.create("Hcal_MaxTiming", "", 20, -0.5, 19.5, "HCAL timing of the max PE hit", 35, 0.0, 35.0);
  histograms_.create("Hcal_MaxSector", "", 20, -0.5, 19.5, "", 5, -0.5, 4.5);

  histograms_.create("BDTDiscrVsHcalPE_PreS", "HCAL PE", 500, -0.5, 500.5, "BDT discriminating score", 100, 0.0, 1.0);
  histograms_.create("BDTDiscrLogVsHcalPE_PreS", "HCAL PE", 500, -0.5, 500.5, "-log(1-BDT discriminating score)", 100, 0.0, 5.0);
  histograms_.create("BDTDiscrVsHcalPE_PostS", "HCAL PE", 500, -0.5, 500.5, "BDT discriminating score", 100, 0.0, 1.0);
  histograms_.create("BDTDiscrLogVsHcalPE_PostS", "HCAL PE", 500, -0.5, 500.5, "-log(1-BDT discriminating score)", 100, 0.0, 5.0);

  histograms_.create("AltCutFlow_RecoilX", "", 20, -0.5, 19.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
  histograms_.create("StdCutFlow_RecoilX", "", 20, -0.5, 19.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
  histograms_.create("BDTCutFlow_RecoilX", "", 20, -0.5, 19.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
  histograms_.create("LinRegCutFlow_RecoilX", "", 20, -0.5, 19.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
  histograms_.create("LinRegCutFlowHcal_RecoilX", "", 20, -0.5, 19.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
 
  histograms_.create("TrackingCutFlow_RecoilX", "", 20, -0.5, 19.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
  histograms_.create("Tracking_TaggerP", "", 20, -0.5, 19.5, "Tagger p [MeV]", 800, 0.0, 10000.0);
  histograms_.create("Tracking_RecoilN", "", 20, -0.5, 19.5, "N_{recoil}", 10, -0.5, 9.5);
  histograms_.create("Tracking_RecoilD0", "", 20, -0.5, 19.5, "d_{0} [mm]", 100, -50.0, 50.0);
  histograms_.create("Tracking_RecoilZ0", "", 20, -0.5, 19.5, "z_{0} [mm]", 100, -50.0, 50.0);
  histograms_.create("TrackingCutFlowHcal_RecoilX", "", 20, -0.5, 19.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
  histograms_.create("TrackingHcal_TaggerP", "", 20, -0.5, 19.5, "Tagger p [MeV]", 800, 0.0, 10000.0);
  histograms_.create("TrackingHcal_RecoilN", "", 20, -0.5, 19.5, "N_{recoil}", 10, -0.5, 9.5);
  histograms_.create("TrackingHcal_RecoilD0", "", 20, -0.5, 19.5, "d_{0} [mm]", 100, -50.0, 50.0);
  histograms_.create("TrackingHcal_RecoilZ0", "", 20, -0.5, 19.5, "z_{0} [mm]", 100, -50.0, 50.0);

// Reverse direction
  histograms_.create("Rev_RecoilX", "", 20, -0.5, 19.5, "RecoilX [mm]", 90, -450.0, 450.0);
  histograms_.create("Rev_AvgLayerHit", "", 20, -0.5, 19.5, "Avg hit layer", 34, -0.5, 33.5);
  histograms_.create("Rev_DeepestLayerHit", "", 20, -0.5, 19.5, "Deepest hit layer", 34, -0.5, 33.5);
  histograms_.create("Rev_EcalBackEnergy", "", 20, -0.5, 19.5, "Ecal back energy [MeV]", 100, 0.0, 3000.0);
  histograms_.create("Rev_EpAng", "", 20, -0.5, 19.5, "EpAng", 100, 0.0, 90.0);
  histograms_.create("Rev_EpSep", "", 20, -0.5, 19.5, "EpSep", 100, 0.0, 1000.0);
  histograms_.create("Rev_FirstNearPhLayer", "", 20, -0.5, 19.5, "First near PhLayer", 34, -0.5, 33.5);
  histograms_.create("Rev_MaxCellDep", "", 20, -0.5, 19.5, "Max cell deposition [MeV]", 100, 0.0, 800.0);
  histograms_.create("Rev_NReadoutHits", "", 20, -0.5, 19.5, "#Readout hits", 150, -0.5, 149.5);
  histograms_.create("Rev_StdLayerHit", "", 20, -0.5, 19.5, "Std of hit layers", 70, -0.5, 34.5);
  histograms_.create("Rev_Straight", "", 20, -0.5, 19.5, "Straight tracks", 15, -0.5, 14.5);
  histograms_.create("Rev_LinRegNew", "", 20, -0.5, 19.5, "Linear regression tracks", 15, -0.5, 14.5);
  histograms_.create("Rev_SummedDet", "", 20, -0.5, 19.5, "Summed ECAL energy [MeV]", 100, 0.0, 10000.0);
  histograms_.create("Rev_SummedTightIso", "", 20, -0.5, 19.5, "Summed ECAL energy with tight iso [MeV]", 100, 0.0, 10000.0);
  histograms_.create("Rev_ShowerRMS", "", 20, -0.5, 19.5, "Shower RMS [mm]", 100, 0.0, 250.0);
  histograms_.create("Rev_XStd", "", 20, -0.5, 19.5, "Shower RMS_{X} [mm]", 100, 0.0, 250.0);
  histograms_.create("Rev_YStd", "", 20, -0.5, 19.5, "Shower RMS_{Y} [mm]", 100, 0.0, 250.0);
  histograms_.create("Rev_Hcal_MaxPE", "", 20, -0.5, 19.5, "HCAL max photo-electron hits", 65, -0.5, 64.5);
  histograms_.create("Rev_Hcal_TotalPE", "", 20, -0.5, 19.5, "HCAL total photo-electron hits", 100, -0.5, 200.5);
  histograms_.create("Rev_Hcal_MaxTiming", "", 20, -0.5, 19.5, "HCAL timing of the max PE hit", 35, 0.0, 35.0);
  histograms_.create("Rev_Hcal_MaxSector", "", 20, -0.5, 19.5, "", 5, -0.5, 4.5);

  // N-1 plots
  histograms_.create("N1_EcalBackEnergy", "", 20, -0.5, 19.5, "Ecal back energy [MeV]", 100, 0.0, 3000.0);
  histograms_.create("N1_MaxCellDep", "", 20, -0.5, 19.5, "Max cell deposition [MeV]", 100, 0.0, 800.0);
  histograms_.create("N1_NReadoutHits", "", 20, -0.5, 19.5, "#Readout hits", 150, -0.5, 149.5);
  histograms_.create("N1_StdLayerHit", "", 20, -0.5, 19.5, "Std of hit layers", 70, -0.5, 34.5);
  histograms_.create("N1_Straight", "", 20, -0.5, 19.5, "Straight tracks", 15, -0.5, 14.5);
  histograms_.create("N1_SummedDet", "", 20, -0.5, 19.5, "Summed ECAL energy [MeV]", 100, 0.0, 10000.0);
  histograms_.create("N1_SummedTightIso", "", 20, -0.5, 19.5, "Summed ECAL energy with tight iso [MeV]", 100, 0.0, 10000.0);
  histograms_.create("N1_ShowerRMS", "", 20, -0.5, 19.5, "Shower RMS [mm]", 100, 0.0, 250.0);
  histograms_.create("N1_YStd", "", 20, -0.5, 19.5, "Shower RMS_{Y} [mm]", 100, 0.0, 250.0);
  histograms_.create("N1_Hcal_MaxPE", "", 20, -0.5, 19.5, "HCAL max photo-electron hits", 65, -0.5, 64.5);
  histograms_.create("N1_Hcal_TotalPE", "", 20, -0.5, 19.5, "HCAL total photo-electron hits", 100, -0.5, 200.5);
  histograms_.create("N1_Hcal_MaxTiming", "", 20, -0.5, 19.5, "HCAL timing of the max PE hit", 35, 0.0, 35.0);
  histograms_.create("N1_Hcal_MaxSector", "", 20, -0.5, 19.5, "", 5, -0.5, 4.5);

  // histograms_.create("N1_EcalBackEnergy", "", 20, -0.5, 19.5, "Ecal back energy [MeV]", 100, 0.0, 3000.0);

  std::vector<std::string> labels = {
    "All",
    "Triggerred",            // 1
    "Fiducial",              // 2
    "E_{sum} < 3500",       // 3
    "E_{SumTight} < 800",        // 4
    "E_{back} < 250",             // 5
    "N_{hits} < 70",              // 6
    "RMS_{shower} < 110",         // 7
    "RMS_{shower,Y} < 70",        // 8
    "E_{cell,max} < 300",     // 9
    "RMS_{Layer,hit} < 5",        // 10
    "N_{straight} < 3",           // 11
    "PASS_{HcalVeto}",            // 12
    ""};
    
  if (!fiducial_) labels.at(2) = "Non-fiducial";

  std::vector<TH1 *> hists = {
  histograms_.get("AvgLayerHit"),
  histograms_.get("DeepestLayerHit"),
  histograms_.get("EcalBackEnergy"),
  histograms_.get("EpAng"),
  histograms_.get("EpSep"),
  histograms_.get("FirstNearPhLayer"),
  histograms_.get("MaxCellDep"),
  histograms_.get("NReadoutHits"),
  histograms_.get("StdLayerHit"),
  histograms_.get("Straight"),
  histograms_.get("LinRegNew"),
  histograms_.get("SummedDet"),
  histograms_.get("SummedTightIso"),
  histograms_.get("ShowerRMS"),
  histograms_.get("XStd"),
  histograms_.get("YStd"),
  histograms_.get("BDTDiscr"),
  histograms_.get("BDTDiscrLog"),
  histograms_.get("StdCutFlow_RecoilX"),

  histograms_.get("RecoilX"),
  histograms_.get("RecoilPT"),
  histograms_.get("RecoilPZ"),
  histograms_.get("RecoilP"),
  histograms_.get("RecoilXAtTarget"),
  histograms_.get("RecoilPTAtTarget"),
  histograms_.get("RecoilPZAtTarget"),
  histograms_.get("RecoilPAtTarget"),
  histograms_.get("RecoilTheta"),
  histograms_.get("RecoilPhi"),
  histograms_.get("Hcal_MaxPE"),
  histograms_.get("Hcal_TotalPE"),
  histograms_.get("Hcal_MaxTiming"),
  histograms_.get("Hcal_MaxSector"),

  histograms_.get("N1_EcalBackEnergy"),
  histograms_.get("N1_MaxCellDep"),
  histograms_.get("N1_NReadoutHits"),
  histograms_.get("N1_StdLayerHit"),
  histograms_.get("N1_Straight"),
  histograms_.get("N1_SummedDet"),
  histograms_.get("N1_SummedTightIso"),
  histograms_.get("N1_ShowerRMS"),
  histograms_.get("N1_YStd"),
  histograms_.get("N1_Hcal_MaxPE"),
  histograms_.get("N1_Hcal_MaxTiming"),
  histograms_.get("N1_Hcal_MaxSector")
  };

  for (int ilabel{1}; ilabel < labels.size(); ++ilabel) {
    for (auto &hist : hists) {
      hist->GetXaxis()->SetBinLabel(ilabel, labels[ilabel - 1].c_str());
    }
  }

  std::vector<std::string> labels_Rev = {
    "All",
    "Triggerred",            // 14
    "PASS_{HcalVeto}",           // 13
    "N_{straight} < 3",           // 11
    "RMS_{Layer,hit} < 5",        // 10
    "E_{cell,max} < 300",         // 9
    "RMS_{shower,Y} < 70",        // 8
    "RMS_{shower} < 110",         // 7
    "N_{hits} < 70",              // 6
    "E_{back} < 250",             // 5
    "E_{SumTight} < 800",        // 4
    "E_{sum} < 3500",       // 3
    "Fiducial",              // 2
    ""};

  if (!fiducial_) labels_Rev.at(13) = "Non-fiducial";

  std::vector<TH1 *> hists_Rev = {
    histograms_.get("Rev_AvgLayerHit"),
    histograms_.get("Rev_DeepestLayerHit"),
    histograms_.get("Rev_EcalBackEnergy"),
    histograms_.get("Rev_EpAng"),
    histograms_.get("Rev_EpSep"),
    histograms_.get("Rev_FirstNearPhLayer"),
    histograms_.get("Rev_MaxCellDep"),
    histograms_.get("Rev_NReadoutHits"),
    histograms_.get("Rev_StdLayerHit"),
    histograms_.get("Rev_Straight"),
    histograms_.get("Rev_LinRegNew"),
    histograms_.get("Rev_SummedDet"),
    histograms_.get("Rev_SummedTightIso"),
    histograms_.get("Rev_ShowerRMS"),
    histograms_.get("Rev_XStd"),
    histograms_.get("Rev_YStd"),
    histograms_.get("Rev_Hcal_MaxPE"),
    histograms_.get("Rev_Hcal_MaxTiming"),
    histograms_.get("Rev_Hcal_MaxSector")
  };

  for (int ilabel{1}; ilabel < labels_Rev.size(); ++ilabel) {
    for (auto &hist : hists_Rev) {
      hist->GetXaxis()->SetBinLabel(ilabel, labels_Rev[ilabel - 1].c_str());
    }
  }

  std::vector<TH1 *> hists_HCALsector = {
    histograms_.get("Hcal_MaxSector"),
    histograms_.get("N1_Hcal_MaxSector"),
    histograms_.get("Rev_Hcal_MaxSector")
  };

  // enum HcalSection { BACK = 0, TOP = 1, BOTTOM = 2, RIGHT = 3, LEFT = 4 };
    std::vector<std::string> labels_HCALsector = {
    "HCAL BACK",         // 1
    "HCAL TOP",          // 2
    "HCAL BOTTOM",       // 3
    "HCAL RIGHT",        // 4
    "HCAL LEFT",         // 5
    ""};

  for (int ilabel{1}; ilabel < labels_HCALsector.size(); ++ilabel) {
    for (auto &hist : hists_HCALsector) {
      hist->GetYaxis()->SetBinLabel(ilabel, labels_HCALsector[ilabel - 1].c_str());
    }
  }

  std::vector<TH1 *> hists_AltCutFlow = {
    histograms_.get("AltCutFlow_RecoilX"),
  };

  // enum HcalSection { BACK = 0, TOP = 1, BOTTOM = 2, RIGHT = 3, LEFT = 4 };
    std::vector<std::string> labels_AltCutFlow = {
    "All",
    "Fiducial",              // 1
    "Triggerred",            // 2
    "E_{sum} < 3500",       // 3
    "E_{SumTight} < 800",        // 4
    "E_{back} < 250",             // 5
    "N_{hits} < 70",              // 6
    "RMS_{shower} < 110",         // 7
    "RMS_{shower,Y} < 70",        // 8
    "E_{cell,max} < 300",     // 9
    "RMS_{Layer,hit} < 5",        // 10
    "N_{straight} < 3",           // 11
    "PASS_{HcalVeto}",            // 12
    ""};

  if (!fiducial_) labels_AltCutFlow.at(1) = "Non-fiducial";

  for (int ilabel{1}; ilabel < labels_AltCutFlow.size(); ++ilabel) {
    for (auto &hist : hists_AltCutFlow) {
      hist->GetXaxis()->SetBinLabel(ilabel, labels_AltCutFlow[ilabel - 1].c_str());
    }
  }

  // CutFlow labels for BDT
  std::vector<TH1 *> hists_BDTCutFlow = {
    histograms_.get("BDTCutFlow_RecoilX"),
  };

    std::vector<std::string> labels_BDTCutFlow = {
    "All",
    "Triggerred",            // 1
    "Fiducial",              // 2
    "BDT passed",            // 3
    "N_{straight} < 3",      // 4
    "PASS_{HcalVeto}",       // 5
    "N_{straight} = 0",      // 6
    "Angle_{e,ph} > 3.",     // 7
    ""};

  if (!fiducial_) labels_BDTCutFlow.at(2) = "Non-fiducial";

  for (int ilabel{1}; ilabel < labels_BDTCutFlow.size(); ++ilabel) {
    for (auto &hist : hists_BDTCutFlow) {
      hist->GetXaxis()->SetBinLabel(ilabel, labels_BDTCutFlow[ilabel - 1].c_str());
    }
  }

  // CutFlow labels for BDT with new lin-reg
   std::vector<TH1 *> hists_LinRegCutFlow = {
    histograms_.get("LinRegCutFlow_RecoilX"),
  };

    std::vector<std::string> labels_LinRegCutFlow = {
    "All",
    "Triggerred",            // 1
    "Fiducial",              // 2
    "BDT passed",            // 3
    "N_{straight} < 3",      // 4
    "PASS_{HcalVeto}",       // 5
    "N_{straight} = 0",      // 6
    "N_{lin-reg} = 0",      // 7
    "Angle_{e,ph} > 3.",     // 8
    ""};

  if (!fiducial_) labels_LinRegCutFlow.at(2) = "Non-fiducial";

  for (int ilabel{1}; ilabel < labels_LinRegCutFlow.size(); ++ilabel) {
    for (auto &hist : hists_LinRegCutFlow) {
      hist->GetXaxis()->SetBinLabel(ilabel, labels_LinRegCutFlow[ilabel - 1].c_str());
    }
  }

  // CutFlow labels for BDT with new lin-reg
   std::vector<TH1 *> hists_LinRegCutFlowHcal = {
    histograms_.get("LinRegCutFlowHcal_RecoilX"),
  };

    std::vector<std::string> labels_LinRegCutFlowHcal = {
    "All",
    "Triggerred",            // 1
    "PASS_{HcalVeto}",       // 2
    "Fiducial",              // 3
    "BDT passed",            // 4
    "N_{straight} = 0",      // 5
    "N_{lin-reg} = 0",      // 6
    "Angle_{e,ph} > 3.",     // 7
    ""};

  if (!fiducial_) labels_LinRegCutFlowHcal.at(3) = "Non-fiducial";

  for (int ilabel{1}; ilabel < labels_LinRegCutFlowHcal.size(); ++ilabel) {
    for (auto &hist : hists_LinRegCutFlowHcal) {
      hist->GetXaxis()->SetBinLabel(ilabel, labels_LinRegCutFlowHcal[ilabel - 1].c_str());
    }
  }

  // CutFlow labels for BDT with tracking
   std::vector<TH1 *> hists_TrackingCutFlow = {
    histograms_.get("TrackingCutFlow_RecoilX"),
    histograms_.get("Tracking_TaggerP"),
    histograms_.get("Tracking_RecoilN"),
    histograms_.get("Tracking_RecoilD0"),
    histograms_.get("Tracking_RecoilZ0"),
  };

    std::vector<std::string> labels_TrackingCutFlow = {
    "All",
    "Triggerred",            // 1
    "p_{tagger} < 5600",     // 2
    "N_{recoil} = 1",  		 // 3
    "|d_{0}| < 10",  		 // 4
    "|z_{0}| < 40",  		 // 5
    "Fiducial",              // 6
    "BDT passed",            // 7
    "N_{straight} < 3",      // 8
    "PASS_{HcalVeto}",       // 9
    "N_{straight} = 0",      // 10
    "Angle_{e,ph} > 3.",     // 11
    ""};

  if (!fiducial_) labels_TrackingCutFlow.at(6) = "Non-fiducial";

  for (int ilabel{1}; ilabel < labels_TrackingCutFlow.size(); ++ilabel) {
    for (auto &hist : hists_TrackingCutFlow) {
      hist->GetXaxis()->SetBinLabel(ilabel, labels_TrackingCutFlow[ilabel - 1].c_str());
    }
  }

// CutFlow labels for BDT with tracking with Hcal first
   std::vector<TH1 *> hists_TrackingCutFlowHcal = {
    histograms_.get("TrackingCutFlowHcal_RecoilX"),
    histograms_.get("TrackingHcal_TaggerP"),
    histograms_.get("TrackingHcal_RecoilN"),
    histograms_.get("TrackingHcal_RecoilD0"),
    histograms_.get("TrackingHcal_RecoilZ0"),
  };

    std::vector<std::string> labels_TrackingCutFlowHcal = {
    "All",
    "Triggerred",            // 1
    "PASS_{HcalVeto}",       // 2
    "Fiducial",              // 3
    "BDT passed",            // 4
    "N_{straight} = 0",      // 5
    "Angle_{e,ph} > 3.",     // 6
    "p_{tagger} < 5600",     // 7
    "N_{recoil} = 1",  		 // 8
    "|d_{0}| < 10",  		 // 9
    "|z_{0}| < 40",  		 // 10
    ""};

  if (!fiducial_) labels_TrackingCutFlowHcal.at(4) = "Non-fiducial";

  for (int ilabel{1}; ilabel < labels_TrackingCutFlowHcal.size(); ++ilabel) {
    for (auto &hist : hists_TrackingCutFlowHcal) {
      hist->GetXaxis()->SetBinLabel(ilabel, labels_TrackingCutFlowHcal[ilabel - 1].c_str());
    }
  }

}

void CutBasedDM::analyze(const framework::Event& event) {
  // std::cout << " ---------------------------------------------" << std::endl;
  auto vetoNew{event.getObject<ldmx::EcalVetoResult>("EcalVetoNew","")};
  auto trigResult{event.getObject<ldmx::TriggerResult>(trigger_collName_, trigger_passName_)};
  auto hcalVeto{event.getObject<ldmx::HcalVetoResult>("HcalVeto","cutbased")};
  auto hcalRecHits{event.getCollection<ldmx::HcalHit>("HcalRecHits", "")};
  auto targetSpHits{event.getCollection<ldmx::SimTrackerHit>("TargetScoringPlaneHits")};
  std::shared_ptr<ldmx::Tracks> taggerTrackCollection_{nullptr};
  std::shared_ptr<ldmx::Tracks> recoilTrackCollection_{nullptr};
  if (event.exists(tagger_track_collection_)) {
    taggerTrackCollection_ = std::make_shared<ldmx::Tracks>(
        event.getCollection<ldmx::Track>(tagger_track_collection_));
  }
  if (event.exists(recoil_track_collection_)) {
    recoilTrackCollection_ = std::make_shared<ldmx::Tracks>(
        event.getCollection<ldmx::Track>(recoil_track_collection_));
  }

  // Take recoil momentum from SIM particles
  float pT{-9999.};
  float pZ{-9999.};
  float totMom{-9999.};
  auto particleMap{event.getMap<int, ldmx::SimParticle>("SimParticles")};
  auto [recoilTrackID, recoilElectron] = CutBasedDM::getRecoilEle(particleMap);
  pT =   sqrt(recoilElectron->getMomentum()[0] * recoilElectron->getMomentum()[0] +  recoilElectron->getMomentum()[1] * recoilElectron->getMomentum()[1]);
  pZ = recoilElectron->getMomentum()[2];
  totMom = sqrt(pT*pT + pZ*pZ);

  // auto phiEle =   (180/M_PI)*std::acos(recoilElectron->getMomentum()[1] /totMom)-90.;
  // auto thetaEle = (180/M_PI)*std::acos(pZ/totMom);

  //  Same but at the target SP
  float XAtTarget{-9999};
  float pTAtTarget{-9999.};
  float pYAtTarget{-9999.};
  float pZAtTarget{-9999.};
  float totMomAtTarget{-9999.};
  for (ldmx::SimTrackerHit &spHit : targetSpHits) {
    ldmx::SimSpecialID hit_id(spHit.getID());
    if (hit_id.plane() != 1 || spHit.getMomentum()[2] <= 0) continue;

    if (spHit.getTrackID() == recoilTrackID) {
      float p_current = sqrt(spHit.getMomentum()[0]*spHit.getMomentum()[0] + spHit.getMomentum()[1]*spHit.getMomentum()[1] + spHit.getMomentum()[2]*spHit.getMomentum()[2]);
      if (p_current > totMomAtTarget) {
        pYAtTarget  = spHit.getMomentum()[1];
        pTAtTarget = sqrt(spHit.getMomentum()[0]*spHit.getMomentum()[0] + spHit.getMomentum()[1]*spHit.getMomentum()[1]);
        pZAtTarget = spHit.getMomentum()[2];
        totMomAtTarget = p_current;
        XAtTarget = spHit.getPosition()[0];
      }
    }
  }
//  std::cout << " pYAtTarget " <<  pYAtTarget << " totMomAtTarget " << totMomAtTarget ;
  auto phiEleAtTarget =   (180/M_PI)*std::acos(pYAtTarget/totMomAtTarget)-90.;
  auto thetaEleAtTarget = (180/M_PI)*std::acos(pZAtTarget/totMomAtTarget);
  // std::cout << " phiEleAtTarget " << phiEleAtTarget << " phiEle " << phiEle
  //  std::cout << " thetaEleAtTarget " << thetaEleAtTarget << " phiEleAtTarget " << phiEleAtTarget << std::endl;

  // Take recoil momentum from ECAL SP
  // pT2 = vetoNew.getRecoilMomentum()[0]*vetoNew.getRecoilMomentum()[0] + vetoNew.getRecoilMomentum()[1]*vetoNew.getRecoilMomentum()[1];
  // pZ =  vetoNew.getRecoilMomentum()[2];


  // HCAL veto calc
  float hcalMaxPE{-9999};
  float hcalMaxTiming{-9999};
  int hcalMaxSector{-1};
  float hcalTotalPe{0};
  
  ldmx::HcalHit defaultMaxHit_;
  defaultMaxHit_.Clear();
  defaultMaxHit_.setPE(-9999);
  defaultMaxHit_.setMinPE(-9999);
  defaultMaxHit_.setSection(-9999);
  defaultMaxHit_.setLayer(-9999);
  defaultMaxHit_.setStrip(-9999);
  defaultMaxHit_.setEnd(-999);
  defaultMaxHit_.setTimeDiff(-9999);
  defaultMaxHit_.setToaPos(-9999);
  defaultMaxHit_.setToaNeg(-9999);
  defaultMaxHit_.setAmplitudePos(-9999);
  defaultMaxHit_.setAmplitudeNeg(-9999);
  const ldmx::HcalHit *maxPEHit{&defaultMaxHit_};

  // Loop on the HCAL hits
  for (const auto hcalHit : hcalRecHits) {
    float maxTime_{50.};
    if (hcalHit.getTime() >= maxTime_) {
      continue;
    }

    float backMinPE_{1.};
    ldmx::HcalID id(hcalHit.getID());
    if ((id.section() == ldmx::HcalID::BACK) &&  (hcalHit.getMinPE() < backMinPE_)) {
      continue;
    }

    float pe = hcalHit.getPE();
    hcalTotalPe += pe;

    if (hcalMaxPE < pe) {
      hcalMaxPE = pe;
      maxPEHit = &hcalHit;
    }
  }

  ldmx::HcalID maxPeId(maxPEHit->getID());
  hcalMaxTiming = maxPEHit->getTime();
  hcalMaxSector =  maxPeId.section();

  // Trigger eff curves
  histograms_.fill("TrigEffVsMissingE", trigResult.passed() , 8000.-vetoNew.getSummedDet() );
  histograms_.fill("TrigEffVsRecoilPTAtTarget", trigResult.passed() , pTAtTarget );

  // CutFlow here
  bool passedCutsArrayCnC[12];
  std::fill(std::begin(passedCutsArrayCnC), std::end(passedCutsArrayCnC),false);
  passedCutsArrayCnC[0]  = (trigResult.passed()) ? true : false;
  passedCutsArrayCnC[1]  = ((fiducial_ && vetoNew.getFiducial()) || (!fiducial_ && !vetoNew.getFiducial())) ? true : false;
  passedCutsArrayCnC[2]  = (vetoNew.getSummedDet() < 3500) ? true : false;
  passedCutsArrayCnC[3]  = (vetoNew.getSummedTightIso() < 800) ? true : false;
  passedCutsArrayCnC[4]  = (vetoNew.getEcalBackEnergy() < 250) ? true : false;
  passedCutsArrayCnC[5]  = (vetoNew.getNReadoutHits() < 70) ? true : false;
  passedCutsArrayCnC[6]  = (vetoNew.getShowerRMS() < 110) ? true : false;
  passedCutsArrayCnC[7]  = (vetoNew.getYStd() < 70) ? true : false;
  passedCutsArrayCnC[8]  = (vetoNew.getMaxCellDep() < 300) ? true : false;
  passedCutsArrayCnC[9]  = (vetoNew.getStdLayerHit() < 5) ? true : false;
  passedCutsArrayCnC[10]  = (vetoNew.getNStraightTracks() < 3) ? true : false;
  passedCutsArrayCnC[11]  = (hcalVeto.passesVeto()) ? true : false;

  // Fill histograms
  histograms_.fill("RecoilX", 0. , vetoNew.getRecoilX() );
  histograms_.fill("AvgLayerHit", 0. , vetoNew.getAvgLayerHit() );
  histograms_.fill("DeepestLayerHit", 0. , vetoNew.getDeepestLayerHit() );
  histograms_.fill("EcalBackEnergy", 0. , vetoNew.getEcalBackEnergy() );
  histograms_.fill("EpAng", 0. , vetoNew.getEPAng() );
  histograms_.fill("EpSep", 0. , vetoNew.getEPSep() );
  histograms_.fill("FirstNearPhLayer", 0. , vetoNew.getFirstNearPhLayer() );
  histograms_.fill("MaxCellDep", 0. , vetoNew.getMaxCellDep() );
  histograms_.fill("NReadoutHits", 0. , vetoNew.getNReadoutHits() );
  histograms_.fill("StdLayerHit", 0. , vetoNew.getStdLayerHit() );
  histograms_.fill("Straight", 0. , vetoNew.getNStraightTracks() );
  histograms_.fill("LinRegNew", 0. , vetoNew.getNLinRegTracks() );
  histograms_.fill("SummedDet", 0. , vetoNew.getSummedDet() );
  histograms_.fill("SummedTightIso", 0. , vetoNew.getSummedTightIso() );
  histograms_.fill("ShowerRMS", 0. , vetoNew.getShowerRMS() );
  histograms_.fill("XStd", 0. , vetoNew.getXStd() );
  histograms_.fill("YStd", 0. , vetoNew.getYStd() );
  histograms_.fill("BDTDiscr", 0. , vetoNew.getDisc() );
  histograms_.fill("BDTDiscrLog", 0. , -log(1-vetoNew.getDisc()) );
  histograms_.fill("StdCutFlow_RecoilX", 0. , vetoNew.getRecoilX() );
  
  histograms_.fill("RecoilPT", 0. ,pT );
  histograms_.fill("RecoilPZ", 0. , pZ );
  histograms_.fill("RecoilP", 0. , totMom );
  histograms_.fill("RecoilXAtTarget", 0. ,XAtTarget );
  histograms_.fill("RecoilPTAtTarget", 0. ,pTAtTarget );
  histograms_.fill("RecoilPZAtTarget", 0. , pZAtTarget );
  histograms_.fill("RecoilPAtTarget", 0. , totMomAtTarget );
  histograms_.fill("RecoilTheta", 0. , thetaEleAtTarget );
  histograms_.fill("RecoilPhi", 0. , phiEleAtTarget );
  histograms_.fill("Hcal_MaxPE", 0. , hcalMaxPE );
  histograms_.fill("Hcal_TotalPE", 0. , hcalTotalPe );
  histograms_.fill("Hcal_MaxTiming", 0. , hcalMaxTiming );
  histograms_.fill("Hcal_MaxSector", 0. , hcalMaxSector );

  histograms_.fill("AltCutFlow_RecoilX", 0 , vetoNew.getRecoilX() );
  histograms_.fill("BDTCutFlow_RecoilX", 0 , vetoNew.getRecoilX() );
  histograms_.fill("LinRegCutFlow_RecoilX", 0 , vetoNew.getRecoilX() );
  histograms_.fill("LinRegCutFlowHcal_RecoilX", 0 , vetoNew.getRecoilX() );
  
  for (size_t i=0;i<sizeof(passedCutsArrayCnC);i++) {
    bool allCutsPassedSoFar = true;
    for (size_t j=0;j<=i;j++) {
      if (!passedCutsArrayCnC[j]) {
        allCutsPassedSoFar = false;
        break;
      }
    }
    if (allCutsPassedSoFar) {
      // std::cout 
      //   << " i-th cut = " << i 
      //   << " trigger = " << trigResult.passed() 
      //   << " getRecoilX = " << vetoNew.getRecoilX() 
      //   << " getSummedDet = " << vetoNew.getSummedDet()
      //   << " getSummedTightIso = " << vetoNew.getSummedTightIso() 
      //   << " getEcalBackEnergy = " << vetoNew.getEcalBackEnergy() 
      //   << " getNReadoutHits = " << vetoNew.getNReadoutHits()
      //   << " getShowerRMS = " << vetoNew.getShowerRMS()
      //   << " getYStd = " << vetoNew.getYStd()
      //   << " getMaxCellDep = " << vetoNew.getMaxCellDep()
      //   << " getStdLayerHit = " << vetoNew.getStdLayerHit()
      //   << " getNStraightTracks = " << vetoNew.getNStraightTracks()
      //   << " hcalVeto = " << hcalVeto.passesVeto()
      // std::cout << "hcalMaxPE = " <<  hcalMaxPE << " hcal total" << hcalTotalPe <<  " maxTime = " << hcalMaxTiming << " where = " << hcalMaxSector
      // << std::endl;

      histograms_.fill("RecoilX", i+1 , vetoNew.getRecoilX() );
      histograms_.fill("AvgLayerHit", i+1 , vetoNew.getAvgLayerHit() );
      histograms_.fill("DeepestLayerHit", i+1 , vetoNew.getDeepestLayerHit() );
      histograms_.fill("EcalBackEnergy", i+1 , vetoNew.getEcalBackEnergy() );
      histograms_.fill("EpAng", i+1 , vetoNew.getEPAng() );
      histograms_.fill("EpSep", i+1 , vetoNew.getEPSep() );
      histograms_.fill("FirstNearPhLayer", i+1 , vetoNew.getFirstNearPhLayer() );
      histograms_.fill("MaxCellDep", i+1 , vetoNew.getMaxCellDep() );
      histograms_.fill("NReadoutHits", i+1 , vetoNew.getNReadoutHits() );
      histograms_.fill("StdLayerHit", i+1 , vetoNew.getStdLayerHit() );
      histograms_.fill("Straight", i+1 , vetoNew.getNStraightTracks() );
      histograms_.fill("LinRegNew", i+1 , vetoNew.getNLinRegTracks() );
      histograms_.fill("SummedDet", i+1 , vetoNew.getSummedDet() );
      histograms_.fill("SummedTightIso", i+1 , vetoNew.getSummedTightIso() );
      histograms_.fill("ShowerRMS", i+1 , vetoNew.getShowerRMS() );
      histograms_.fill("XStd", i+1 , vetoNew.getXStd() );
      histograms_.fill("YStd", i+1 , vetoNew.getYStd() );
      histograms_.fill("BDTDiscr", i+1 , vetoNew.getDisc() );
      histograms_.fill("BDTDiscrLog", i+1 , -log(1-vetoNew.getDisc()) );
      histograms_.fill("StdCutFlow_RecoilX", i+1 , vetoNew.getRecoilX() );
      histograms_.fill("RecoilPT", i+1 ,pT );
      histograms_.fill("RecoilPZ", i+1 , pZ );
      histograms_.fill("RecoilP", i+1 , totMom );
      histograms_.fill("RecoilXAtTarget", i+1 ,XAtTarget );
      histograms_.fill("RecoilPTAtTarget", i+1 ,pTAtTarget );
      histograms_.fill("RecoilPZAtTarget", i+1 , pZAtTarget );
      histograms_.fill("RecoilPAtTarget", i+1 , totMomAtTarget );
      histograms_.fill("RecoilTheta", i+1. , thetaEleAtTarget );
      histograms_.fill("RecoilPhi", i+1 , phiEleAtTarget );
      histograms_.fill("Hcal_MaxPE", i+1. , hcalMaxPE );
      histograms_.fill("Hcal_TotalPE", i+1. , hcalTotalPe );
      histograms_.fill("Hcal_MaxTiming", i+1. , hcalMaxTiming );
      histograms_.fill("Hcal_MaxSector", i+1. , hcalMaxSector );

      if (i==1) {
        histograms_.fill("BDTDiscrVsHcalPE_PreS", hcalMaxPE , vetoNew.getDisc() );
        histograms_.fill("BDTDiscrLogVsHcalPE_PreS", hcalMaxPE , -log(1-vetoNew.getDisc()) );
      } 
      if (i==10) {
        histograms_.fill("BDTDiscrVsHcalPE_PostS", hcalMaxPE , vetoNew.getDisc() );
        histograms_.fill("BDTDiscrLogVsHcalPE_PostS", hcalMaxPE , -log(1-vetoNew.getDisc()) );
      }
    }
  }

  // Alternative cutFlow here
  bool passedCutsArrayAlt[12];
  std::fill(std::begin(passedCutsArrayAlt), std::end(passedCutsArrayAlt),false);
  passedCutsArrayAlt[0]  = ((fiducial_ && vetoNew.getFiducial()) || (!fiducial_ && !vetoNew.getFiducial())) ? true : false;
  passedCutsArrayAlt[1]  = (trigResult.passed()) ? true : false;
  passedCutsArrayAlt[2]  = (vetoNew.getSummedDet() < 3500) ? true : false;
  passedCutsArrayAlt[3]  = (vetoNew.getSummedTightIso() < 800) ? true : false;
  passedCutsArrayAlt[4]  = (vetoNew.getEcalBackEnergy() < 250) ? true : false;
  passedCutsArrayAlt[5]  = (vetoNew.getNReadoutHits() < 70) ? true : false;
  passedCutsArrayAlt[6]  = (vetoNew.getShowerRMS() < 110) ? true : false;
  passedCutsArrayAlt[7]  = (vetoNew.getYStd() < 70) ? true : false;
  passedCutsArrayAlt[8]  = (vetoNew.getMaxCellDep() < 300) ? true : false;
  passedCutsArrayAlt[9]  = (vetoNew.getStdLayerHit() < 5) ? true : false;
  passedCutsArrayAlt[10]  = (vetoNew.getNStraightTracks() < 3) ? true : false;
  passedCutsArrayAlt[11]  = (hcalVeto.passesVeto()) ? true : false;

  for (size_t i=0;i<sizeof(passedCutsArrayAlt);i++) {
    bool allCutsPassedSoFar = true;
    for (size_t j=0;j<=i;j++) {
      if (!passedCutsArrayAlt[j]) {
        allCutsPassedSoFar = false;
        break;
      }
    }
    if (allCutsPassedSoFar) {
      histograms_.fill("AltCutFlow_RecoilX", i+1 , vetoNew.getRecoilX() );
    }
  }

  // BDT based cutFlow here
  bool passedCutsArrayBDT[7];
  std::fill(std::begin(passedCutsArrayBDT), std::end(passedCutsArrayBDT),false);
  passedCutsArrayBDT[0]  = (trigResult.passed()) ? true : false;
  passedCutsArrayBDT[1]  = ((fiducial_ && vetoNew.getFiducial()) || (!fiducial_ && !vetoNew.getFiducial())) ? true : false;
  passedCutsArrayBDT[2]  = (vetoNew.getDisc() > 0.99741) ? true : false;
  passedCutsArrayBDT[3]  = (vetoNew.getNStraightTracks() < 3) ? true : false;
  passedCutsArrayBDT[4]  = (hcalVeto.passesVeto()) ? true : false;
  passedCutsArrayBDT[5]  = (vetoNew.getNStraightTracks() == 0) ? true : false;
  passedCutsArrayBDT[6]  = ((vetoNew.getEPAng() > 3) && (fiducial_ && vetoNew.getEPAng()  < 999) || (!fiducial_ )) ? true : false;


  for (size_t i=0;i<sizeof(passedCutsArrayBDT);i++) {
    bool allCutsPassedSoFar = true;
    for (size_t j=0;j<=i;j++) {
      if (!passedCutsArrayBDT[j]) {
        allCutsPassedSoFar = false;
        break;
      }
    }
    if (allCutsPassedSoFar) {
      histograms_.fill("BDTCutFlow_RecoilX", i+1 , vetoNew.getRecoilX() );
    }
  }

  // BDT based cutFlow here with linreg
  bool passedCutsArrayLinReg[8];
  std::fill(std::begin(passedCutsArrayLinReg), std::end(passedCutsArrayLinReg),false);
  passedCutsArrayLinReg[0]  = (trigResult.passed()) ? true : false;
  passedCutsArrayLinReg[1]  = ((fiducial_ && vetoNew.getFiducial()) || (!fiducial_ && !vetoNew.getFiducial())) ? true : false;
  passedCutsArrayLinReg[2]  = (vetoNew.getDisc() > 0.99741) ? true : false;
  passedCutsArrayLinReg[3]  = (vetoNew.getNStraightTracks() < 3) ? true : false;
  passedCutsArrayLinReg[4]  = (hcalVeto.passesVeto()) ? true : false;
  passedCutsArrayLinReg[5]  = (vetoNew.getNStraightTracks() == 0) ? true : false;
  passedCutsArrayLinReg[6]  = (vetoNew.getNLinRegTracks() == 0) ? true : false;
  passedCutsArrayLinReg[7]  = ((vetoNew.getEPAng() > 3) && (fiducial_ && vetoNew.getEPAng()  < 999) || (!fiducial_ )) ? true : false;

  for (size_t i=0;i<sizeof(passedCutsArrayLinReg);i++) {
    bool allCutsPassedSoFar = true;
    for (size_t j=0;j<=i;j++) {
      if (!passedCutsArrayLinReg[j]) {
        allCutsPassedSoFar = false;
        break;
      }
    }
    if (allCutsPassedSoFar) {
      histograms_.fill("LinRegCutFlow_RecoilX", i+1 , vetoNew.getRecoilX() );
    }
  }

    // BDT based cutFlow here with linreg, starting with Hcal
  bool passedCutsArrayLinRegHcal[6];
  std::fill(std::begin(passedCutsArrayLinRegHcal), std::end(passedCutsArrayLinRegHcal),false);
  passedCutsArrayLinRegHcal[0]  = (trigResult.passed()) ? true : false;
  passedCutsArrayLinRegHcal[1]  = (hcalVeto.passesVeto()) ? true : false;
  passedCutsArrayLinRegHcal[2]  = ((fiducial_ && vetoNew.getFiducial()) || (!fiducial_ && !vetoNew.getFiducial())) ? true : false;
  passedCutsArrayLinRegHcal[3]  = (vetoNew.getDisc() > 0.99741) ? true : false;
  passedCutsArrayLinRegHcal[4]  = (vetoNew.getNStraightTracks() == 0) ? true : false;
  passedCutsArrayLinRegHcal[5]  = (vetoNew.getNLinRegTracks() == 0) ? true : false;


  for (size_t i=0;i<sizeof(passedCutsArrayLinRegHcal);i++) {
    bool allCutsPassedSoFar = true;
    for (size_t j=0;j<=i;j++) {
      if (!passedCutsArrayLinRegHcal[j]) {
        allCutsPassedSoFar = false;
        break;
      }
    }
    if (allCutsPassedSoFar) {

      histograms_.fill("LinRegCutFlowHcal_RecoilX", i+1 , vetoNew.getRecoilX() );
    }
  }
  // --------------------------------------------------------------------------
  // Calculate tracking variables if tracking is available
  float taggerP{0.0}; // Make sure this is in MeV!!
  int recoilN{0};
  float recoilD0{-9999.};
  float recoilZ0{-9999.};
  if (acts_tracking_) {
    // Start with tagger tracks
    auto taggerN = taggerTrackCollection_->size();
    if (taggerN == 1) {
      for (auto& trk : *(taggerTrackCollection_)) {
        auto QoP = trk.getQoP();
        taggerP = 1000. / abs(QoP);
      }
    }

    // Recoil tracks now
    recoilN = recoilTrackCollection_->size();
    if (recoilN == 1) {
      for (auto& trk : *(recoilTrackCollection_)) {
        recoilD0 = trk.getD0();
        recoilZ0 = trk.getZ0();
      }
    }
  }
  // --------------------------------------------------------------------------
  // BDT based cutFlow with tracking
  histograms_.fill("TrackingCutFlow_RecoilX", 0 , vetoNew.getRecoilX() );
  histograms_.fill("Tracking_TaggerP", 0 , taggerP);
  histograms_.fill("Tracking_RecoilN", 0 , recoilN);
  histograms_.fill("Tracking_RecoilD0", 0, recoilD0);
  histograms_.fill("Tracking_RecoilZ0", 0 , recoilZ0);

  bool passedCutsArrayTracking[11];
  std::fill(std::begin(passedCutsArrayTracking), std::end(passedCutsArrayTracking),false);
  passedCutsArrayTracking[0]  = (trigResult.passed()) ? true : false;
  passedCutsArrayTracking[1]  = (taggerP > 5600) ? true : false;
  passedCutsArrayTracking[2]  = (recoilN == 1) ? true : false;
  passedCutsArrayTracking[3]  = (abs(recoilD0) < 10) ? true : false;
  passedCutsArrayTracking[4]  = (abs(recoilZ0) < 40) ? true : false;
  passedCutsArrayTracking[5]  = ((fiducial_ && vetoNew.getFiducial()) || (!fiducial_ && !vetoNew.getFiducial())) ? true : false;
  passedCutsArrayTracking[6]  = (vetoNew.getDisc() > 0.99741) ? true : false;
  passedCutsArrayTracking[7]  = (vetoNew.getNStraightTracks() < 3) ? true : false;
  passedCutsArrayTracking[8]  = (hcalVeto.passesVeto()) ? true : false;
  passedCutsArrayTracking[9]  = (vetoNew.getNStraightTracks() == 0) ? true : false;
  passedCutsArrayTracking[10]  = ((vetoNew.getEPAng() > 3) && (fiducial_ && vetoNew.getEPAng()  < 999) || (!fiducial_ )) ? true : false;


  for (size_t i=0;i<sizeof(passedCutsArrayTracking);i++) {
    bool allCutsPassedSoFar = true;
    for (size_t j=0;j<=i;j++) {
      if (!passedCutsArrayTracking[j]) {
        allCutsPassedSoFar = false;
        break;
      }
    }
    if (allCutsPassedSoFar) {
      histograms_.fill("TrackingCutFlow_RecoilX", i+1 , vetoNew.getRecoilX() );
      histograms_.fill("Tracking_TaggerP", i+1 , taggerP);
      histograms_.fill("Tracking_RecoilN", i+1 , recoilN);
      histograms_.fill("Tracking_RecoilD0", i+1 , recoilD0);
      histograms_.fill("Tracking_RecoilZ0", i+1 , recoilZ0);
    }
  }

  // BDT based cutFlow with tracking starting with Hcal and Ecal veto
  histograms_.fill("TrackingCutFlowHcal_RecoilX", 0 , vetoNew.getRecoilX() );
  histograms_.fill("TrackingHcal_TaggerP", 0 , taggerP);
  histograms_.fill("TrackingHcal_RecoilN", 0 , recoilN);
  histograms_.fill("TrackingHcal_RecoilD0", 0 , recoilD0);
  histograms_.fill("TrackingHcal_RecoilZ0", 0 , recoilZ0);

  bool passedCutsArrayTrackingHcal[10];
  std::fill(std::begin(passedCutsArrayTrackingHcal), std::end(passedCutsArrayTrackingHcal),false);
  passedCutsArrayTrackingHcal[0]  = (trigResult.passed()) ? true : false;
  passedCutsArrayTrackingHcal[1]  = (hcalVeto.passesVeto()) ? true : false;
  passedCutsArrayTrackingHcal[2]  = ((fiducial_ && vetoNew.getFiducial()) || (!fiducial_ && !vetoNew.getFiducial())) ? true : false;
  passedCutsArrayTrackingHcal[3]  = (vetoNew.getDisc() > 0.99741) ? true : false;
  passedCutsArrayTrackingHcal[4]  = (vetoNew.getNStraightTracks() == 0) ? true : false;
  passedCutsArrayTrackingHcal[5]  = ((vetoNew.getEPAng() > 3) && (fiducial_ && vetoNew.getEPAng()  < 999) || (!fiducial_ )) ? true : false;
  passedCutsArrayTrackingHcal[6]  = (taggerP > 5600) ? true : false;
  passedCutsArrayTrackingHcal[7]  = (recoilN == 1) ? true : false;
  passedCutsArrayTrackingHcal[8]  = (abs(recoilD0) < 10) ? true : false;
  passedCutsArrayTrackingHcal[9]  = (abs(recoilZ0) < 40) ? true : false;


  for (size_t i=0;i<sizeof(passedCutsArrayTrackingHcal);i++) {
    bool allCutsPassedSoFar = true;
    for (size_t j=0;j<=i;j++) {
      if (!passedCutsArrayTrackingHcal[j]) {
        allCutsPassedSoFar = false;
        break;
      }
    }
    if (allCutsPassedSoFar) {
      histograms_.fill("TrackingCutFlowHcal_RecoilX", i+1 , vetoNew.getRecoilX() );
      histograms_.fill("TrackingHcal_TaggerP", i+1 , taggerP);
      histograms_.fill("TrackingHcal_RecoilN", i+1 , recoilN);
      histograms_.fill("TrackingHcal_RecoilD0", i+1 , recoilD0);
      histograms_.fill("TrackingHcal_RecoilZ0", i+1 , recoilZ0);
    }
  }

  // --------------------------------------------------------------------------
   // Reverse cutflow, i.e. start with the last cut from the original cutflow
histograms_.fill("Rev_AvgLayerHit", 0. , vetoNew.getAvgLayerHit() );
histograms_.fill("Rev_DeepestLayerHit", 0. , vetoNew.getDeepestLayerHit() );
histograms_.fill("Rev_EcalBackEnergy", 0. , vetoNew.getEcalBackEnergy() );
histograms_.fill("Rev_EpAng", 0. , vetoNew.getEPAng() );
histograms_.fill("Rev_EpSep", 0. , vetoNew.getEPSep() );
histograms_.fill("Rev_FirstNearPhLayer", 0. , vetoNew.getFirstNearPhLayer() );
histograms_.fill("Rev_MaxCellDep", 0. , vetoNew.getMaxCellDep() );
histograms_.fill("Rev_NReadoutHits", 0. , vetoNew.getNReadoutHits() );
histograms_.fill("Rev_StdLayerHit", 0. , vetoNew.getStdLayerHit() );
histograms_.fill("Rev_Straight", 0. , vetoNew.getNStraightTracks() );
histograms_.fill("Rev_LinRegNew", 0. , vetoNew.getNLinRegTracks() );
histograms_.fill("Rev_SummedDet", 0. , vetoNew.getSummedDet() );
histograms_.fill("Rev_SummedTightIso", 0. , vetoNew.getSummedTightIso() );
histograms_.fill("Rev_ShowerRMS", 0. , vetoNew.getShowerRMS() );
histograms_.fill("Rev_XStd", 0. , vetoNew.getXStd() );
histograms_.fill("Rev_YStd", 0. , vetoNew.getYStd() );
histograms_.fill("Rev_Hcal_MaxPE", 0. , hcalMaxPE );
histograms_.fill("Rev_Hcal_TotalPE", 0. , hcalTotalPe );
histograms_.fill("Rev_Hcal_MaxTiming", 0. , hcalMaxTiming );
histograms_.fill("Rev_Hcal_MaxSector", 0. , hcalMaxSector );

bool passedCutsArrayReverse[12];
std::fill(std::begin(passedCutsArrayReverse), std::end(passedCutsArrayReverse),false);
passedCutsArrayReverse[0]  = (trigResult.passed()) ? true : false;
passedCutsArrayReverse[1]  = (hcalVeto.passesVeto()) ? true : false;
passedCutsArrayReverse[2]  = (vetoNew.getNStraightTracks() < 3) ? true : false;
passedCutsArrayReverse[3]  = (vetoNew.getStdLayerHit() < 5) ? true : false;
passedCutsArrayReverse[4]  = (vetoNew.getMaxCellDep() < 300) ? true : false;
passedCutsArrayReverse[5]  = (vetoNew.getYStd() < 70) ? true : false;
passedCutsArrayReverse[6]  = (vetoNew.getShowerRMS() < 110) ? true : false;
passedCutsArrayReverse[7]  = (vetoNew.getNReadoutHits() < 70) ? true : false;
passedCutsArrayReverse[8]  = (vetoNew.getEcalBackEnergy() < 250) ? true : false;
passedCutsArrayReverse[9]  = (vetoNew.getSummedTightIso() < 800) ? true : false;
passedCutsArrayReverse[10]  = (vetoNew.getSummedDet() < 3500) ? true : false;
passedCutsArrayReverse[11]  = ((fiducial_ && vetoNew.getFiducial()) || (!fiducial_ && !vetoNew.getFiducial())) ? true : false;

for (size_t i=0;i<sizeof(passedCutsArrayReverse);i++) {
  bool allCutsPassedSoFar = true;
  for (size_t j=0;j<=i;j++) {
    if (!passedCutsArrayReverse[j]) {
      allCutsPassedSoFar = false;
    }
  }
  if (allCutsPassedSoFar) {
    histograms_.fill("Rev_AvgLayerHit", i+1 , vetoNew.getAvgLayerHit() );
    histograms_.fill("Rev_DeepestLayerHit", i+1 , vetoNew.getDeepestLayerHit() );
    histograms_.fill("Rev_EcalBackEnergy", i+1 , vetoNew.getEcalBackEnergy() );
    histograms_.fill("Rev_EpAng", i+1 , vetoNew.getEPAng() );
    histograms_.fill("Rev_EpSep", i+1 , vetoNew.getEPSep() );
    histograms_.fill("Rev_FirstNearPhLayer", i+1 , vetoNew.getFirstNearPhLayer() );
    histograms_.fill("Rev_MaxCellDep", i+1 , vetoNew.getMaxCellDep() );
    histograms_.fill("Rev_NReadoutHits", i+1 , vetoNew.getNReadoutHits() );
    histograms_.fill("Rev_StdLayerHit", i+1 , vetoNew.getStdLayerHit() );
    histograms_.fill("Rev_Straight", i+1 , vetoNew.getNStraightTracks() );
    histograms_.fill("Rev_LinRegNew", i+1 , vetoNew.getNLinRegTracks() );
    histograms_.fill("Rev_SummedDet", i+1 , vetoNew.getSummedDet() );
    histograms_.fill("Rev_SummedTightIso", i+1 , vetoNew.getSummedTightIso() );
    histograms_.fill("Rev_ShowerRMS", i+1 , vetoNew.getShowerRMS() );
    histograms_.fill("Rev_XStd", i+1 , vetoNew.getXStd() );
    histograms_.fill("Rev_YStd", i+1 , vetoNew.getYStd() );
    histograms_.fill("Rev_Hcal_MaxPE", i+1. , hcalMaxPE );
    histograms_.fill("Rev_Hcal_TotalPE", i+1. , hcalTotalPe );
    histograms_.fill("Rev_Hcal_MaxTiming", i+1. , hcalMaxTiming );
    histograms_.fill("Rev_Hcal_MaxSector", i+1. , hcalMaxSector );
  }
}
  
// N-1 plots
// << "      >> Doing N1 plots";
// i=0 is trigger, i=1 is fiducial
 for (size_t i=2;i<sizeof(passedCutsArrayCnC);i++) {
   bool allOtherCutsPassed = true;
   for (size_t j=2;j<sizeof(passedCutsArrayCnC);j++) {
     if (i==j) continue;
     if (!passedCutsArrayCnC[j]) {
       allOtherCutsPassed = false;
         // We found a cut that's not passed, no point in looking into the rest of them
       break;
     }
   }

   if (allOtherCutsPassed && trigResult.passed() && ((fiducial_ && vetoNew.getFiducial()) || (!fiducial_ && !vetoNew.getFiducial()))) {
    if (i==2) histograms_.fill("N1_SummedDet", i+1 , vetoNew.getSummedDet() );
    if (i==3) histograms_.fill("N1_SummedTightIso", i+1 , vetoNew.getSummedTightIso() );
    if (i==4) histograms_.fill("N1_EcalBackEnergy", i+1 , vetoNew.getEcalBackEnergy() );
    if (i==5) histograms_.fill("N1_NReadoutHits", i+1 , vetoNew.getNReadoutHits() );
    if (i==6) histograms_.fill("N1_ShowerRMS", i+1 , vetoNew.getShowerRMS() );
    if (i==7) histograms_.fill("N1_YStd", i+1 , vetoNew.getYStd() );
    if (i==8) histograms_.fill("N1_MaxCellDep", i+1 , vetoNew.getMaxCellDep() );
    if (i==9) histograms_.fill("N1_StdLayerHit", i+1 , vetoNew.getStdLayerHit() );
    if (i==10) histograms_.fill("N1_Straight", i+1 , vetoNew.getNStraightTracks() );
    if (i==11) {
      histograms_.fill("N1_Hcal_MaxPE", i+1. , hcalMaxPE );
      histograms_.fill("N1_Hcal_TotalPE", i+1. , hcalTotalPe );
      histograms_.fill("N1_Hcal_MaxTiming", i+1. , hcalMaxTiming );
      histograms_.fill("N1_Hcal_MaxSector", i+1. , hcalMaxSector );
    }
    
    // histograms_.fill("N1_LinRegNew", i+1 , vetoNew.getNLinRegTracks() );
    // histograms_.fill("N1_EpAng", i+1 , vetoNew.getEPAng() );
    // histograms_.fill("N1_EpSep", i+1 , vetoNew.getEPSep() );
    // histograms_.fill("N1_FirstNearPhLayer", i+1 , vetoNew.getFirstNearPhLayer() );
    // histograms_.fill("N1_XStd", i+1 , vetoNew.getXStd() );
    // histograms_.fill("N1_AvgLayerHit", i+1 , vetoNew.getAvgLayerHit() );
    // histograms_.fill("N1_DeepestLayerHit", i+1 , vetoNew.getDeepestLayerHit() );
   }
 }
}

template <typename T, size_t n>
bool CutBasedDM::passPreselection(T (&passedCutsArray)[n], bool verbose) {
  std::map<int, std::string> namingMap;
  namingMap[0] = "Trigger";
  namingMap[1] = "p_{T}";
  namingMap[2] = "#eta";
  namingMap[3] = "N_{no-L1 pixel hits}";
  namingMap[4] = "f_{valid/all hits}";
  namingMap[5] = "N_{dEdx hits}";
  namingMap[6] = "HighPurity";
  namingMap[7] = "#chi^{2} / N_{dof}";
  namingMap[8] = "d_{z}";
  namingMap[9] = "d_{xy}";
  namingMap[10] = "MiniRelIsoAll";
  namingMap[11] = "MiniRelTkIso";
  namingMap[12] = "E/p";
  namingMap[13] = "#sigma_{p_{T}} / p_{T}^{2}";
  namingMap[14] = "F_{i}";
  
// Return false in the function if a given cut is not passed
  for (size_t i=0;i<sizeof(T) * n;i++) {
    if (passedCutsArray[i]) {
    } else {
      if (verbose) ldmx_log(debug) << "        >> Preselection not passed for the " <<  namingMap[i];
      return false;
    }
  }
}

std::tuple<int, const ldmx::SimParticle *> CutBasedDM::getRecoilEle(
    const std::map<int, ldmx::SimParticle> &particleMap) {
  // The recoil electron is "produced" in the dark brem geneartion
  for (const auto &[trackID, particle] : particleMap) {
    if (particle.getPdgID() == 11 and particle.getProcessType() == ldmx::SimParticle::ProcessType::eDarkBrem) {
      return {trackID, &particle};
    }   
  }

  // // only get here if recoil electron was not "produced" by dark brem
  // //   in this case (bkgd), we interpret the primary electron as also the recoil
  // //   electron
  // ldmx::SimParticle::ProcessType::Primary
  return {1, &(particleMap.at(1))};
}

DECLARE_ANALYZER(CutBasedDM);

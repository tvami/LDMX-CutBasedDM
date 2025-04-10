#include "Framework/EventProcessor.h"
#include "Ecal/Event/EcalVetoResult.h"
#include "Recon/Event/TriggerResult.h"
#include "Hcal/Event/HcalVetoResult.h"
#include "DetDescr/HcalID.h"
#include "SimCore/Event/SimParticle.h"
#include "DetDescr/SimSpecialID.h"
#include "SimCore/Event/SimTrackerHit.h"
#include "Tracking/include/Tracking/Event/Track.h"
#include "Recon/Event/FiducialFlag.h"


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
// v15: Acceptance plot, and add Acceptance to all cutflows
// v16: Add ignore_fiducial_analysis_, move to recoil from tracking
// v17: Running on resim sample, 17b adding HCAL plots, print out for surviving everything
// v18: ldmx-sw v4.2.19, adding CnCwithTracking


class CutBasedDM : public framework::Analyzer {
public:
  CutBasedDM(const std::string& name, framework::Process& p)
  : framework::Analyzer(name, p) {}
  ~CutBasedDM() = default;
  void configure(framework::config::Parameters &ps);
  void setHistLabels(const std::string &name, const std::vector<std::string> &labels);
  void setHistLabelsY(const std::string &name, const std::vector<std::string> &labels);

  void onProcessStart();
  void analyze(const framework::Event& event) final;
  template <typename T, size_t n>
  bool passPreselection(T (&passedCutsArray)[n], bool verbose);
  std::tuple<int, const ldmx::SimParticle *> getRecoilEle(
    const std::map<int, ldmx::SimParticle> &particleMap);
  std::string trigger_collName_;
  std::string trigger_passName_;
  std::string sp_pass_name_;
  std::string track_pass_name_;
  // std::string tagger_track_collection_;
  std::string recoil_track_collection_;
  bool fiducial_analysis_;
  bool ignore_fiducial_analysis_;
  bool ignore_tagger_analysis_;
  bool signal_;
};


void CutBasedDM::configure(framework::config::Parameters &ps) {
  trigger_collName_ = ps.getParameter<std::string>("trigger_name");
  trigger_passName_ = ps.getParameter<std::string>("trigger_pass");
  sp_pass_name_ = ps.getParameter<std::string>("sp_pass_name");
  track_pass_name_ = ps.getParameter<std::string>("track_pass_name","");
  recoil_track_collection_ = ps.getParameter<std::string>("recoil_track_collection","");
  fiducial_analysis_ = ps.getParameter<bool>("fiducial_analysis");
  ignore_fiducial_analysis_ = ps.getParameter<bool>("ignore_fiducial_analysis");
  ignore_tagger_analysis_ = ps.getParameter<bool>("ignore_tagger_analysis",false);
  signal_ = ps.getParameter<bool>("signal", true);

  return;
}

void CutBasedDM::setHistLabels(const std::string &name,
  const std::vector<std::string> &labels) {
    auto histo{histograms_.get(name)};
    for (std::size_t ibin{1}; ibin <= labels.size(); ibin++) {
      histo->GetXaxis()->SetBinLabel(ibin, labels[ibin - 1].c_str());
    }
}

void CutBasedDM::setHistLabelsY(const std::string &name,
  const std::vector<std::string> &labels) {
    auto histo{histograms_.get(name)};
    for (std::size_t ibin{1}; ibin <= labels.size(); ibin++) {
      histo->GetYaxis()->SetBinLabel(ibin, labels[ibin - 1].c_str());
    }
}

void CutBasedDM::onProcessStart(){
  getHistoDirectory();
  
  histograms_.create("Acceptance", "", 6, -0.5, 5.5, "", 90, -450.0, 450.0);
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
   
  histograms_.create("RecoilPT", "", 20, -0.5, 19.5, "Recoil p_{T} [MeV]", 200, 0.0, 1000.0);
  histograms_.create("RecoilPZ", "", 20, -0.5, 19.5, "Recoil p_{Z} [MeV]", 800, -10.0, 8010.0);
  histograms_.create("RecoilP", "", 20, -0.5, 19.5, "Recoil p [MeV]", 2000, 0.0, 10000.0);
  histograms_.create("RecoilPTAtTarget", "", 20, -0.5, 19.5, "Recoil p_{T} @Target [MeV]", 400, 0.0, 4000.0);
  histograms_.create("RecoilPZAtTarget", "", 20, -0.5, 19.5, "Recoil p_{Z} @Target [MeV]", 800, -10.0, 8010.0);
  histograms_.create("RecoilPAtTarget", "", 20, -0.5, 19.5, "Recoil p @Target [MeV]", 800, 0.0, 10000.0);
  histograms_.create("RecoilTheta", "", 20, -0.5, 19.5, "Recoil theta @Target", 90, 0.0, 90.0);
  histograms_.create("RecoilPhi", "", 20, -0.5, 19.5, "Recoil phi @Target", 360, -180.0, 180.0);


  histograms_.create("Hcal_MaxPE", "", 20, -0.5, 19.5, "HCAL max photo-electron hits", 65, -0.5, 64.5);
  histograms_.create("Hcal_MaxPE_Extended", "", 20, -0.5, 19.5, "HCAL max photo-electron hits", 120, -0.5, 600.5);
  histograms_.create("Hcal_TotalPE", "", 20, -0.5, 19.5, "HCAL total photo-electron hits", 100, -0.5, 200.5);
  histograms_.create("Hcal_TotalPE_AboveMax8PE", "", 20, -0.5, 19.5, "HCAL total photo-electron hits (for MaxPE > 8)", 200, -0.5, 400.5);
  histograms_.create("Hcal_MaxTiming", "", 20, -0.5, 19.5, "HCAL timing of the max PE hit", 35, 0.0, 35.0);
  histograms_.create("Hcal_MaxSector", "", 20, -0.5, 19.5, "", 5, -0.5, 4.5);

  histograms_.create("BDTDiscrVsHcalPE_PreS", "HCAL PE", 500, -0.5, 500.5, "BDT discriminating score", 100, 0.0, 1.0);
  histograms_.create("BDTDiscrLogVsHcalPE_PreS", "HCAL PE", 500, -0.5, 500.5, "-log(1-BDT discriminating score)", 100, 0.0, 5.0);
  histograms_.create("BDTDiscrVsHcalPE_PostS", "HCAL PE", 500, -0.5, 500.5, "BDT discriminating score", 100, 0.0, 1.0);
  histograms_.create("BDTDiscrLogVsHcalPE_PostS", "HCAL PE", 500, -0.5, 500.5, "-log(1-BDT discriminating score)", 100, 0.0, 5.0);

  histograms_.create("AltCutFlow_RecoilX", "", 18, -0.5, 17.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
  histograms_.create("StdCutFlow_RecoilX", "", 18, -0.5, 17.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
  histograms_.create("StdCutFlowWithTracking_RecoilX", "", 20, -0.5, 19.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
  histograms_.create("BDTCutFlow_RecoilX", "", 18, -0.5, 17.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
  // histograms_.create("LinRegCutFlow_RecoilX", "",  18, -0.5, 17.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
  // histograms_.create("LinRegCutFlowHcal_RecoilX", "",  18, -0.5, 17.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
 
  histograms_.create("TrackingCutFlow_RecoilX", "", 18, -0.5, 17.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
  histograms_.create("Tracking_TaggerP", "", 18, -0.5, 17.5, "Tagger p [MeV]", 800, 0.0, 10000.0);
  histograms_.create("Tracking_RecoilN", "", 18, -0.5, 17.5, "N_{recoil}", 10, -0.5, 9.5);
  histograms_.create("Tracking_RecoilP", "", 18, -0.5, 17.5, "Recoil p [MeV]", 2000, 0.0, 10000.0);
  histograms_.create("Tracking_RecoilPt", "", 18, -0.5, 17.5, "Recoil p_{T} [MeV]", 200, 0.0, 1000.0);
  histograms_.create("Tracking_RecoilD0", "", 18, -0.5, 17.5, "d_{0} [mm]", 100, -50.0, 50.0);
  histograms_.create("Tracking_RecoilZ0", "", 18, -0.5, 17.5, "z_{0} [mm]", 100, -50.0, 50.0);
  histograms_.create("TrackingCutFlowHcal_RecoilX", "", 18, -0.5, 17.5, "RecoilX @Ecal [mm]", 90, -450.0, 450.0);
  histograms_.create("TrackingHcal_TaggerP", "", 18, -0.5, 17.5, "Tagger p [MeV]", 2000, 0.0, 10000.0);
  histograms_.create("TrackingHcal_RecoilN", "", 18, -0.5, 17.5, "N_{recoil}", 10, -0.5, 9.5);
  histograms_.create("TrackingHcal_RecoilD0", "", 18, -0.5, 17.5, "d_{0} [mm]", 100, -50.0, 50.0);
  histograms_.create("TrackingHcal_RecoilZ0", "", 18, -0.5, 17.5, "z_{0} [mm]", 100, -50.0, 50.0);

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
  // histograms_.create("Rev_LinRegNew", "", 20, -0.5, 19.5, "Linear regression tracks", 15, -0.5, 14.5);
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
  // histograms_.create("N1_EcalBackEnergy", "", 20, -0.5, 19.5, "Ecal back energy [MeV]", 100, 0.0, 3000.0);
  // histograms_.create("N1_MaxCellDep", "", 20, -0.5, 19.5, "Max cell deposition [MeV]", 100, 0.0, 800.0);
  // histograms_.create("N1_NReadoutHits", "", 20, -0.5, 19.5, "#Readout hits", 150, -0.5, 149.5);
  // histograms_.create("N1_StdLayerHit", "", 20, -0.5, 19.5, "Std of hit layers", 70, -0.5, 34.5);
  // histograms_.create("N1_Straight", "", 20, -0.5, 19.5, "Straight tracks", 15, -0.5, 14.5);
  // histograms_.create("N1_SummedDet", "", 20, -0.5, 19.5, "Summed ECAL energy [MeV]", 100, 0.0, 10000.0);
  // histograms_.create("N1_SummedTightIso", "", 20, -0.5, 19.5, "Summed ECAL energy with tight iso [MeV]", 100, 0.0, 10000.0);
  // histograms_.create("N1_ShowerRMS", "", 20, -0.5, 19.5, "Shower RMS [mm]", 100, 0.0, 250.0);
  // histograms_.create("N1_YStd", "", 20, -0.5, 19.5, "Shower RMS_{Y} [mm]", 100, 0.0, 250.0);
  // histograms_.create("N1_Hcal_MaxPE", "", 20, -0.5, 19.5, "HCAL max photo-electron hits", 65, -0.5, 64.5);
  // histograms_.create("N1_Hcal_TotalPE", "", 20, -0.5, 19.5, "HCAL total photo-electron hits", 100, -0.5, 200.5);
  // histograms_.create("N1_Hcal_MaxTiming", "", 20, -0.5, 19.5, "HCAL timing of the max PE hit", 35, 0.0, 35.0);
  // histograms_.create("N1_Hcal_MaxSector", "", 20, -0.5, 19.5, "", 5, -0.5, 4.5);

  // histograms_.create("N1_EcalBackEnergy", "", 20, -0.5, 19.5, "Ecal back energy [MeV]", 100, 0.0, 3000.0);

  std::vector<std::string> labels = {
    "All / Acceptance",
    "Fiducial",              // 1
    "Triggerred",            // 2
    "E_{sum} < 3500",        // 3
    "E_{SumTight} < 800",        // 4
    "E_{back} < 250",             // 5
    "N_{hits} < 70",              // 6
    "RMS_{shower} < 110",         // 7
    "RMS_{shower,Y} < 70",        // 8
    "E_{cell,max} < 300",         // 9
    "RMS_{Layer,hit} < 5",        // 10
    "N_{straight} < 3",           // 11
    "PE_{HCal,max} < 8",            // 12
    };
    
  if (!fiducial_analysis_) labels.at(1) = "Non-fiducial";

  setHistLabels("AvgLayerHit", labels);
  setHistLabels("DeepestLayerHit", labels);
  setHistLabels("EcalBackEnergy", labels);
  setHistLabels("EpAng", labels);
  setHistLabels("EpSep", labels);
  setHistLabels("FirstNearPhLayer", labels);
  setHistLabels("MaxCellDep", labels);
  setHistLabels("NReadoutHits", labels);
  setHistLabels("StdLayerHit", labels);
  setHistLabels("Straight", labels);
  setHistLabels("LinRegNew", labels);
  setHistLabels("SummedDet", labels);
  setHistLabels("SummedTightIso", labels);
  setHistLabels("ShowerRMS", labels);
  setHistLabels("XStd", labels);
  setHistLabels("YStd", labels);
  setHistLabels("BDTDiscr", labels);
  setHistLabels("BDTDiscrLog", labels);
  setHistLabels("StdCutFlow_RecoilX", labels);

  setHistLabels("RecoilX", labels);
  setHistLabels("RecoilPT", labels);
  setHistLabels("RecoilPZ", labels);
  setHistLabels("RecoilP", labels);
  setHistLabels("RecoilXAtTarget", labels);
  setHistLabels("RecoilPTAtTarget", labels);
  setHistLabels("RecoilPZAtTarget", labels);
  setHistLabels("RecoilPAtTarget", labels);
  setHistLabels("RecoilTheta", labels);
  setHistLabels("RecoilPhi", labels);
  setHistLabels("Hcal_MaxPE", labels);
  setHistLabels("Hcal_MaxPE_Extended", labels);
  setHistLabels("Hcal_TotalPE", labels);
  setHistLabels("Hcal_TotalPE_AboveMax8PE", labels);
  setHistLabels("Hcal_MaxTiming", labels);
  setHistLabels("Hcal_MaxSector", labels);

  // setHistLabels("N1_EcalBackEnergy", labels);
  // setHistLabels("N1_MaxCellDep", labels);
  // setHistLabels("N1_NReadoutHits", labels);
  // setHistLabels("N1_StdLayerHit", labels);
  // setHistLabels("N1_Straight", labels);
  // setHistLabels("N1_SummedDet", labels);
  // setHistLabels("N1_SummedTightIso", labels);
  // setHistLabels("N1_ShowerRMS", labels);
  // setHistLabels("N1_YStd", labels);
  // setHistLabels("N1_Hcal_MaxPE", labels);
  // setHistLabels("N1_Hcal_MaxTiming", labels);
  // setHistLabels("N1_Hcal_MaxSector", labels);

  std::vector<std::string> labelsWithTracking = {
    "All / Acceptance",
    "Fiducial",              // 1
    "Triggerred",            // 2
    "p_{tagger} > 5600",     //
    "N_{recoil} = 1",
    "|d_{0}| < 10",
    "|z_{0}| < 40",
    "E_{sum} < 3500",        //
    "E_{SumTight} < 800",        // 4
    "E_{back} < 250",             // 5
    "N_{hits} < 70",              // 6
    "RMS_{shower} < 110",         // 7
    "RMS_{shower,Y} < 70",        // 8
    "E_{cell,max} < 300",         // 9
    "RMS_{Layer,hit} < 5",        // 10
    "N_{straight} < 3",           // 11
    "PE_{HCal,max} < 8",            // 12
    "N_{straight} = 0",           //
    };
    
  if (!fiducial_analysis_) labelsWithTracking.at(1) = "Non-fiducial";

  setHistLabels("StdCutFlowWithTracking_RecoilX", labelsWithTracking);
  
  std::vector<std::string> labels_Rev = {
    "All / Acceptance",      // 0
    "Fiducial",              // 1
    "Triggerred",            // 2
    "PE_{HCal,max} < 8",           // 13
    "N_{straight} < 3",           // 11
    "RMS_{Layer,hit} < 5",        // 10
    "E_{cell,max} < 300",         // 9
    "RMS_{shower,Y} < 70",        // 8
    "RMS_{shower} < 110",         // 7
    "N_{hits} < 70",              // 6
    "E_{back} < 250",             // 5
    "E_{SumTight} < 800",        // 4
    "E_{sum} < 3500",       // 3
    };

  if (!fiducial_analysis_) labels_Rev.at(1) = "Non-fiducial";
  setHistLabels("Rev_AvgLayerHit", labels_Rev);
  setHistLabels("Rev_DeepestLayerHit", labels_Rev);
  setHistLabels("Rev_EcalBackEnergy", labels_Rev);
  setHistLabels("Rev_EpAng", labels_Rev);
  setHistLabels("Rev_EpSep", labels_Rev);
  setHistLabels("Rev_FirstNearPhLayer", labels_Rev);
  setHistLabels("Rev_MaxCellDep", labels_Rev);
  setHistLabels("Rev_NReadoutHits", labels_Rev);
  setHistLabels("Rev_StdLayerHit", labels_Rev);
  setHistLabels("Rev_Straight", labels_Rev);
  // setHistLabels("Rev_LinRegNew", labels_Rev);
  setHistLabels("Rev_SummedDet", labels_Rev);
  setHistLabels("Rev_SummedTightIso", labels_Rev);
  setHistLabels("Rev_ShowerRMS", labels_Rev);
  setHistLabels("Rev_XStd", labels_Rev);
  setHistLabels("Rev_YStd", labels_Rev);
  setHistLabels("Rev_Hcal_MaxPE", labels_Rev);
  setHistLabels("Rev_Hcal_MaxTiming", labels_Rev);
  setHistLabels("Rev_Hcal_MaxSector", labels_Rev);
 

  // enum HcalSection { BACK = 0, TOP = 1, BOTTOM = 2, RIGHT = 3, LEFT = 4 };
  std::vector<std::string> labels_HCALsector = {
    "HCAL BACK",         // 0
    "HCAL TOP",          // 1
    "HCAL BOTTOM",       // 2
    "HCAL RIGHT",        // 3
    "HCAL LEFT",         // 4
    };

  setHistLabelsY("Hcal_MaxSector", labels_HCALsector);
  // setHistLabelsY("N1_Hcal_MaxSector", labels_HCALsector);
  setHistLabelsY("Rev_Hcal_MaxSector", labels_HCALsector);

  std::vector<std::string> labels_accpt = {
      "All",         // 0
      "Min energy",  // 1
      "Min tk hits", // 2
      "Ecal hit",    // 3
      "Hcal hit",    // 4
      "Acceptance"   // 5
      };
  setHistLabels("Acceptance", labels_accpt);

  // enum HcalSection { BACK = 0, TOP = 1, BOTTOM = 2, RIGHT = 3, LEFT = 4 };
  std::vector<std::string> labels_AltCutFlow = {
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
  "PE_{HCal,max} < 8",            // 12
  };
  if (!fiducial_analysis_) labels_AltCutFlow.at(2) = "Non-fiducial";

  setHistLabels("AltCutFlow_RecoilX",labels_AltCutFlow);


  // CutFlow labels for BDT
  std::vector<std::string> labels_BDTCutFlow = {
  "All / Acceptance",      // 0
  "Fiducial",              // 1
  "Triggerred",            // 2
  "Ecal BDT",              // 3
  "N_{straight} < 3",      // 4
  "PE_{HCal,max} < 8",       // 5
  "N_{straight} = 0",      // 6
  "Angle_{e,ph} > 3.",     // 7
  };

  if (!fiducial_analysis_) labels_BDTCutFlow.at(1) = "Non-fiducial";
  setHistLabels("BDTCutFlow_RecoilX",labels_BDTCutFlow);

  // // CutFlow labels for BDT with new lin-reg
  // std::vector<std::string> labels_LinRegCutFlow = {
  // "All",
  // "Triggerred",            // 1
  // "Fiducial",              // 2
  // "Ecal BDT",            // 3
  // "N_{straight} < 3",      // 4
  // "PE_{HCal,max} < 8",       // 5
  // "N_{straight} = 0",      // 6
  // "N_{lin-reg} = 0",      // 7
  // "Angle_{e,ph} > 3.",     // 8
  // ""};

  // if (!fiducial_analysis_) labels_LinRegCutFlow.at(2) = "Non-fiducial";
  // setHistLabels("LinRegCutFlow_RecoilX",labels_LinRegCutFlow);

  // // CutFlow labels for BDT with new lin-reg
  // std::vector<std::string> labels_LinRegCutFlowHcal = {
  // "All",
  // "Triggerred",            // 1
  // "PE_{HCal,max} < 8",       // 2
  // "Fiducial",              // 3
  // "Ecal BDT",            // 4
  // "N_{straight} = 0",      // 5
  // "N_{lin-reg} = 0",      // 6
  // "Angle_{e,ph} > 3.",     // 7
  // ""};

  // if (!fiducial_analysis_) labels_LinRegCutFlowHcal.at(3) = "Non-fiducial";
  // setHistLabels("LinRegCutFlowHcal_RecoilX",labels_LinRegCutFlowHcal);

  // CutFlow labels for BDT with tracking
  std::vector<std::string> labels_TrackingCutFlow = {
  "All / Acceptance",      // 0
  "Fiducial",              // 1
  "Triggerred",            // 2
  "p_{tagger} > 5600",     // 3
  "N_{recoil} = 1",
  "|d_{0}| < 10",
  "|z_{0}| < 40",
  "Ecal BDT",
  "N_{straight} < 3",
  "PE_{HCal,max} < 8",
  "N_{straight} = 0",
  "Angle_{e,ph} > 3.",     // 11
  ""};

  if (!fiducial_analysis_) labels_TrackingCutFlow.at(6) = "Non-fiducial";

  setHistLabels("TrackingCutFlow_RecoilX",labels_TrackingCutFlow);
  setHistLabels("Tracking_TaggerP",labels_TrackingCutFlow);
  setHistLabels("Tracking_RecoilN",labels_TrackingCutFlow);
  setHistLabels("Tracking_RecoilP",labels_TrackingCutFlow);
  setHistLabels("Tracking_RecoilPt",labels_TrackingCutFlow);
  setHistLabels("Tracking_RecoilD0",labels_TrackingCutFlow);
  setHistLabels("Tracking_RecoilZ0",labels_TrackingCutFlow);

  // CutFlow labels for BDT with tracking with Hcal first
  std::vector<std::string> labels_TrackingCutFlowHcal = {
  "All / Acceptance",
  "Fiducial",              // 1
  "Triggerred",            // 2
  "PE_{HCal,max} < 8",       // 3
  "Ecal BDT",            // 4
  "N_{straight} = 0",      // 5
  "Angle_{e,ph} > 3.",     // 6
  "p_{tagger} > 5600",     // 7
  "N_{recoil} = 1",  		 // 8
  "|d_{0}| < 10",  		 // 9
  "|z_{0}| < 40",  		 // 10
  ""};
  if (!fiducial_analysis_) labels_TrackingCutFlowHcal.at(1) = "Non-fiducial";
  setHistLabels("TrackingCutFlowHcal_RecoilX",labels_TrackingCutFlowHcal);
  setHistLabels("TrackingHcal_TaggerP",labels_TrackingCutFlowHcal);
  setHistLabels("TrackingHcal_RecoilN",labels_TrackingCutFlowHcal);
  setHistLabels("TrackingHcal_RecoilD0",labels_TrackingCutFlowHcal);
  setHistLabels("TrackingHcal_RecoilZ0",labels_TrackingCutFlowHcal);

} 

void CutBasedDM::analyze(const framework::Event& event) {
  //std::cout << " ---------------------------------------------" << std::endl;
  auto vetoNew{event.getObject<ldmx::EcalVetoResult>("EcalVetoNew","")};
  auto trigResult{event.getObject<ldmx::TriggerResult>(trigger_collName_, trigger_passName_)};
  auto hcalVeto{event.getObject<ldmx::HcalVetoResult>("HcalVeto","cutbased")};
  auto hcalRecHits{event.getCollection<ldmx::HcalHit>("HcalRecHits", "")};
  auto targetSpHits{event.getCollection<ldmx::SimTrackerHit>("TargetScoringPlaneHits",sp_pass_name_)};
  auto recoilTrackCollection{event.getCollection<ldmx::Track>(recoil_track_collection_)};
  std::vector<ldmx::Track> taggerTrackCollection;
  if (!ignore_tagger_analysis_) {
    taggerTrackCollection = event.getCollection<ldmx::Track>("TaggerTracks","cutbased");
  }

  bool acceptance{true};
  int fiducial_analysis_flag{-1};
  if (signal_) {
    auto acceptanceChecks{event.getObject<ldmx::FiducialFlag>("RecoilTruthFiducialFlags")};
    acceptance =  acceptanceChecks.isFiducial();
    fiducial_analysis_flag = acceptanceChecks.getFiducialFlag();
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
  //std::cout << " pYAtTarget " <<  pYAtTarget << " totMomAtTarget " << totMomAtTarget << std::endl ;
  auto phiEleAtTarget =   (180/M_PI)*std::acos(pYAtTarget/totMomAtTarget)-90.;
  auto thetaEleAtTarget = (180/M_PI)*std::acos(pZAtTarget/totMomAtTarget);
  // //std::cout << " phiEleAtTarget " << phiEleAtTarget << " phiEle " << phiEle
  //  //std::cout << " thetaEleAtTarget " << thetaEleAtTarget << " phiEleAtTarget " << phiEleAtTarget << std::endl;

  // Take recoil momentum from ECAL SP
  // pT2 = vetoNew.getRecoilMomentum()[0]*vetoNew.getRecoilMomentum()[0] + vetoNew.getRecoilMomentum()[1]*vetoNew.getRecoilMomentum()[1];
  // pZ =  vetoNew.getRecoilMomentum()[2];


  // HCAL veto calc
  float hcalMaxPE{-9999};
  float hcalMaxTiming{-9999};
  int hcalMaxSector{-1};
  float hcalTotalPe{0};
  float hcalTotalPeAbove8PE{0};
  
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
  //std::cout << " Num of hcalRecHits = " << hcalRecHits.size() << std::endl;
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
    if (pe > 8) {
      hcalTotalPeAbove8PE  += pe;
    }

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


  // std::cout << "Fiducial = " << vetoNew.getFiducial() << std::endl;

  // CutFlow here
  bool passedCutsArrayCnC[13];
  //std::cout << " CnC cutflow = " << std::endl;
  std::fill(std::begin(passedCutsArrayCnC), std::end(passedCutsArrayCnC),false);
  passedCutsArrayCnC[0]  = (acceptance) ? true : false;
  passedCutsArrayCnC[1]  = (ignore_fiducial_analysis_ || (fiducial_analysis_ && vetoNew.getFiducial()) || (!fiducial_analysis_ && !vetoNew.getFiducial())) ? true : false;
  passedCutsArrayCnC[2]  = (trigResult.passed()) ? true : false;
  passedCutsArrayCnC[3]  = (vetoNew.getSummedDet() < 3500) ? true : false;
  passedCutsArrayCnC[4]  = (vetoNew.getSummedTightIso() < 800) ? true : false;
  passedCutsArrayCnC[5]  = (vetoNew.getEcalBackEnergy() < 250) ? true : false;
  passedCutsArrayCnC[6]  = (vetoNew.getNReadoutHits() < 70) ? true : false;
  passedCutsArrayCnC[7]  = (vetoNew.getShowerRMS() < 110) ? true : false;
  passedCutsArrayCnC[8]  = (vetoNew.getYStd() < 70) ? true : false;
  passedCutsArrayCnC[9]  = (vetoNew.getMaxCellDep() < 300) ? true : false;
  passedCutsArrayCnC[10]  = (vetoNew.getStdLayerHit() < 5) ? true : false;
  passedCutsArrayCnC[11]  = (vetoNew.getNStraightTracks() < 3) ? true : false;
  passedCutsArrayCnC[12]  = (hcalVeto.passesVeto()) ? true : false;

  // Fill histograms

  bool has_min_energy       = fiducial_analysis_flag & (1 << 0);
  bool has_min_tracker_hits = fiducial_analysis_flag & (1 << 1);
  bool has_ecal_hit         = fiducial_analysis_flag & (1 << 2);
  bool has_hcal_hit         = fiducial_analysis_flag & (1 << 3);
  if (fiducial_analysis_flag > 0) {
    histograms_.fill("Acceptance", 0. , vetoNew.getRecoilX() );
    if (has_min_energy) histograms_.fill("Acceptance", 1. , vetoNew.getRecoilX() );
    if (has_min_tracker_hits) histograms_.fill("Acceptance", 2. , vetoNew.getRecoilX() );
    if (has_ecal_hit) histograms_.fill("Acceptance", 3. , vetoNew.getRecoilX() );
    if (has_hcal_hit) histograms_.fill("Acceptance", 4. , vetoNew.getRecoilX() );
    if (acceptance) histograms_.fill("Acceptance", 5. , vetoNew.getRecoilX() );
  }
  

  for (size_t i=0;i<sizeof(passedCutsArrayCnC);i++) {
    bool allCutsPassedSoFar = true;
    for (size_t j=0;j<=i;j++) {
      if (!passedCutsArrayCnC[j]) {
        allCutsPassedSoFar = false;
        break;
      }
    }
    if (allCutsPassedSoFar) {
      // //std::cout 
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
      // //std::cout << "hcalMaxPE = " <<  hcalMaxPE << " hcal total" << hcalTotalPe <<  " maxTime = " << hcalMaxTiming << " where = " << hcalMaxSector
      // << std::endl;

      histograms_.fill("RecoilX", i, vetoNew.getRecoilX() );
      histograms_.fill("AvgLayerHit", i, vetoNew.getAvgLayerHit() );
      histograms_.fill("DeepestLayerHit", i, vetoNew.getDeepestLayerHit() );
      histograms_.fill("EcalBackEnergy", i, vetoNew.getEcalBackEnergy() );
      histograms_.fill("EpAng", i, vetoNew.getEPAng() );
      histograms_.fill("EpSep", i, vetoNew.getEPSep() );
      histograms_.fill("FirstNearPhLayer", i, vetoNew.getFirstNearPhLayer() );
      histograms_.fill("MaxCellDep", i, vetoNew.getMaxCellDep() );
      histograms_.fill("NReadoutHits", i, vetoNew.getNReadoutHits() );
      histograms_.fill("StdLayerHit", i, vetoNew.getStdLayerHit() );
      histograms_.fill("Straight", i, vetoNew.getNStraightTracks() );
      histograms_.fill("LinRegNew", i, vetoNew.getNLinRegTracks() );
      histograms_.fill("SummedDet", i, vetoNew.getSummedDet() );
      histograms_.fill("SummedTightIso", i, vetoNew.getSummedTightIso() );
      histograms_.fill("ShowerRMS", i, vetoNew.getShowerRMS() );
      histograms_.fill("XStd", i, vetoNew.getXStd() );
      histograms_.fill("YStd", i, vetoNew.getYStd() );
      histograms_.fill("BDTDiscr", i, vetoNew.getDisc() );
      histograms_.fill("BDTDiscrLog", i, -log(1-vetoNew.getDisc()) );
      histograms_.fill("StdCutFlow_RecoilX", i, vetoNew.getRecoilX() );
      histograms_.fill("RecoilPT", i, pT );
      histograms_.fill("RecoilPZ", i, pZ );
      histograms_.fill("RecoilP", i, totMom );
      histograms_.fill("RecoilXAtTarget", i,XAtTarget );
      histograms_.fill("RecoilPTAtTarget", i,pTAtTarget );
      histograms_.fill("RecoilPZAtTarget", i, pZAtTarget );
      histograms_.fill("RecoilPAtTarget", i, totMomAtTarget );
      histograms_.fill("RecoilTheta", i, thetaEleAtTarget );
      histograms_.fill("RecoilPhi", i, phiEleAtTarget );
      histograms_.fill("Hcal_MaxPE", i, hcalMaxPE );
      histograms_.fill("Hcal_MaxPE_Extended", i, hcalMaxPE );
      histograms_.fill("Hcal_TotalPE", i, hcalTotalPe );
      histograms_.fill("Hcal_TotalPE_AboveMax8PE", i, hcalTotalPeAbove8PE );
      histograms_.fill("Hcal_MaxTiming", i, hcalMaxTiming );
      histograms_.fill("Hcal_MaxSector", i, hcalMaxSector );

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
  bool passedCutsArrayAlt[13];
  //std::cout << " Alternative cutFlow here " << std::endl;
  std::fill(std::begin(passedCutsArrayAlt), std::end(passedCutsArrayAlt),false);
  passedCutsArrayAlt[0]  = (acceptance) ? true : false;
  passedCutsArrayAlt[1]  = (ignore_fiducial_analysis_ || (fiducial_analysis_ && vetoNew.getFiducial()) || (!fiducial_analysis_ && !vetoNew.getFiducial())) ? true : false;
  passedCutsArrayAlt[2]  = (trigResult.passed()) ? true : false;
  passedCutsArrayAlt[3]  = (vetoNew.getSummedDet() < 3500) ? true : false;
  passedCutsArrayAlt[4]  = (vetoNew.getSummedTightIso() < 800) ? true : false;
  passedCutsArrayAlt[5]  = (vetoNew.getEcalBackEnergy() < 250) ? true : false;
  passedCutsArrayAlt[6]  = (vetoNew.getNReadoutHits() < 70) ? true : false;
  passedCutsArrayAlt[7]  = (vetoNew.getShowerRMS() < 110) ? true : false;
  passedCutsArrayAlt[8]  = (vetoNew.getYStd() < 70) ? true : false;
  passedCutsArrayAlt[9]  = (vetoNew.getMaxCellDep() < 300) ? true : false;
  passedCutsArrayAlt[10]  = (vetoNew.getStdLayerHit() < 5) ? true : false;
  passedCutsArrayAlt[11]  = (vetoNew.getNStraightTracks() < 3) ? true : false;
  passedCutsArrayAlt[12]  = (hcalVeto.passesVeto()) ? true : false;

  for (size_t i=0;i<sizeof(passedCutsArrayAlt);i++) {
    bool allCutsPassedSoFar = true;
    for (size_t j=0;j<=i;j++) {
      if (!passedCutsArrayAlt[j]) {
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
      //   << " hcalVeto = " << hcalVeto.passesVeto() << std::endl;
      // std::cout << "hcalMaxPE = " <<  hcalMaxPE << " hcal total" << hcalTotalPe <<  " maxTime = " << hcalMaxTiming << " where = " << hcalMaxSector << std::endl;
      histograms_.fill("AltCutFlow_RecoilX", i, vetoNew.getRecoilX() );
    }
  }

  // BDT based cutFlow here
  bool passedCutsArrayBDT[8];
  //std::cout << " BDT cutFlow here " << std::endl;
  std::fill(std::begin(passedCutsArrayBDT), std::end(passedCutsArrayBDT),false);
  passedCutsArrayBDT[0]  = (acceptance) ? true : false;
  passedCutsArrayBDT[1]  = (ignore_fiducial_analysis_ || (fiducial_analysis_ && vetoNew.getFiducial()) || (!fiducial_analysis_ && !vetoNew.getFiducial())) ? true : false;
  passedCutsArrayBDT[2]  = (trigResult.passed()) ? true : false;
  passedCutsArrayBDT[3]  = (vetoNew.getDisc() > 0.99741) ? true : false;
  passedCutsArrayBDT[4]  = (vetoNew.getNStraightTracks() < 3) ? true : false;
  passedCutsArrayBDT[5]  = (hcalVeto.passesVeto()) ? true : false;
  passedCutsArrayBDT[6]  = (vetoNew.getNStraightTracks() == 0) ? true : false;
  if (fiducial_analysis_) {
    passedCutsArrayBDT[7]  = (vetoNew.getEPAng() > 3)  ? true : false;
  } else {
    passedCutsArrayBDT[7]  = true;
  }


  for (size_t i=0;i<sizeof(passedCutsArrayBDT);i++) {
    bool allCutsPassedSoFar = true;
    for (size_t j=0;j<=i;j++) {
      if (!passedCutsArrayBDT[j]) {
        allCutsPassedSoFar = false;
        break;
      }
    }
    if (allCutsPassedSoFar) {
      // std::cout << " vetoNew.getEPAng() = " << vetoNew.getEPAng() << "  i = " << i << std::endl;
      histograms_.fill("BDTCutFlow_RecoilX", i, vetoNew.getRecoilX() );
    }
  }

  // // BDT based cutFlow here with linreg
  // bool passedCutsArrayLinReg[8];
  // std::fill(std::begin(passedCutsArrayLinReg), std::end(passedCutsArrayLinReg),false);
  // passedCutsArrayLinReg[0]  = (trigResult.passed()) ? true : false;
  // passedCutsArrayLinReg[1]  = ((fiducial_analysis_ && vetoNew.getFiducial()) || (!fiducial_analysis_ && !vetoNew.getFiducial())) ? true : false;
  // passedCutsArrayLinReg[2]  = (vetoNew.getDisc() > 0.99741) ? true : false;
  // passedCutsArrayLinReg[3]  = (vetoNew.getNStraightTracks() < 3) ? true : false;
  // passedCutsArrayLinReg[4]  = (hcalVeto.passesVeto()) ? true : false;
  // passedCutsArrayLinReg[5]  = (vetoNew.getNStraightTracks() == 0) ? true : false;
  // passedCutsArrayLinReg[6]  = (vetoNew.getNLinRegTracks() == 0) ? true : false;
  // passedCutsArrayLinReg[7]  = ((vetoNew.getEPAng() > 3) && (fiducial_analysis_ && vetoNew.getEPAng()  < 999) || (!fiducial_analysis_ )) ? true : false;

  // for (size_t i=0;i<sizeof(passedCutsArrayLinReg);i++) {
  //   bool allCutsPassedSoFar = true;
  //   for (size_t j=0;j<=i;j++) {
  //     if (!passedCutsArrayLinReg[j]) {
  //       allCutsPassedSoFar = false;
  //       break;
  //     }
  //   }
  //   if (allCutsPassedSoFar) {
  //     histograms_.fill("LinRegCutFlow_RecoilX", i, vetoNew.getRecoilX() );
  //   }
  // }

  //   // BDT based cutFlow here with linreg, starting with Hcal
  // bool passedCutsArrayLinRegHcal[6];
  // std::fill(std::begin(passedCutsArrayLinRegHcal), std::end(passedCutsArrayLinRegHcal),false);
  // passedCutsArrayLinRegHcal[0]  = (trigResult.passed()) ? true : false;
  // passedCutsArrayLinRegHcal[1]  = (hcalVeto.passesVeto()) ? true : false;
  // passedCutsArrayLinRegHcal[2]  = ((fiducial_analysis_ && vetoNew.getFiducial()) || (!fiducial_analysis_ && !vetoNew.getFiducial())) ? true : false;
  // passedCutsArrayLinRegHcal[3]  = (vetoNew.getDisc() > 0.99741) ? true : false;
  // passedCutsArrayLinRegHcal[4]  = (vetoNew.getNStraightTracks() == 0) ? true : false;
  // passedCutsArrayLinRegHcal[5]  = (vetoNew.getNLinRegTracks() == 0) ? true : false;


  // for (size_t i=0;i<sizeof(passedCutsArrayLinRegHcal);i++) {
  //   bool allCutsPassedSoFar = true;
  //   for (size_t j=0;j<=i;j++) {
  //     if (!passedCutsArrayLinRegHcal[j]) {
  //       allCutsPassedSoFar = false;
  //       break;
  //     }
  //   }
  //   if (allCutsPassedSoFar) {

  //     histograms_.fill("LinRegCutFlowHcal_RecoilX", i, vetoNew.getRecoilX() );
  //   }
  // }
  // --------------------------------------------------------------------------
  // Calculate tracking variables if tracking is available
  //std::cout << " Tracking variables = " << std::endl;
  float taggerP{0.0}; // Make sure this is in MeV!!
  // Start with tagger tracks
  auto taggerN = taggerTrackCollection.size();
  // std::cout << " taggerN = " << taggerN << std::endl;
  if (taggerN == 1) {
    for (const auto trk : taggerTrackCollection) {
      auto QoP = trk.getQoP();
      taggerP = 1000. / std::abs(QoP);
    }
  }

  // Recoil tracks now
  float recoilP{0.0}; // Make sure this is in MeV!!
  float recoilPt{0.0}; // Make sure this is in MeV!!
  float recoilD0{-9999.};
  float recoilZ0{-9999.};
  auto recoilN = recoilTrackCollection.size();
  //std::cout << " recoilN = " << recoilN << std::endl;
  if (recoilN == 1) {
    for (const auto trk : recoilTrackCollection) {
      recoilD0 = trk.getD0();
      recoilZ0 = trk.getZ0();
      auto QoP = trk.getQoP();
      recoilP = 1000. / std::abs(QoP);
      auto trk_mom = trk.getMomentum();
      recoilPt = 1000 * std::sqrt(trk_mom[1] * trk_mom[1] + trk_mom[2] * trk_mom[2]);
    }
  }

  // --------------------------------------------------------------------------
  // CnC based cutFlow with tracking
  // CutFlow here
  bool passedCutsArrayCnCWithTracking[18];
  //std::cout << " CnC cutflow = " << std::endl;
  std::fill(std::begin(passedCutsArrayCnCWithTracking), std::end(passedCutsArrayCnCWithTracking),false);
  passedCutsArrayCnCWithTracking[0]  = (acceptance) ? true : false;
  passedCutsArrayCnCWithTracking[1]  = (ignore_fiducial_analysis_ || (fiducial_analysis_ && vetoNew.getFiducial()) || (!fiducial_analysis_ && !vetoNew.getFiducial())) ? true : false;
  passedCutsArrayCnCWithTracking[2]  = (trigResult.passed()) ? true : false;
  passedCutsArrayCnCWithTracking[3]  = (ignore_tagger_analysis_ || (taggerP > 5600)) ? true : false;
  passedCutsArrayCnCWithTracking[4]  = (recoilN == 1) ? true : false;
  passedCutsArrayCnCWithTracking[5]  = (std::abs(recoilD0) < 10.) ? true : false;
  passedCutsArrayCnCWithTracking[6]  = (std::abs(recoilZ0) < 40.) ? true : false; 
  passedCutsArrayCnCWithTracking[7]  = (vetoNew.getSummedDet() < 3500) ? true : false;
  passedCutsArrayCnCWithTracking[8]  = (vetoNew.getSummedTightIso() < 800) ? true : false;
  passedCutsArrayCnCWithTracking[9]  = (vetoNew.getEcalBackEnergy() < 250) ? true : false;
  passedCutsArrayCnCWithTracking[10]  = (vetoNew.getNReadoutHits() < 70) ? true : false;
  passedCutsArrayCnCWithTracking[11]  = (vetoNew.getShowerRMS() < 110) ? true : false;
  passedCutsArrayCnCWithTracking[12]  = (vetoNew.getYStd() < 70) ? true : false;
  passedCutsArrayCnCWithTracking[13]  = (vetoNew.getMaxCellDep() < 300) ? true : false;
  passedCutsArrayCnCWithTracking[14]  = (vetoNew.getStdLayerHit() < 5) ? true : false;
  passedCutsArrayCnCWithTracking[15]  = (vetoNew.getNStraightTracks() < 3) ? true : false;
  passedCutsArrayCnCWithTracking[16]  = (hcalVeto.passesVeto()) ? true : false;
  passedCutsArrayCnCWithTracking[17]  = (vetoNew.getNStraightTracks() == 0) ? true : false;

  for (size_t i=0;i<sizeof(passedCutsArrayCnCWithTracking);i++) {
    bool allCutsPassedSoFar = true;
    for (size_t j=0;j<=i;j++) {
      if (!passedCutsArrayCnCWithTracking[j]) {
        allCutsPassedSoFar = false;
        break;
      }
    }
    if (allCutsPassedSoFar) {
      histograms_.fill("StdCutFlowWithTracking_RecoilX", i, vetoNew.getRecoilX() );
    }
  }

  // BDT based cutFlow with tracking
  bool passedCutsArrayTracking[12];
  std::fill(std::begin(passedCutsArrayTracking), std::end(passedCutsArrayTracking),false);
  passedCutsArrayTracking[0]  = (acceptance) ? true : false;
  passedCutsArrayTracking[1]  = (ignore_fiducial_analysis_ || (fiducial_analysis_ && vetoNew.getFiducial()) || (!fiducial_analysis_ && !vetoNew.getFiducial())) ? true : false;
  passedCutsArrayTracking[2]  = (trigResult.passed()) ? true : false;
  passedCutsArrayTracking[3]  = (ignore_tagger_analysis_ || (taggerP > 5600)) ? true : false;
  passedCutsArrayTracking[4]  = (recoilN == 1) ? true : false;
  passedCutsArrayTracking[5]  = (std::abs(recoilD0) < 10.) ? true : false;
  passedCutsArrayTracking[6]  = (std::abs(recoilZ0) < 40.) ? true : false;
  passedCutsArrayTracking[7]  = (vetoNew.getDisc() > 0.99741) ? true : false;
  passedCutsArrayTracking[8]  = (vetoNew.getNStraightTracks() < 3) ? true : false;
  passedCutsArrayTracking[9]  = (hcalVeto.passesVeto()) ? true : false;
  passedCutsArrayTracking[10]  = (vetoNew.getNStraightTracks() == 0) ? true : false;
  if (fiducial_analysis_) {
    passedCutsArrayTracking[11]  = (vetoNew.getEPAng() > 3)  ? true : false;
  } else {
    passedCutsArrayTracking[11]  = true;
  }

  for (size_t i=0;i<sizeof(passedCutsArrayTracking);i++) {
    bool allCutsPassedSoFar = true;
    for (size_t j=0;j<=i;j++) {
      if (!passedCutsArrayTracking[j]) {
        allCutsPassedSoFar = false;
        break;
      }
    }
    if (allCutsPassedSoFar) {
      histograms_.fill("TrackingCutFlow_RecoilX", i, vetoNew.getRecoilX() );
      if (i == (sizeof(passedCutsArrayTracking)-1) && !signal_) {
        std::cout << " This bkg event survived all the cuts!!!" << std::endl;
      }
      histograms_.fill("Tracking_TaggerP", i, taggerP);
      histograms_.fill("Tracking_RecoilN", i, recoilN);
      histograms_.fill("Tracking_RecoilP", i, recoilP);
      histograms_.fill("Tracking_RecoilPt", i, recoilPt);
      histograms_.fill("Tracking_RecoilD0", i, recoilD0);
      histograms_.fill("Tracking_RecoilZ0", i, recoilZ0);
    }
  }

  // BDT based cutFlow with tracking starting with Hcal and Ecal veto
  bool passedCutsArrayTrackingHcal[11];
  std::fill(std::begin(passedCutsArrayTrackingHcal), std::end(passedCutsArrayTrackingHcal),false);
  passedCutsArrayTrackingHcal[0]  = (acceptance) ? true : false;
  passedCutsArrayTrackingHcal[1]  = (ignore_fiducial_analysis_ || (fiducial_analysis_ && vetoNew.getFiducial()) || (!fiducial_analysis_ && !vetoNew.getFiducial())) ? true : false;
  passedCutsArrayTrackingHcal[2]  = (trigResult.passed()) ? true : false;
  passedCutsArrayTrackingHcal[3]  = (hcalVeto.passesVeto()) ? true : false;
  passedCutsArrayTrackingHcal[4]  = (vetoNew.getDisc() > 0.99741) ? true : false;
  passedCutsArrayTrackingHcal[5]  = (vetoNew.getNStraightTracks() == 0) ? true : false;
  if (fiducial_analysis_) {
    passedCutsArrayTrackingHcal[6]  = ((vetoNew.getEPAng() > 3))  ? true : false;
  } else {
    passedCutsArrayTrackingHcal[6]  = true;
  }
  passedCutsArrayTrackingHcal[7]  = (taggerP > 5600) ? true : false;
  passedCutsArrayTrackingHcal[8]  = (recoilN == 1) ? true : false;
  passedCutsArrayTrackingHcal[9]  = (std::abs(recoilD0) < 10) ? true : false;
  passedCutsArrayTrackingHcal[10]  = (std::abs(recoilZ0) < 40) ? true : false;


  for (size_t i=0;i<sizeof(passedCutsArrayTrackingHcal);i++) {
    bool allCutsPassedSoFar = true;
    for (size_t j=0;j<=i;j++) {
      if (!passedCutsArrayTrackingHcal[j]) {
        allCutsPassedSoFar = false;
        break;
      }
    }
    if (allCutsPassedSoFar) {
      histograms_.fill("TrackingCutFlowHcal_RecoilX", i, vetoNew.getRecoilX() );
      histograms_.fill("TrackingHcal_TaggerP", i, taggerP);
      histograms_.fill("TrackingHcal_RecoilN", i, recoilN);
      histograms_.fill("TrackingHcal_RecoilD0", i, recoilD0);
      histograms_.fill("TrackingHcal_RecoilZ0", i, recoilZ0);
    }
  }

  // --------------------------------------------------------------------------
  // Reverse cutflow, i.e. start with the last cut from the original cutflow
  bool passedCutsArrayReverse[13];
  std::fill(std::begin(passedCutsArrayReverse), std::end(passedCutsArrayReverse),false);
  passedCutsArrayReverse[0]  = (acceptance) ? true : false;
  passedCutsArrayReverse[1]  = (ignore_fiducial_analysis_ || (fiducial_analysis_ && vetoNew.getFiducial()) || (!fiducial_analysis_ && !vetoNew.getFiducial())) ? true : false;
  passedCutsArrayReverse[2]  = (trigResult.passed()) ? true : false;
  passedCutsArrayReverse[3]  = (hcalVeto.passesVeto()) ? true : false;
  passedCutsArrayReverse[4]  = (vetoNew.getNStraightTracks() < 3) ? true : false;
  passedCutsArrayReverse[5]  = (vetoNew.getStdLayerHit() < 5) ? true : false;
  passedCutsArrayReverse[6]  = (vetoNew.getMaxCellDep() < 300) ? true : false;
  passedCutsArrayReverse[7]  = (vetoNew.getYStd() < 70) ? true : false;
  passedCutsArrayReverse[8]  = (vetoNew.getShowerRMS() < 110) ? true : false;
  passedCutsArrayReverse[9]  = (vetoNew.getNReadoutHits() < 70) ? true : false;
  passedCutsArrayReverse[10]  = (vetoNew.getEcalBackEnergy() < 250) ? true : false;
  passedCutsArrayReverse[11]  = (vetoNew.getSummedTightIso() < 800) ? true : false;
  passedCutsArrayReverse[12]  = (vetoNew.getSummedDet() < 3500) ? true : false;

  for (size_t i=0;i<sizeof(passedCutsArrayReverse);i++) {
    bool allCutsPassedSoFar = true;
    for (size_t j=0;j<=i;j++) {
      if (!passedCutsArrayReverse[j]) {
        allCutsPassedSoFar = false;
      }
    }
    if (allCutsPassedSoFar) {
      histograms_.fill("Rev_AvgLayerHit", i, vetoNew.getAvgLayerHit() );
      histograms_.fill("Rev_DeepestLayerHit", i, vetoNew.getDeepestLayerHit() );
      histograms_.fill("Rev_EcalBackEnergy", i, vetoNew.getEcalBackEnergy() );
      histograms_.fill("Rev_EpAng", i, vetoNew.getEPAng() );
      histograms_.fill("Rev_EpSep", i, vetoNew.getEPSep() );
      histograms_.fill("Rev_FirstNearPhLayer", i, vetoNew.getFirstNearPhLayer() );
      histograms_.fill("Rev_MaxCellDep", i, vetoNew.getMaxCellDep() );
      histograms_.fill("Rev_NReadoutHits", i, vetoNew.getNReadoutHits() );
      histograms_.fill("Rev_StdLayerHit", i, vetoNew.getStdLayerHit() );
      histograms_.fill("Rev_Straight", i, vetoNew.getNStraightTracks() );
      // histograms_.fill("Rev_LinRegNew", i, vetoNew.getNLinRegTracks() );
      histograms_.fill("Rev_SummedDet", i, vetoNew.getSummedDet() );
      histograms_.fill("Rev_SummedTightIso", i, vetoNew.getSummedTightIso() );
      histograms_.fill("Rev_ShowerRMS", i, vetoNew.getShowerRMS() );
      histograms_.fill("Rev_XStd", i, vetoNew.getXStd() );
      histograms_.fill("Rev_YStd", i, vetoNew.getYStd() );
      histograms_.fill("Rev_Hcal_MaxPE", i, hcalMaxPE );
      histograms_.fill("Rev_Hcal_TotalPE", i, hcalTotalPe );
      histograms_.fill("Rev_Hcal_MaxTiming", i, hcalMaxTiming );
      histograms_.fill("Rev_Hcal_MaxSector", i, hcalMaxSector );
    }
  }
  
    // // N-1 plots
    // // << "      >> Doing N1 plots";
    // // i=0 is trigger, i=1 is fiducial
    //  for (size_t i=2;i<sizeof(passedCutsArrayCnC);i++) {
    //    bool allOtherCutsPassed = true;
    //    for (size_t j=2;j<sizeof(passedCutsArrayCnC);j++) {
    //      if (i==j) continue;
    //      if (!passedCutsArrayCnC[j]) {
    //        allOtherCutsPassed = false;
    //          // We found a cut that's not passed, no point in looking into the rest of them
    //        break;
    //      }
    //    }

    //    if (allOtherCutsPassed && trigResult.passed() && ((fiducial_analysis_ && vetoNew.getFiducial()) || (!fiducial_analysis_ && !vetoNew.getFiducial()))) {
    //     if (i==2) histograms_.fill("N1_SummedDet", i, vetoNew.getSummedDet() );
    //     if (i==3) histograms_.fill("N1_SummedTightIso", i, vetoNew.getSummedTightIso() );
    //     if (i==4) histograms_.fill("N1_EcalBackEnergy", i, vetoNew.getEcalBackEnergy() );
    //     if (i==5) histograms_.fill("N1_NReadoutHits", i, vetoNew.getNReadoutHits() );
    //     if (i==6) histograms_.fill("N1_ShowerRMS", i, vetoNew.getShowerRMS() );
    //     if (i==7) histograms_.fill("N1_YStd", i, vetoNew.getYStd() );
    //     if (i==8) histograms_.fill("N1_MaxCellDep", i, vetoNew.getMaxCellDep() );
    //     if (i==9) histograms_.fill("N1_StdLayerHit", i, vetoNew.getStdLayerHit() );
    //     if (i==10) histograms_.fill("N1_Straight", i, vetoNew.getNStraightTracks() );
    //     if (i==11) {
    //       histograms_.fill("N1_Hcal_MaxPE", i, hcalMaxPE );
    //       histograms_.fill("N1_Hcal_TotalPE", i, hcalTotalPe );
    //       histograms_.fill("N1_Hcal_MaxTiming", i, hcalMaxTiming );
    //       histograms_.fill("N1_Hcal_MaxSector", i, hcalMaxSector );
    //     }
    
    // histograms_.fill("N1_LinRegNew", i, vetoNew.getNLinRegTracks() );
    // histograms_.fill("N1_EpAng", i, vetoNew.getEPAng() );
    // histograms_.fill("N1_EpSep", i, vetoNew.getEPSep() );
    // histograms_.fill("N1_FirstNearPhLayer", i, vetoNew.getFirstNearPhLayer() );
    // histograms_.fill("N1_XStd", i, vetoNew.getXStd() );
    // histograms_.fill("N1_AvgLayerHit", i, vetoNew.getAvgLayerHit() );
    // histograms_.fill("N1_DeepestLayerHit", i, vetoNew.getDeepestLayerHit() );
  //  }
  // }
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

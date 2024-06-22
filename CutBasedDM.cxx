#include "Framework/EventProcessor.h"
#include "Ecal/Event/EcalVetoResult.h"
#include "Recon/Event/TriggerResult.h"
#include "Hcal/Event/HcalVetoResult.h"

// v0: Just a few cuts, establish minimal scenario
// v1: More strict cuts
// v2: Fix for N-1 plots, fix NReadoutHits binning
// v3: Max cell dep is not depth but deposition, non-fid option

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
  std::string trigger_collName_;
  std::string trigger_passName_;
  bool fiducial_;
};


void CutBasedDM::configure(framework::config::Parameters &ps) {
  trigger_collName_ = ps.getParameter<std::string>("trigger_name");
  trigger_passName_ = ps.getParameter<std::string>("trigger_pass");
  fiducial_ = ps.getParameter<bool>("fiducial");
  
  return;
}

void CutBasedDM::onProcessStart(){
  getHistoDirectory();
  
  histograms_.create("RecoilX", "", 20, -0.5, 19.5, "RecoilX [mm]", 90, -450.0, 450.0);
  histograms_.create("AvgLayerHit", "", 20, -0.5, 19.5, "AvgLayerHit", 35, -0.5, 34.5);
  histograms_.create("DeepestLayerHit", "", 20, -0.5, 19.5, "DeepestLayerHit", 35, -0.5, 34.5);
  histograms_.create("EcalBackEnergy", "", 20, -0.5, 19.5, "EcalBackEnergy", 100, 0.0, 3000.0);
  histograms_.create("EpAng", "", 20, -0.5, 19.5, "EpAng", 100, 0.0, 90.0);
  histograms_.create("EpSep", "", 20, -0.5, 19.5, "EpSep", 100, 0.0, 1000.0);
  histograms_.create("FirstNearPhLayerZ", "", 20, -0.5, 19.5, "FirstNearPhLayerZ", 100, 0.0, 800.0);
  histograms_.create("MaxCellDep", "", 20, -0.5, 19.5, "MaxCellDep", 100, 0.0, 800.0);
  histograms_.create("NReadoutHits", "", 20, -0.5, 19.5, "NReadoutHits", 150, -0.5, 149.5);
  histograms_.create("StdLayerHit", "", 20, -0.5, 19.5, "StdLayerHit", 70, -0.5, 34.5);
  histograms_.create("Straight", "", 20, -0.5, 19.5, "StraightTracks", 15, -0.5, 14.5);
  histograms_.create("LinRegNew", "", 20, -0.5, 19.5, "LinearRegressionTracks", 15, -0.5, 14.5);
  histograms_.create("SummedDet", "", 20, -0.5, 19.5, "SummedDet", 100, 0.0, 8000.0);
  histograms_.create("SummedTightIso", "", 20, -0.5, 19.5, "SummedTightIso", 100, 0.0, 8000.0);
  histograms_.create("ShowerRMS", "", 20, -0.5, 19.5, "ShowerRMS [mm]", 100, 0.0, 250.0);
  histograms_.create("XStd", "", 20, -0.5, 19.5, "XStd [mm]", 100, 0.0, 250.0);
  histograms_.create("YStd", "", 20, -0.5, 19.5, "YStd [mm]", 100, 0.0, 250.0);
  histograms_.create("RecoilPT", "", 20, -0.5, 19.5, "Recoil p_{T} [MeV]", 800, 0.0, 8000.0);
  histograms_.create("RecoilPZ", "", 20, -0.5, 19.5, "Recoil p_{Z} [MeV]", 800, -10.0, 8010.0);
  histograms_.create("RecoilP", "", 20, -0.5, 19.5, "Recoil p [MeV]", 800, 0.0, 8000.0);

  histograms_.create("Rev_AvgLayerHit", "", 20, -0.5, 19.5, "AvgLayerHit", 35, -0.5, 34.5);
  histograms_.create("Rev_DeepestLayerHit", "", 20, -0.5, 19.5, "DeepestLayerHit", 35, -0.5, 34.5);
  histograms_.create("Rev_EcalBackEnergy", "", 20, -0.5, 19.5, "EcalBackEnergy", 100, 0.0, 3000.0);
  histograms_.create("Rev_EpAng", "", 20, -0.5, 19.5, "EpAng", 100, 0.0, 90.0);
  histograms_.create("Rev_EpSep", "", 20, -0.5, 19.5, "EpSep", 100, 0.0, 1000.0);
  histograms_.create("Rev_FirstNearPhLayerZ", "", 20, -0.5, 19.5, "FirstNearPhLayerZ", 100, 0.0, 800.0);
  histograms_.create("Rev_MaxCellDep", "", 20, -0.5, 19.5, "MaxCellDep", 100, 0.0, 800.0);
  histograms_.create("Rev_NReadoutHits", "", 20, -0.5, 19.5, "NReadoutHits", 150, -0.5, 149.5);
  histograms_.create("Rev_StdLayerHit", "", 20, -0.5, 19.5, "StdLayerHit", 70, -0.5, 34.5);
  histograms_.create("Rev_Straight", "", 20, -0.5, 19.5, "StraightTracks", 15, -0.5, 14.5);
  histograms_.create("Rev_LinRegNew", "", 20, -0.5, 19.5, "LinearRegressionTracks", 15, -0.5, 14.5);
  histograms_.create("Rev_SummedDet", "", 20, -0.5, 19.5, "SummedDet", 100, 0.0, 8000.0);
  histograms_.create("Rev_SummedTightIso", "", 20, -0.5, 19.5, "SummedTightIso", 100, 0.0, 8000.0);
  histograms_.create("Rev_ShowerRMS", "", 20, -0.5, 19.5, "ShowerRMS [mm]", 100, 0.0, 250.0);
  histograms_.create("Rev_XStd", "", 20, -0.5, 19.5, "XStd [mm]", 100, 0.0, 250.0);
  histograms_.create("Rev_YStd", "", 20, -0.5, 19.5, "YStd [mm]", 100, 0.0, 250.0);

  // histograms_.create("N1_AvgLayerHit", "", 20, -0.5, 19.5, "AvgLayerHit", 35, -0.5, 34.5);
  // histograms_.create("N1_DeepestLayerHit", "", 20, -0.5, 19.5, "DeepestLayerHit", 35, -0.5, 34.5);
  histograms_.create("N1_EcalBackEnergy", "", 20, -0.5, 19.5, "EcalBackEnergy", 100, 0.0, 3000.0);
  // histograms_.create("N1_EpAng", "", 20, -0.5, 19.5, "EpAng", 100, 0.0, 90.0);
  // histograms_.create("N1_EpSep", "", 20, -0.5, 19.5, "EpSep", 100, 0.0, 1000.0);
  // histograms_.create("N1_FirstNearPhLayerZ", "", 20, -0.5, 19.5, "FirstNearPhLayerZ", 100, 0.0, 800.0);
  histograms_.create("N1_MaxCellDep", "", 20, -0.5, 19.5, "MaxCellDep", 100, 0.0, 800.0);
  histograms_.create("N1_NReadoutHits", "", 20, -0.5, 19.5, "NReadoutHits", 150, -0.5, 149.5);
  histograms_.create("N1_StdLayerHit", "", 20, -0.5, 19.5, "StdLayerHit", 70, -0.5, 34.5);
  histograms_.create("N1_Straight", "", 20, -0.5, 19.5, "StraightTracks", 15, -0.5, 14.5);
  // histograms_.create("N1_LinRegNew", "", 20, -0.5, 19.5, "LinearRegressionTracks", 15, -0.5, 14.5);
  histograms_.create("N1_SummedDet", "", 20, -0.5, 19.5, "SummedDet", 100, 0.0, 8000.0);
  histograms_.create("N1_SummedTightIso", "", 20, -0.5, 19.5, "SummedTightIso", 100, 0.0, 8000.0);
  histograms_.create("N1_ShowerRMS", "", 20, -0.5, 19.5, "ShowerRMS [mm]", 100, 0.0, 250.0);
  // histograms_.create("N1_XStd", "", 20, -0.5, 19.5, "XStd [mm]", 100, 0.0, 250.0);
  histograms_.create("N1_YStd", "", 20, -0.5, 19.5, "YStd [mm]", 100, 0.0, 250.0);

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
    "N_{straight} < 6",           // 11
    "PASS_{HCalVeto}",            // 12
    ""};

  std::vector<TH1 *> hists = {
  histograms_.get("AvgLayerHit"),
  histograms_.get("DeepestLayerHit"),
  histograms_.get("EcalBackEnergy"),
  histograms_.get("EpAng"),
  histograms_.get("EpSep"),
  histograms_.get("FirstNearPhLayerZ"),
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
  histograms_.get("RecoilX"),
  histograms_.get("RecoilPT"),
  histograms_.get("RecoilPZ"),

  // histograms_.get("N1_AvgLayerHit"),
  // histograms_.get("N1_DeepestLayerHit"),
  histograms_.get("N1_EcalBackEnergy"),
  // histograms_.get("N1_EpAng"),
  // histograms_.get("N1_EpSep"),
  // histograms_.get("N1_FirstNearPhLayerZ"),
  histograms_.get("N1_MaxCellDep"),
  histograms_.get("N1_NReadoutHits"),
  histograms_.get("N1_StdLayerHit"),
  histograms_.get("N1_Straight"),
  // histograms_.get("N1_LinRegNew"),
  histograms_.get("N1_SummedDet"),
  histograms_.get("N1_SummedTightIso"),
  histograms_.get("N1_ShowerRMS"),
  // histograms_.get("N1_XStd"),
  histograms_.get("N1_YStd"),
  };

  for (int ilabel{1}; ilabel < labels.size(); ++ilabel) {
    for (auto &hist : hists) {
      hist->GetXaxis()->SetBinLabel(ilabel, labels[ilabel - 1].c_str());
    }
  }

  std::vector<std::string> labels_Rev = {
    "All",
    "PASS_{HCalVeto}",           // 12
    "N_{straight} < 6",           // 11
    "RMS_{Layer,hit} < 5",        // 10
    "E_{cell,max} < 300",         // 9
    "RMS_{shower,Y} < 70",        // 8
    "RMS_{shower} < 110",         // 7
    "N_{hits} < 70",              // 6
    "E_{back} < 250",             // 5
    "E_{SumTight} < 800",        // 4
    "E_{sum} < 3500",       // 3
    "Fiducial",              // 2
    "Triggerred",            // 1
    ""};

  std::vector<TH1 *> hists_Rev = {
    histograms_.get("Rev_AvgLayerHit"),
    histograms_.get("Rev_DeepestLayerHit"),
    histograms_.get("Rev_EcalBackEnergy"),
    histograms_.get("Rev_EpAng"),
    histograms_.get("Rev_EpSep"),
    histograms_.get("Rev_FirstNearPhLayerZ"),
    histograms_.get("Rev_MaxCellDep"),
    histograms_.get("Rev_NReadoutHits"),
    histograms_.get("Rev_StdLayerHit"),
    histograms_.get("Rev_Straight"),
    histograms_.get("Rev_LinRegNew"),
    histograms_.get("Rev_SummedDet"),
    histograms_.get("Rev_SummedTightIso"),
    histograms_.get("Rev_ShowerRMS"),
    histograms_.get("Rev_XStd"),
    histograms_.get("Rev_YStd")
  };

    for (int ilabel{1}; ilabel < labels_Rev.size(); ++ilabel) {
    for (auto &hist : hists_Rev) {
      hist->GetXaxis()->SetBinLabel(ilabel, labels_Rev[ilabel - 1].c_str());
    }
  }

}

void CutBasedDM::analyze(const framework::Event& event) {
  // std::cout << " ---------------------------------------------" << std::endl;
  auto vetoNew{event.getObject<ldmx::EcalVetoResult>("EcalVetoNew","")};
  auto trigResult{event.getObject<ldmx::TriggerResult>(trigger_collName_, trigger_passName_)};
  auto hcalVeto{event.getObject<ldmx::HcalVetoResult>("HcalVeto","")};

  float pT2{-9999.};
  float pT{-9999.};
  float pZ{-9999.};
  float totMom2{-9999.};
  float totMom{-9999.};

  pT2 = vetoNew.getRecoilMomentum()[0]*vetoNew.getRecoilMomentum()[0] + vetoNew.getRecoilMomentum()[1]*vetoNew.getRecoilMomentum()[1];
  if (pT2 >= 0) pT = sqrt(pT2);
  pZ =  vetoNew.getRecoilMomentum()[2];
  if (pZ < 0 && pZ != -9999) pZ = -5;
  if (pT2 >= 0) totMom2 = pT2 + pZ*pZ;
  if (totMom2 >= 0) totMom = sqrt(totMom2);
  

  // CutFlow here
  bool passedCutsArray[12];
  std::fill(std::begin(passedCutsArray), std::end(passedCutsArray),false);
  passedCutsArray[0]  = (trigResult.passed()) ? true : false;

  passedCutsArray[1]  = ((fiducial_ && vetoNew.getRecoilX() !=-9999) || (!fiducial_ && vetoNew.getRecoilX() ==-9999)) ? true : false;
  passedCutsArray[2]  = (vetoNew.getSummedDet() < 3500) ? true : false;
  passedCutsArray[3]  = (vetoNew.getSummedTightIso() < 800) ? true : false;
  passedCutsArray[4]  = (vetoNew.getEcalBackEnergy() < 250) ? true : false;
  passedCutsArray[5]  = (vetoNew.getNReadoutHits() < 70) ? true : false;
  passedCutsArray[6]  = (vetoNew.getShowerRMS() < 110) ? true : false;
  passedCutsArray[7]  = (vetoNew.getYStd() < 70) ? true : false;
  passedCutsArray[8]  = (vetoNew.getMaxCellDep() < 300) ? true : false;
  passedCutsArray[9]  = (vetoNew.getStdLayerHit() < 5) ? true : false;
  passedCutsArray[10]  = (vetoNew.getNStraightTracks() < 6) ? true : false;
  passedCutsArray[11]  = (hcalVeto.passesVeto()) ? true : false;

  histograms_.fill("RecoilX", 0. , vetoNew.getRecoilX() );
  histograms_.fill("AvgLayerHit", 0. , vetoNew.getAvgLayerHit() );
  histograms_.fill("DeepestLayerHit", 0. , vetoNew.getDeepestLayerHit() );
  histograms_.fill("EcalBackEnergy", 0. , vetoNew.getEcalBackEnergy() );
  histograms_.fill("EpAng", 0. , vetoNew.getEPAng() );
  histograms_.fill("EpSep", 0. , vetoNew.getEPSep() );
  histograms_.fill("FirstNearPhLayerZ", 0. , vetoNew.getFirstNearPhLayer() );
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
  histograms_.fill("RecoilPT", 0. ,pT );
  histograms_.fill("RecoilPZ", 0. , pZ );
  histograms_.fill("RecoilP", 0. , totMom );
  
  for (size_t i=0;i<sizeof(passedCutsArray);i++) {
    bool allCutsPassedSoFar = true;
    for (size_t j=0;j<=i;j++) {
      if (!passedCutsArray[j]) {
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
      // << std::endl;

      histograms_.fill("RecoilX", i+1 , vetoNew.getRecoilX() );
      histograms_.fill("AvgLayerHit", i+1 , vetoNew.getAvgLayerHit() );
      histograms_.fill("DeepestLayerHit", i+1 , vetoNew.getDeepestLayerHit() );
      histograms_.fill("EcalBackEnergy", i+1 , vetoNew.getEcalBackEnergy() );
      histograms_.fill("EpAng", i+1 , vetoNew.getEPAng() );
      histograms_.fill("EpSep", i+1 , vetoNew.getEPSep() );
      histograms_.fill("FirstNearPhLayerZ", i+1 , vetoNew.getFirstNearPhLayer() );
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
      histograms_.fill("RecoilPT", i+1 ,pT );
      histograms_.fill("RecoilPZ", i+1 , pZ );
      histograms_.fill("RecoilP", i+1 , totMom );
    }
  }

  
   // Reverse cutflow, i.e. start with the last cut from the original cutflow
histograms_.fill("Rev_AvgLayerHit", 0. , vetoNew.getAvgLayerHit() );
histograms_.fill("Rev_DeepestLayerHit", 0. , vetoNew.getDeepestLayerHit() );
histograms_.fill("Rev_EcalBackEnergy", 0. , vetoNew.getEcalBackEnergy() );
histograms_.fill("Rev_EpAng", 0. , vetoNew.getEPAng() );
histograms_.fill("Rev_EpSep", 0. , vetoNew.getEPSep() );
histograms_.fill("Rev_FirstNearPhLayerZ", 0. , vetoNew.getFirstNearPhLayer() );
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

bool passedCutsArrayReverse[12];
std::reverse_copy(std::begin(passedCutsArray), std::end(passedCutsArray), std::begin(passedCutsArrayReverse));
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
    histograms_.fill("Rev_FirstNearPhLayerZ", i+1 , vetoNew.getFirstNearPhLayer() );
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
  }
}
  

  
   // N-1 plots
  // << "      >> Doing N1 plots";
  // i=0 is trigger, i=1 is fiducial
 for (size_t i=2;i<sizeof(passedCutsArray);i++) {
   bool allOtherCutsPassed = true;
   for (size_t j=2;j<sizeof(passedCutsArray);j++) {
     if (i==j) continue;
     if (!passedCutsArray[j]) {
       allOtherCutsPassed = false;
         // We found a cut that's not passed, no point in looking into the rest of them
       break;
     }
   }

   if (allOtherCutsPassed && trigResult.passed() && ((fiducial_ && vetoNew.getRecoilX() !=-9999) || (!fiducial_ && vetoNew.getRecoilX() ==-9999))) {
    if (i==2) histograms_.fill("N1_SummedDet", i+1 , vetoNew.getSummedDet() );
    if (i==3) histograms_.fill("N1_SummedTightIso", i+1 , vetoNew.getSummedTightIso() );
    if (i==4) histograms_.fill("N1_EcalBackEnergy", i+1 , vetoNew.getEcalBackEnergy() );
    if (i==5) histograms_.fill("N1_NReadoutHits", i+1 , vetoNew.getNReadoutHits() );
    if (i==6) histograms_.fill("N1_ShowerRMS", i+1 , vetoNew.getShowerRMS() );
    if (i==7) histograms_.fill("N1_YStd", i+1 , vetoNew.getYStd() );
    if (i==8) histograms_.fill("N1_MaxCellDep", i+1 , vetoNew.getMaxCellDep() );
    if (i==9) histograms_.fill("N1_StdLayerHit", i+1 , vetoNew.getStdLayerHit() );
    if (i==10) histograms_.fill("N1_Straight", i+1 , vetoNew.getNStraightTracks() );
    // histograms_.fill("N1_LinRegNew", i+1 , vetoNew.getNLinRegTracks() );
    // histograms_.fill("N1_EpAng", i+1 , vetoNew.getEPAng() );
    // histograms_.fill("N1_EpSep", i+1 , vetoNew.getEPSep() );
    // histograms_.fill("N1_FirstNearPhLayerZ", i+1 , vetoNew.getFirstNearPhLayer() );
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

DECLARE_ANALYZER(CutBasedDM);

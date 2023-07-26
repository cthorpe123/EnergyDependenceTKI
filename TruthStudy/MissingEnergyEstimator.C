R__LOAD_LIBRARY(STV_Tools_cxx.so);

#include "STV_Tools.h"
#include "Funcs.h"

// Calculate the missing hadronic energy

void MissingEnergyEstimator(){

  double data_pot = 6e20;
  bool save_weights = false;

  // load the MC
  TFile* f_in = TFile::Open("/pnfs/uboone/persistent/users/davidc/searchingfornues/v08_00_00_43/0928/prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root"); 
  const double mc_pot = 1.31E+21;

  double pot_weight = data_pot/mc_pot;

  // Load the branches needed
  TTree* t_in = (TTree*)f_in->Get("nuselection/NeutrinoSelectionFilter");
  float nu_e; 
  vector<int>* mc_pdg=0;
  vector<float>* mc_E=0;
  vector<float>* mc_px=0;
  vector<float>* mc_py=0;
  vector<float>* mc_pz=0;
  map<string,vector<double>>* weights=0;
  vector<unsigned short>* weightsFlux=0;
  vector<unsigned short>* weightsGenie=0;
  vector<unsigned short>* weightsReint=0;
  float weightSplineTimesTune;

  t_in->SetBranchAddress("nu_e",&nu_e);
  t_in->SetBranchAddress("mc_pdg",&mc_pdg);
  t_in->SetBranchAddress("mc_E",&mc_E);
  t_in->SetBranchAddress("mc_px",&mc_px);
  t_in->SetBranchAddress("mc_py",&mc_py);
  t_in->SetBranchAddress("mc_pz",&mc_pz);
  t_in->SetBranchAddress("weightSplineTimesTune",&weightSplineTimesTune);
  if(save_weights){
    //t_in->SetBranchAddress("weights",&weights);
    t_in->SetBranchAddress("weightsFlux",&weightsFlux);
    t_in->SetBranchAddress("weightsGenie",&weightsGenie);
    //t_in->SetBranchAddress("weightsReint",&weightsReint);
  }

  std::vector<float>* mc_p = new vector<float>;

  // store everything important in a TTree for future analysis   
  TTree* t_out = new TTree("t_out","TKI Tree");

  Float_t missing_e;
  Float_t total_e;
  Float_t muon_e;
  Float_t muon_costheta;
  Float_t had_e;
  Float_t deltapt;
  Float_t deltaalphat;
  Float_t deltaphit; 

  t_out->Branch("nu_e",&nu_e);
  t_out->Branch("missing_e",&missing_e);
  t_out->Branch("total_e",&total_e);
  t_out->Branch("muon_e",&muon_e);
  t_out->Branch("muon_costheta",&muon_costheta);
  t_out->Branch("had_e",&had_e);
  t_out->Branch("deltapt",&deltapt);
  t_out->Branch("deltaalphat",&deltaalphat);
  t_out->Branch("deltaphit",&deltaphit);
  t_out->Branch("weightSplineTimesTune",&weightSplineTimesTune);
  if(save_weights){
    t_out->Branch("weights",&weights);
    t_out->Branch("weightsFlux",&weightsFlux);
    t_out->Branch("weightsGenie",&weightsGenie);
    //t_out->Branch("weightsReint",&weightsReint);
  }

  t_out->SetDirectory(0);

  int nentries = t_in->GetEntries();
  for(int ientry=0;ientry<nentries;ientry++){
    if(ientry % 50000 == 0) std::cout << "Entry " << ientry << "/" << nentries << std::endl;
    t_in->GetEntry(ientry);
    mc_p->clear();
    had_e = 0.0;
    if(std::isinf(weightSplineTimesTune) || std::isnan(weightSplineTimesTune)) continue;

    for(size_t i_mc=0;i_mc<mc_pdg->size();i_mc++) mc_p->push_back(TVector3(mc_px->at(i_mc),mc_py->at(i_mc),mc_pz->at(i_mc)).Mag());
    bool IsSignal = Is1muNp(mc_pdg,mc_p);
    if(!IsSignal) continue; 

    missing_e = MissingEnergy(mc_pdg,mc_p,mc_E);
    total_e = TotalFSEnergy(mc_pdg,mc_p,mc_E);

    // Get the three variables previously used in the constraint
    TVector3 Muon3Momentum(0,0,0),Hadron3Momentum(0,0,0);
    double MuonEnergy=0,HadronEnergy=0;
    for(size_t i_mc=0;i_mc<mc_pdg->size();i_mc++){
      if(abs(mc_pdg->at(i_mc)) == 13){
        muon_e = mc_E->at(i_mc);
        muon_costheta = cos(TVector3(mc_px->at(i_mc),mc_py->at(i_mc),mc_pz->at(i_mc)).Theta());
        Muon3Momentum = TVector3(mc_px->at(i_mc),mc_py->at(i_mc),mc_pz->at(i_mc));
        MuonEnergy = mc_E->at(i_mc);
      } 
      if(mc_pdg->at(i_mc) == 2212 && mc_p->at(i_mc) > 0.3){
        had_e += mc_E->at(i_mc) - mass_proton + 0.0086;
        Hadron3Momentum += TVector3(mc_px->at(i_mc),mc_py->at(i_mc),mc_pz->at(i_mc));
        HadronEnergy += mc_E->at(i_mc);
      }
      if(abs(mc_pdg->at(i_mc)) == 211 && mc_p->at(i_mc) > 0.1){
        had_e += mc_E->at(i_mc);
        Hadron3Momentum += TVector3(mc_px->at(i_mc),mc_py->at(i_mc),mc_pz->at(i_mc));
        HadronEnergy += mc_E->at(i_mc);
      }
    }

    // Calculate TKI Variables

    STV_Tools stv_tools(Muon3Momentum,Hadron3Momentum,MuonEnergy,HadronEnergy);
    deltapt = TVector3(Muon3Momentum.X() - Hadron3Momentum.X(),Muon3Momentum.Y() - Hadron3Momentum.Y(),0).Mag();
    deltaalphat = stv_tools.ReturnDeltaAlphaT();
    deltaphit = stv_tools.ReturnDeltaPhiT();

    //std::cout << std::endl;
    //std::cout << "Muon:  " <<  Muon3Momentum.X() << " " << Muon3Momentum.Y() << " " << Muon3Momentum.Z() << std::endl;
    //std::cout << "Hadron:  " <<  Hadron3Momentum.X() << " " << Hadron3Momentum.Y() << " " << Hadron3Momentum.Z() << std::endl;
    //std::cout << "deltapt = " << deltapt << "  deltaalphat = " << deltaalphat << "  deltaphit = " << deltaphit << std::endl;   

    t_out->Fill();

  }

  f_in->Close();

  std::cout << "Writing events to output" << std::endl;
  gSystem->Exec("mkdir -p rootfiles/");
  string name = save_weights ? "rootfiles/out_3.root" : "rootfiles/out_noweights.root";
  TFile* f_out = TFile::Open(name.c_str(),"RECREATE");
  t_out->Write();
  f_out->Close();
}

#include "Funcs.h"

void DrawDists(){

  double data_pot = 6e20;
  bool save_weights = true;

  string name = save_weights ? "rootfiles/out.root" : "rootfiles/out_noweights.root";
  TFile* f_in = TFile::Open(name.c_str());
  const double mc_pot = 1.31E+21;

  // store everything important in a TTree for future analysis   
  TTree* t_in = (TTree*)f_in->Get("t_out");

  Float_t nu_e; 
  Float_t missing_e;
  Float_t total_e;
  Float_t muon_e;
  Float_t muon_costheta;
  Float_t had_e;
  Float_t deltapt;
  Float_t deltaalphat;
  Float_t deltaphit; 
  map<string,vector<double>>* weights=0;
  vector<unsigned short>* weightsFlux=0;
  vector<unsigned short>* weightsGenie=0;
  vector<unsigned short>* weightsReint=0;
  float weightSplineTimesTune;

  t_in->SetBranchAddress("nu_e",&nu_e);
  t_in->SetBranchAddress("missing_e",&missing_e);
  t_in->SetBranchAddress("total_e",&total_e);
  t_in->SetBranchAddress("muon_e",&muon_e);
  t_in->SetBranchAddress("muon_costheta",&muon_costheta);
  t_in->SetBranchAddress("had_e",&had_e);
  t_in->SetBranchAddress("deltapt",&deltapt);
  t_in->SetBranchAddress("deltaalphat",&deltaalphat);
  t_in->SetBranchAddress("deltaphit",&deltaphit);
  t_in->SetBranchAddress("weightSplineTimesTune",&weightSplineTimesTune);
  if(save_weights){
    //t_in->SetBranchAddress("weights",&weights);
    t_in->SetBranchAddress("weightsFlux",&weightsFlux);
    t_in->SetBranchAddress("weightsGenie",&weightsGenie);
    //t_in->SetBranchAddress("weightsReint",&weightsReint);
  }

  TH1D* h_muon_e = new TH1D("h_muon_e",";Muon Energy (GeV);Events",40,0.0,2.0);
  TH1D* h_muon_costheta = new TH1D("h_muon_costheta",";Muon Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_had_e = new TH1D("h_had_e",";Hadronic Energy (GeV);Events",40,0.0,1.1);
  TH1D* h_deltapt = new TH1D("h_deltapt",";#delta p_{T} (GeV);Events",40,0.0,2.0);
  TH1D* h_deltaalphat = new TH1D("h_deltaalphat",";#delta #alpha_{T};Events",40,0.0,180.0);
  TH1D* h_deltaphit = new TH1D("h_deltaphit",";#delta #phi_{T};Events",40,0.0,180);
  TH1D* h_missing_e = new TH1D("h_missing_e",";Missing Energy (GeV);Events",40,0.0,0.3);
  TH1D* h_missing_e_frac = new TH1D("h_missing_e_frac",";Missing Energy/E_{#nu};Events",40,0.0,0.5);

  t_in->GetEntry(0);

  std::vector<TH1D*> h_muon_e_v;
  std::vector<TH1D*> h_muon_costheta_v;
  std::vector<TH1D*> h_had_e_v;
  std::vector<TH1D*> h_deltapt_v;
  std::vector<TH1D*> h_deltaalphat_v;
  std::vector<TH1D*> h_deltaphit_v;
  std::vector<TH1D*> h_missing_e_v;
  std::vector<TH1D*> h_missing_e_frac_v;
  if(save_weights){
    for(size_t i_w=0;i_w<weightsGenie->size();i_w++){
      h_muon_e_v.push_back((TH1D*)h_muon_e->Clone(("h_muon_e_v_"+std::to_string(i_w)).c_str()));
      h_muon_costheta_v.push_back((TH1D*)h_muon_costheta->Clone(("h_muon_costheta_v_"+std::to_string(i_w)).c_str()));
      h_had_e_v.push_back((TH1D*)h_had_e->Clone(("h_had_e_v_"+std::to_string(i_w)).c_str()));
      h_deltapt_v.push_back((TH1D*)h_deltapt->Clone(("h_deltapt_v_"+std::to_string(i_w)).c_str()));
      h_deltaalphat_v.push_back((TH1D*)h_deltaalphat->Clone(("h_deltaalphat_v_"+std::to_string(i_w)).c_str()));
      h_deltaphit_v.push_back((TH1D*)h_deltaphit->Clone(("h_deltaphit_v_"+std::to_string(i_w)).c_str()));
      h_missing_e_v.push_back((TH1D*)h_missing_e->Clone(("h_missing_e_v_"+std::to_string(i_w)).c_str()));
      h_missing_e_frac_v.push_back((TH1D*)h_missing_e_frac->Clone(("h_missing_e_frac_v_"+std::to_string(i_w)).c_str()));
    }
  }

  double pot_weight = data_pot/mc_pot;

  for(int ientry=0;ientry<t_in->GetEntries();ientry++){
    if(ientry % 10000 == 0) std::cout << "Entry " << ientry << "/" << t_in->GetEntries() << std::endl;
    t_in->GetEntry(ientry);
    double weight = pot_weight*weightSplineTimesTune;
    h_muon_e->Fill(muon_e,weight);
    h_muon_costheta->Fill(muon_costheta,weight);
    h_had_e->Fill(had_e,weight);
    h_deltapt->Fill(deltapt,weight);
    h_deltaalphat->Fill(deltaalphat,weight);
    h_deltaphit->Fill(deltaphit,weight);
    h_missing_e->Fill(missing_e,weight);
    h_missing_e_frac->Fill(missing_e/nu_e,weight);

    if(save_weights){
      for(size_t i_w=0;i_w<weightsGenie->size();i_w++){
        //std::cout << weight*weightsGenie->at(i_w)/1000 << std::endl;
        double weight2 = pot_weight*weightsGenie->at(i_w)/1000;
        h_muon_e_v.at(i_w)->Fill(muon_e,weight2);
        h_muon_costheta_v.at(i_w)->Fill(muon_costheta,weight2);
        h_had_e_v.at(i_w)->Fill(had_e,weight2);
        h_deltapt_v.at(i_w)->Fill(deltapt,weight2);
        h_deltaalphat_v.at(i_w)->Fill(deltaalphat,weight2);
        h_deltaphit_v.at(i_w)->Fill(deltaphit,weight2);
        h_missing_e_v.at(i_w)->Fill(missing_e,weight2);
        h_missing_e_frac_v.at(i_w)->Fill(missing_e/nu_e,weight2);
      }
    }

  }  

  DrawVariations("muon_e",h_muon_e,h_muon_e_v);
  DrawVariations("muon_costheta",h_muon_costheta,h_muon_costheta_v);
  DrawVariations("had_e",h_had_e,h_had_e_v);
  DrawVariations("deltapt",h_deltapt,h_deltapt_v);
  DrawVariations("deltaalphat",h_deltaalphat,h_deltaalphat_v);
  DrawVariations("deltaphit",h_deltaphit,h_deltaphit_v);
  DrawVariations("missing_e",h_missing_e,h_missing_e_v);
  DrawVariations("missing_e_frac",h_missing_e_frac,h_missing_e_frac_v);

  DrawCovarianceMatrix("missing_e_frac_muon_e",{h_missing_e_frac,h_muon_e},{h_missing_e_frac_v,h_muon_e_v});
  DrawCovarianceMatrix("missing_e_frac_muon_costheta",{h_missing_e_frac,h_muon_costheta},{h_missing_e_frac_v,h_muon_costheta_v});
  DrawCovarianceMatrix("missing_e_frac_had_e",{h_missing_e_frac,h_had_e},{h_missing_e_frac_v,h_had_e_v});
  DrawCovarianceMatrix("missing_e_frac_deltapt",{h_missing_e_frac,h_deltapt},{h_missing_e_frac_v,h_deltapt_v});
  DrawCovarianceMatrix("missing_e_frac_deltaalphat",{h_missing_e_frac,h_deltaalphat},{h_missing_e_frac_v,h_deltaalphat_v});
  DrawCovarianceMatrix("missing_e_frac_deltaphit",{h_missing_e_frac,h_deltaphit},{h_missing_e_frac_v,h_deltaphit_v});

}

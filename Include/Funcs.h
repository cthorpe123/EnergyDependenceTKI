#ifndef _Funcs_h_
#define _Funcs_h_

const double mass_proton = 0.9383;
const double mass_neutron = 0.9396;
const double mass_pion = 0.1396;
const double mass_pizero = 0.1350;


// Function to identify which events are true 1muNp

bool Is1muNp(const std::vector<int>* mc_pdg,const std::vector<float>* mc_p){

  bool has_muon=false,has_proton=false,has_pion=false;
  for(size_t i_mc=0;i_mc<mc_pdg->size();i_mc++){
    if(abs(mc_pdg->at(i_mc)) == 13 && mc_p->at(i_mc) > 0.1) has_muon = true;
    if(abs(mc_pdg->at(i_mc)) == 2212 && mc_p->at(i_mc) > 0.3) has_proton = true;
    if((abs(mc_pdg->at(i_mc)) == 211 && mc_p->at(i_mc) > 0.1) || mc_pdg->at(i_mc) == 111) has_pion = true;// TODO: what's our detection threshold for pi0?
  }

  return has_muon && has_proton && !has_pion;
}

// Total final state energy (proxy for neutrino energy)

double TotalFSEnergy(const std::vector<int>* mc_pdg,const std::vector<float>* mc_p,const std::vector<float>* mc_E){

  double energy = 0;

  for(size_t i_mc=0;i_mc<mc_pdg->size();i_mc++){
    if(mc_pdg->at(i_mc) == 2212 || mc_pdg->at(i_mc) == 2112) energy += mc_E->at(i_mc) - mass_proton + 0.0086; // nucleons add the KE
    if(abs(mc_pdg->at(i_mc)) == 13 || abs(mc_pdg->at(i_mc)) == 211 || mc_pdg->at(i_mc) == 111) energy += mc_E->at(i_mc); // muons and pions add the total energy
  }

  return energy;
}

// Truth based estimator of missing energy. Add together the energies of the 
// invisible hadrons (neutrons and hadrons below threshold)

double MissingEnergy(const std::vector<int>* mc_pdg,const std::vector<float>* mc_p,const std::vector<float>* mc_E){

  double energy = 0;

  for(size_t i_mc=0;i_mc<mc_pdg->size();i_mc++){
    if(mc_pdg->at(i_mc) == 2112) energy += mc_E->at(i_mc) - mass_neutron + 0.0086; // nucleons add the KE
    if(abs(mc_pdg->at(i_mc)) == 211 && mc_p->at(i_mc) < 0.1) energy += mc_E->at(i_mc);
    if(mc_pdg->at(i_mc) == 111) energy += mc_E->at(i_mc); // muons and pions add the total energy
  }

  return energy;

}

// Draw some plots with systematic variations

void DrawVariations(string title,TH1D* h_CV,vector<TH1D*> h_Alt={}){

  TCanvas* c = new TCanvas("c","c");  
  string xaxistitle = h_CV->GetXaxis()->GetTitle();
  string yaxistitle = h_CV->GetYaxis()->GetTitle();
  string plottitle = ";" + xaxistitle + ";" + yaxistitle;
  THStack* hs = new THStack("hs",plottitle.c_str());

  for(TH1D* h : h_Alt){
    h->SetLineColor(3);
    hs->Add(h);
  }
  h_CV->SetLineColor(1);
  h_CV->SetLineWidth(2);
  hs->Add(h_CV);

  hs->Draw("nostack HIST");

  gSystem->Exec("mkdir -p  Plots/");
  c->Print(("Plots/" + title + ".png").c_str());
  c->Close();
}

double CalcFracCovarianceMultisim(double CV1,double CV2,std::vector<double> Alt1,std::vector<double> Alt2){

  if(Alt1.size() != Alt2.size()) throw std::invalid_argument("CalcFracCovarianceMultisim: Alt1.size() != Alt2.size()");

  double x = 0.0;
  for(size_t i=0;i<Alt1.size();i++) x += (CV1 - Alt1.at(i))*(CV2 - Alt2.at(i));
  x /= Alt1.size()*CV1*CV2;

  return x;
}

// Calculate the covariance matrix for multiple variables

void DrawCovarianceMatrix(string title,std::vector<TH1D*> h_CV,std::vector<std::vector<TH1D*>> h_Alt){

  // check sizes are consistent first
  if(h_CV.size() != h_Alt.size()) throw std::invalid_argument("DrawCovarianceMatrix: h_CV.size() != h_Alt.size()");

  std::vector<std::vector<TH2D*>> Matrices(h_CV.size(),std::vector<TH2D*>(h_CV.size(),nullptr)); 

  for(size_t i=0;i<Matrices.size();i++){
    for(size_t j=0;j<Matrices.size();j++){
      Matrices.at(i).at(j) = new TH2D(("h_Matrix_" + title + "_" + std::to_string(10*i+j)).c_str(),"",h_CV.at(i)->GetNbinsX(),0.5,h_CV.at(i)->GetNbinsX()+0.5,h_CV.at(j)->GetNbinsX(),0.5,h_CV.at(j)->GetNbinsX()+0.5);
      for(int i_b=0;i_b<h_CV.at(i)->GetNbinsX()+1;i_b++){
        for(int j_b=0;j_b<h_CV.at(j)->GetNbinsX()+1;j_b++){
          double cv1 = h_CV.at(i)->GetBinContent(i_b);
          double cv2 = h_CV.at(j)->GetBinContent(j_b);
          std::vector<double> alt1,alt2;
          for(size_t i_w=0;i_w<h_Alt.at(i).size();i_w++) alt1.push_back(h_Alt.at(i).at(i_w)->GetBinContent(i_b));        
          for(size_t j_w=0;j_w<h_Alt.at(j).size();j_w++) alt2.push_back(h_Alt.at(j).at(j_w)->GetBinContent(j_b));        
          Matrices.at(i).at(j)->SetBinContent(i_b,j_b,CalcFracCovarianceMultisim(cv1,cv2,alt1,alt2));
        }
      }
      Matrices.at(i).at(j)->GetXaxis()->SetTitle(h_CV.at(i)->GetXaxis()->GetTitle());
      Matrices.at(i).at(j)->GetYaxis()->SetTitle(h_CV.at(j)->GetXaxis()->GetTitle());
    }  
  }

  // Draw the plots
  TCanvas* c = new TCanvas("c","c");
  gSystem->Exec("mkdir -p Plots/");
  for(size_t i=0;i<Matrices.size();i++){
    for(size_t j=0;j<Matrices.size();j++){
      Matrices.at(i).at(j)->Draw("colz");
      Matrices.at(i).at(j)->SetStats(0);
      string plotname = "Plots/FCov_" + title + "_" + std::to_string(i) + std::to_string(j) + ".png";
      c->Print(plotname.c_str());
    }
  }

}

#endif

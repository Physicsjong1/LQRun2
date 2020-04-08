#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "Math/Vector4Dfwd.h"
#include "TStyle.h"
#include <vector>

using namespace ROOT::VecOps;
using rvec_f = RVec<float>;
using rvec_b = RVec<bool>;
using rvec_i = RVec<int>;

template <typename T>
void plot(T hist, TString name){
  TCanvas * c = new TCanvas("c",Form("c_%s", name.Data()));
  hist->Write();
  hist->DrawClone();
  //c->Print(Form("hist_sig/%s.pdf",name.Data()));
  c->Print(Form("hist_bkg/%s.pdf",name.Data()));
  //c->Print(Form("%s.pdf",name.Data()));
}

vector<int> goodmuons_idx(rvec_i g){
  vector<int> out;
  for(int i = 0; i < g.size(); i++){
    if( g[i] ) out.push_back( i );
  }
  return out; 
}

void ana_lq(){

  //ROOT::RDataFrame df("Events", "/xrootd/store/data/Run2018B/SingleMuon/NANOAOD/Nano14Dec2018-v1/10000/BCC1B466-EF27-9D40-A0B7-9FE64F456E13.root");
  //ROOT::RDataFrame df("Events", "signal_nano.root");
  //ROOT::RDataFrame df("Events", "/xrootd/store/user/ljw1015/LQ_Signals/LQ_2017_nano_v1.root");
  //ROOT::RDataFrame df("Events","/cms/ldap_home/ljw1015/public/LQ_Signals/LQ_2016_nano.root");
  ROOT::RDataFrame df("Events","/xrootd/store/mc/RunIISummer16NanoAODv6/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/NANOAODSIM/PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v2/30000/A7753E9B-A80C-9547-8729-17F953C5B5DE.root");
  //ROOT::RDataFrame df("Events", "/xrootd/store/user/ljw1015/LQ_Signals/LQ_2016_nano.root");

  // Run Selection
  int year = 2016;
  //int year = 2017;
  //int year = 2018;
  
  //good muon selection

  auto df_S1_muon = df.Filter("nMuon >= 1", "Events with one lepton")
                       .Define("goodmuons","Muon_pt > 20  && abs(Muon_eta) < 2.4 && Muon_tightId && Muon_pfRelIso04_all < 0.15")
                       .Define("goodtaus","Tau_pt > 30 && abs(Tau_eta) < 2.4 && (Tau_idMVAoldDM & 2)")
                       .Define("goodjets","Jet_pt > 20 && abs(Jet_eta) < 2.4 && Jet_jetId >= 1");
  auto d2 = df_S1_muon.Display({"nTau","Tau_idMVAoldDM","goodtaus"},100);
  d2->Print();
  
  auto df_S1_goodmuon = df_S1_muon.Filter("Sum( goodmuons ) >=1 ","Events with at least a goood muon")
                                  .Define("goodmuons_idx",goodmuons_idx,{"goodmuons"});
  
  auto df_S1_goodtau = df_S1_goodmuon;
  if(year == 2016){
    df_S1_goodtau = df_S1_goodmuon.Filter("Sum( goodtaus ) >=1 ","Events with at least a goood tau")
                                  .Define("CvsB","Jet_btagDeepFlavC/(Jet_btagDeepFlavC+Jet_btagDeepFlavB)")
                                  .Define("goodcjets","Jet_pt > 20 && abs(Jet_eta) < 2.4 && CvsB > 0.4")
                                  .Define("goodbjets_l","Jet_pt > 20 && abs(Jet_eta) < 2.4 && Jet_btagDeepFlavB > 0.0614")
                                  .Define("goodbjets_m","Jet_pt > 20 && abs(Jet_eta) < 2.4 && Jet_btagDeepFlavB > 0.3093")
                                  .Define("goodbjets_t","Jet_pt > 20 && abs(Jet_eta) < 2.4 && Jet_btagDeepFlavB > 0.7221");
  }
  else if(year == 2017){
    df_S1_goodtau = df_S1_goodmuon.Filter("Sum( goodtaus ) >=1 ","Events with at least a goood tau")
                                  .Define("CvsB","Jet_btagDeepFlavC/(Jet_btagDeepFlavC+Jet_btagDeepFlavB)")
                                  .Define("goodcjets","Jet_pt > 20 && abs(Jet_eta) < 2.4 && CvsB > 0.4")
                                  .Define("goodbjets_l","Jet_pt > 20 && abs(Jet_eta) < 2.4 && Jet_btagDeepFlavB > 0.0521")
                                  .Define("goodbjets_m","Jet_pt > 20 && abs(Jet_eta) < 2.4 && Jet_btagDeepFlavB > 0.3033")
                                  .Define("goodbjets_t","Jet_pt > 20 && abs(Jet_eta) < 2.4 && Jet_btagDeepFlavB > 0.7489");
  }
  else if(year == 2018){
    df_S1_goodtau = df_S1_goodmuon.Filter("Sum( goodtaus ) >=1 ","Events with at least a goood tau")
                                  .Define("CvsB","Jet_btagDeepFlavC/(Jet_btagDeepFlavC+Jet_btagDeepFlavB)")
                                  .Define("goodcjets","Jet_pt > 20 && abs(Jet_eta) < 2.4 && CvsB > 0.4")
                                  .Define("goodbjets_l","Jet_pt > 20 && abs(Jet_eta) < 2.4 && Jet_btagDeepFlavB > 0.0494")
                                  .Define("goodbjets_m","Jet_pt > 20 && abs(Jet_eta) < 2.4 && Jet_btagDeepFlavB > 0.277")
                                  .Define("goodbjets_t","Jet_pt > 20 && abs(Jet_eta) < 2.4 && Jet_btagDeepFlavB > 0.7264");
  }
  auto df_S1_bjets = df_S1_goodtau.Filter("Sum( goodbjets_m ) >=1","Events with at least a b jet");
  auto df_S1_cjets = df_S1_bjets.Filter("Sum( goodcjets ) >=1","Events with at least a c jet");
 

  //auto d1 = df_S1_ttbar_categorize.Display({"ttbarAdditionalJetID","ttb","ttc","ttbb","ttcc","ttlf"},10000);
  //d1->Print();
  auto df_S1_tagged_jets = df_S1_cjets.Define("cjet_pt","Jet_pt[goodcjets]")
                                      .Define("cjet_eta","Jet_eta[goodcjets]")
                                      .Define("cjet_phi","Jet_phi[goodcjets]")
                                      .Define("cjet_mass","Jet_mass[goodcjets]")
                                      .Define("bjet_pt","Jet_pt[goodbjets_m]")
                                      .Define("bjet_eta","Jet_eta[goodbjets_m]")
                                      .Define("bjet_phi","Jet_phi[goodbjets_m]")
                                      .Define("bjet_mass","Jet_mass[goodbjets_m]")
                                      .Define("good_muon_pt","Muon_pt[goodmuons]")
                                      .Define("good_muon_eta","Muon_eta[goodmuons]")
                                      .Define("good_muon_phi","Muon_phi[goodmuons]")
                                      .Define("good_muon_mass","Muon_mass[goodmuons]")
                                      .Define("good_tau_pt","Tau_pt[goodtaus]")
                                      .Define("good_tau_eta","Tau_eta[goodtaus]")
                                      .Define("good_tau_phi","Tau_phi[goodtaus]")
                                      .Define("good_tau_mass","Tau_mass[goodtaus]")
                                      .Define("nbJet","Sum( goodbjets_m )")
                                      .Define("ncJet","Sum( goodcjets )")
                                      .Define("bjet1_pt","bjet_pt[0]")
                                      .Define("cjet1_pt","cjet_pt[0]");

  auto df_S1_ttbar_categorize = df_S1_tagged_jets.Define("ttbarAdditionalJetID","genTtbarId % 100")
                                                 .Define("nCjetsFromW","(genTtbarId % 100000) / 10000")
                                                 .Define("nBjetsFromW","(genTtbarId % 10000) / 1000")
                                                 .Define("nBjetsFromTop","(genTtbarId % 1000) / 100")
                                                 .Define("ttb","(ttbarAdditionalJetID/10)==5")
                                                 .Define("ttc","(ttbarAdditionalJetID/10)==4")
                                                 .Define("ttbb","(ttb % 10)>=3")
                                                 .Define("ttcc","(ttc % 10)>=3")
                                                 .Define("ttlf","ttbarAdditionalJetID == 0");
  auto calculate_dr = [](rvec_f a_eta, rvec_f b_eta, rvec_f a_phi, rvec_f b_phi) {
    float dr = DeltaR(a_eta[0], b_eta[0], a_phi[0], b_phi[0]);
    //cout<<dr<<"    "<<a_eta[0]<<"    "<<b_eta[0]<<"    "<<a_phi[0]<<"    "<<b_phi[0]<<endl;
    return dr;
  };

  auto calculate_M = [](rvec_f a_pt, rvec_f a_eta, rvec_f a_phi, rvec_f a_mass, rvec_f b_pt, rvec_f b_eta, rvec_f b_phi, rvec_f b_mass){
    TLorentzVector v1, v2;
    v1.SetPtEtaPhiM(a_pt[0],a_eta[0],a_phi[0],a_mass[0]);
    v2.SetPtEtaPhiM(b_pt[0],b_eta[0],b_phi[0],b_mass[0]);
    float invariant_mass = (v1 + v2).M();
    //auto invariant_mass = InvariantMasses(a_pt,a_eta,a_phi,a_mass,b_pt,b_eta,b_phi,b_mass);
    return invariant_mass;
  };

  auto df_S1_kinematics = df_S1_ttbar_categorize.Define("ctau_dR",calculate_dr,{"cjet_eta","good_tau_eta","cjet_phi","good_tau_phi"})
                                                .Define("ctau_M",calculate_M,{"cjet_pt","cjet_eta","cjet_phi","cjet_mass","good_tau_pt","good_tau_eta","good_tau_phi","good_tau_mass"});
  auto df_S1_ctau_dR = df_S1_kinematics.Filter("ctau_dR > 0.4","ctau dR > 0.4");
  auto df_S1_ttbar = df_S1_ctau_dR.Filter("genTtbarId > 0","ttbar events");
  
  //histograms 
  auto h_muon_pt = df_S1_goodmuon.Define("muon_pt","Muon_pt[goodmuons][0]").Histo1D({"h_muon_pt", "h_muon_pt", 100, 0, 100}, "muon_pt");
  auto h_n_selmuon = df_S1_goodmuon.Define("n_selmuon","Sum(goodmuons)").Histo1D({"h_n_selmuon", "h_n_selmuon", 5, 0, 5}, "n_selmuon");
  auto h_n_tau = df_S1_goodmuon.Histo1D({"h_n_tau", "h_n_tau", 5, 0, 5}, "nTau");
  auto h_n_seltau = df_S1_goodmuon.Define("ntaujets", "Sum(goodtaus)").Histo1D({"h_n_seltau", "h_n_seltau", 5, 0, 5}, "ntaujets");
  
  auto h_ctag = df_S1_goodtau.Histo1D({"h_ctag", "h_ctag", 100, 0, 1}, "CvsB");
  auto h_n_cjets = df_S1_goodtau.Define("ncjets","Sum( goodcjets )").Histo1D({"h_n_cjets", "h_n_cjets", 5, 0, 5}, "ncjets");
  auto h_n_bjets_l = df_S1_goodtau.Define("nbjets_l","Sum( goodbjets_l )").Histo1D({"h_n_bjets_l", "h_n_bjets_l", 5, 0, 5}, "nbjets_l");
  auto h_n_bjets_m = df_S1_goodtau.Define("nbjets_m","Sum( goodbjets_m )").Histo1D({"h_n_bjets_m", "h_n_bjets_m", 5, 0, 5}, "nbjets_m");
  auto h_n_bjets_t = df_S1_goodtau.Define("nbjets_t","Sum( goodbjets_t )").Histo1D({"h_n_bjets_t", "h_n_bjets_t", 5, 0, 5}, "nbjets_t");
  
  //tt histogram
  auto h_tag_bjet_pt= df_S1_tagged_jets.Histo1D({"h_tag_bjet1_pt","h_tag_bjet1_pt",60,0,300},"bjet1_pt");
  auto h_tag_cjet_pt= df_S1_tagged_jets.Histo1D({"h_tag_cjet1_pt","h_tag_cjet1_pt",60,0,300},"cjet1_pt");
  auto h_tag_muon_pt= df_S1_tagged_jets.Define("good_muon1_pt","good_muon_pt[0]").Histo1D({"h_tag_muon_pt","h_tag_muon_pt",60,0,300},"good_muon1_pt");
  auto h_tag_tau_pt= df_S1_tagged_jets.Define("good_tau1_pt","good_tau_pt[0]").Histo1D({"h_tag_tau_pt","h_tag_tau_pt",60,0,300},"good_tau1_pt");
  auto h_tag_nJet= df_S1_tagged_jets.Histo1D({"h_tag_nJet","h_tag_nJet",15,0,15},"nJet");
  auto h_tag_nbJet= df_S1_tagged_jets.Histo1D({"h_tag_nbJet","h_tag_nbJet",5,0,5},"nbJet");
  auto h_tag_ncJet= df_S1_tagged_jets.Histo1D({"h_tag_ncJet","h_tag_ncJet",5,0,5},"ncJet");
  auto h_tag_nMuon= df_S1_tagged_jets.Histo1D({"h_tag_nMuon","h_tag_nMuon",5,0,5},"nMuon");
  auto h_tag_nTau= df_S1_tagged_jets.Histo1D({"h_tag_nTau","h_tag_nTau",5,0,5},"nTau");
  
  //ttb histogram
  auto h_ttb_bjet_pt = df_S1_ttbar_categorize.Define("ttb_bjet_pt","Jet_pt[goodbjets_m && ttb][0]")
                                              .Histo1D({"h_ttb_bjet_pt","h_ttb_bjet_pt",60,0,300},"ttb_bjet_pt");
  auto h_ttb_cjet_pt = df_S1_ttbar_categorize.Define("ttb_cjet_pt","Jet_pt[goodcjets && ttb][0]")
                                              .Histo1D({"h_ttb_cjet_pt","h_ttb_cjet_pt",60,0,300},"ttb_cjet_pt");
  auto h_ttb_muon_pt = df_S1_ttbar_categorize.Define("ttb_muon_pt","Muon_pt[goodmuons && ttb][0]")
                                              .Histo1D({"h_ttb_muon_pt","h_ttb_muon_pt",60,0,300},"ttb_muon_pt");
  auto h_ttb_tau_pt = df_S1_ttbar_categorize.Define("ttb_tau_pt","Tau_pt[goodtaus && ttb][0]")
                                             .Histo1D({"h_ttb_tau_pt","h_ttb_tau_pt",60,0,300},"ttb_tau_pt");
  auto h_ttb_nJet = df_S1_ttbar_categorize.Define("ttb_nJet","Sum(goodjets && ttb)")
                                           .Histo1D({"h_ttb_nJet","h_ttb_nJet",15,0,15},"ttb_nJet");
  auto h_ttb_nbJet = df_S1_ttbar_categorize.Define("ttb_nbJet","Sum(goodbjets_m && ttb)")
                                            .Histo1D({"h_ttb_nbJet","h_ttb_nbJet",5,0,5},"ttb_nbJet");
  auto h_ttb_ncJet = df_S1_ttbar_categorize.Define("ttb_ncJet","Sum(goodcjets && ttb)")
                                            .Histo1D({"h_ttb_ncJet","h_ttb_ncJet",5,0,5},"ttb_ncJet");
  auto h_ttb_nMuon = df_S1_ttbar_categorize.Define("ttb_nMuon","Sum(goodmuons && ttb)")
                                            .Histo1D({"h_ttb_nMuon","h_ttb_nMuon",5,0,5},"ttb_nMuon");
  auto h_ttb_nTau = df_S1_ttbar_categorize.Define("ttb_nTau","Sum(goodtaus && ttb)")
                                           .Histo1D({"h_ttb_nTau","h_ttb_nTau",5,0,5},"ttb_nTau");
  
  //ttc histogram
  auto h_ttc_bjet_pt = df_S1_ttbar_categorize.Define("ttc_bjet_pt","Jet_pt[goodbjets_m && ttc][0]")
                                              .Histo1D({"h_ttc_bjet_pt","h_ttc_bjet_pt",60,0,300},"ttc_bjet_pt");
  auto h_ttc_cjet_pt = df_S1_ttbar_categorize.Define("ttc_cjet_pt","Jet_pt[goodcjets && ttc][0]")
                                              .Histo1D({"h_ttc_cjet_pt","h_ttc_cjet_pt",60,0,300},"ttc_cjet_pt");
  auto h_ttc_muon_pt = df_S1_ttbar_categorize.Define("ttc_muon_pt","Muon_pt[goodmuons && ttc][0]")
                                              .Histo1D({"h_ttc_muon_pt","h_ttc_muon_pt",60,0,300},"ttc_muon_pt");
  auto h_ttc_tau_pt = df_S1_ttbar_categorize.Define("ttc_tau_pt","Tau_pt[goodtaus && ttc][0]")
                                             .Histo1D({"h_ttc_tau_pt","h_ttc_tau_pt",60,0,300},"ttc_tau_pt");
  auto h_ttc_nJet = df_S1_ttbar_categorize.Define("ttc_nJet","Sum(goodjets && ttc)")
                                           .Histo1D({"h_ttc_nJet","h_ttc_nJet",15,0,15},"ttc_nJet");
  auto h_ttc_nbJet = df_S1_ttbar_categorize.Define("ttc_nbJet","Sum(goodbjets_m && ttc)")
                                            .Histo1D({"h_ttc_nbJet","h_ttc_nbJet",5,0,5},"ttc_nbJet");
  auto h_ttc_ncJet = df_S1_ttbar_categorize.Define("ttc_ncJet","Sum(goodcjets && ttc)")
                                            .Histo1D({"h_ttc_ncJet","h_ttc_ncJet",5,0,5},"ttc_ncJet");
  auto h_ttc_nMuon = df_S1_ttbar_categorize.Define("ttc_nMuon","Sum(goodmuons && ttc)")
                                            .Histo1D({"h_ttc_nMuon","h_ttc_nMuon",5,0,5},"ttc_nMuon");
  auto h_ttc_nTau = df_S1_ttbar_categorize.Define("ttc_nTau","Sum(goodtaus && ttc)")
                                           .Histo1D({"h_ttc_nTau","h_ttc_nTau",5,0,5},"ttc_nTau");
  
  //ttlf histogram
  auto h_ttlf_bjet_pt = df_S1_ttbar_categorize.Define("ttlf_bjet_pt","Jet_pt[goodbjets_m && ttlf][0]")
                                              .Histo1D({"h_ttlf_bjet_pt","h_ttlf_bjet_pt",60,0,300},"ttlf_bjet_pt");
  auto h_ttlf_cjet_pt = df_S1_ttbar_categorize.Define("ttlf_cjet_pt","Jet_pt[goodcjets && ttlf][0]")
                                              .Histo1D({"h_ttlf_cjet_pt","h_ttlf_cjet_pt",60,0,300},"ttlf_cjet_pt");
  auto h_ttlf_muon_pt = df_S1_ttbar_categorize.Define("ttlf_muon_pt","Muon_pt[goodmuons && ttlf][0]")
                                              .Histo1D({"h_ttlf_muon_pt","h_ttlf_muon_pt",60,0,300},"ttlf_muon_pt");
  auto h_ttlf_tau_pt = df_S1_ttbar_categorize.Define("ttlf_tau_pt","Tau_pt[goodtaus && ttlf][0]")
                                              .Histo1D({"h_ttlf_tau_pt","h_ttlf_tau_pt",60,0,300},"ttlf_tau_pt");
  auto h_ttlf_nJet = df_S1_ttbar_categorize.Define("ttlf_nJet","Sum(goodjets && ttlf)")
                                           .Histo1D({"h_ttlf_nJet","h_ttlf_nJet",15,0,15},"ttlf_nJet");
  auto h_ttlf_nbJet = df_S1_ttbar_categorize.Define("ttlf_nbJet","Sum(goodbjets_m && ttlf)")
                                            .Histo1D({"h_ttlf_nbJet","h_ttlf_nbJet",5,0,5},"ttlf_nbJet");
  auto h_ttlf_ncJet = df_S1_ttbar_categorize.Define("ttlf_ncJet","Sum(goodcjets && ttlf)")
                                            .Histo1D({"h_ttlf_ncJet","h_ttlf_ncJet",5,0,5},"ttlf_ncJet");
  auto h_ttlf_nMuon = df_S1_ttbar_categorize.Define("ttlf_nMuon","Sum(goodmuons && ttlf)")
                                            .Histo1D({"h_ttlf_nMuon","h_ttlf_nMuon",5,0,5},"ttlf_nMuon");
  auto h_ttlf_nTau = df_S1_ttbar_categorize.Define("ttlf_nTau","Sum(goodtaus && ttlf)")
                                           .Histo1D({"h_ttlf_nTau","h_ttlf_nTau",5,0,5},"ttlf_nTau");

  //ctau histogram
  auto h_ctau_dR = df_S1_ctau_dR.Histo1D({"h_ctau_dR","h_ctau_dR",25,0,5},"ctau_dR");
  auto h_ctau_M = df_S1_ctau_dR.Histo1D({"h_ctau_M","h_ctau_M",100,0,500},"ctau_M");
  auto h_ttb_ctau_dR = df_S1_ttbar.Filter("ttb","ttb").Histo1D({"h_ttb_ctau_dR","h_ttb_ctau_dR",25,0,5},"ctau_dR");
  auto h_ttb_ctau_M = df_S1_ttbar.Filter("ttb","ttb").Histo1D({"h_ttb_ctau_M","h_ttb_ctau_M",100,0,500},"ctau_M");
  auto h_ttc_ctau_dR = df_S1_ttbar.Filter("ttc","ttc").Histo1D({"h_ttc_ctau_dR","h_ttc_ctau_dR",25,0,5},"ctau_dR");
  auto h_ttc_ctau_M = df_S1_ttbar.Filter("ttc","ttc").Histo1D({"h_ttc_ctau_M","h_ttc_ctau_M",100,0,500},"ctau_M");
  auto h_ttlf_ctau_dR = df_S1_ttbar.Filter("ttlf","ttlf").Histo1D({"h_ttlf_ctau_dR","h_ttlf_ctau_dR",25,0,5},"ctau_dR");
  auto h_ttlf_ctau_M = df_S1_ttbar.Filter("ttlf","ttlf").Histo1D({"h_ttlf_ctau_M","h_ttlf_ctau_M",100,0,500},"ctau_M");
  
  //df_S1_goodmuon.Snapshot("tree", "f.root");
  TFile f("f.root", "recreate");
  
  plot( h_muon_pt, "h_muon_pt");
  plot( h_n_selmuon, "h_n_selmuon");
  plot( h_n_tau, "h_n_tau");
  plot( h_n_seltau, "h_n_seltau");
  plot( h_ctag, "h_ctag");
  plot( h_n_cjets, "h_n_cjets");
  plot( h_n_bjets_l, "h_n_bjets_l");
  plot( h_n_bjets_m, "h_n_bjets_m");
  plot( h_n_bjets_t, "h_n_bjets_t");
  plot( h_tag_bjet_pt, "h_tag_bjet_pt");
  plot( h_tag_cjet_pt, "h_tag_cjet_pt");
  plot( h_tag_muon_pt, "h_tag_muon_pt");
  plot( h_tag_tau_pt, "h_tag_tau_pt");
  plot( h_tag_nJet, "h_tag_nJet");
  plot( h_tag_nbJet, "h_tag_nbJet");
  plot( h_tag_ncJet, "h_tag_ncJet");
  plot( h_tag_nMuon, "h_tag_nMuon");
  plot( h_tag_nTau, "h_tag_nTau");
  plot( h_ttb_bjet_pt, "h_ttb_bjet_pt");
  plot( h_ttb_cjet_pt, "h_ttb_cjet_pt");
  plot( h_ttb_muon_pt, "h_ttb_muon_pt");
  plot( h_ttb_tau_pt, "h_ttb_tau_pt");
  plot( h_ttb_nJet, "h_ttb_nJet");
  plot( h_ttb_nbJet, "h_ttb_nbJet");
  plot( h_ttb_ncJet, "h_ttb_ncJet");
  plot( h_ttb_nMuon, "h_ttb_nMuon");
  plot( h_ttb_nTau, "h_ttb_nTau");
  plot( h_ttc_bjet_pt, "h_ttc_bjet_pt");
  plot( h_ttc_cjet_pt, "h_ttc_cjet_pt");
  plot( h_ttc_muon_pt, "h_ttc_muon_pt");
  plot( h_ttc_tau_pt, "h_ttc_tau_pt");
  plot( h_ttc_nJet, "h_ttc_nJet");
  plot( h_ttc_nbJet, "h_ttc_nbJet");
  plot( h_ttc_ncJet, "h_ttc_ncJet");
  plot( h_ttc_nMuon, "h_ttc_nMuon");
  plot( h_ttc_nTau, "h_ttc_nTau");
  plot( h_ttlf_bjet_pt, "h_ttlf_bjet_pt");
  plot( h_ttlf_cjet_pt, "h_ttlf_cjet_pt");
  plot( h_ttlf_muon_pt, "h_ttlf_muon_pt");
  plot( h_ttlf_tau_pt, "h_ttlf_tau_pt");
  plot( h_ttlf_nJet, "h_ttlf_nJet");
  plot( h_ttlf_nbJet, "h_ttlf_nbJet");
  plot( h_ttlf_ncJet, "h_ttlf_ncJet");
  plot( h_ttlf_nMuon, "h_ttlf_nMuon");
  plot( h_ttlf_nTau, "h_ttlf_nTau");
  plot( h_ctau_dR, "h_ctau_dR");
  plot( h_ctau_M, "h_ctau_M");
  plot( h_ttb_ctau_dR, "h_ttb_ctau_dR");
  plot( h_ttb_ctau_M, "h_ttb_ctau_M");
  plot( h_ttc_ctau_dR, "h_ttc_ctau_dR");
  plot( h_ttc_ctau_M, "h_ttc_ctau_M");
  plot( h_ttlf_ctau_dR, "h_ttlf_ctau_dR");
  plot( h_ttlf_ctau_M, "h_ttlf_ctau_M");
  f.Close();

  //auto report = df_S1_tagged_jets.Report();
  auto report = df_S1_ttbar.Report();
  report->Print();

}

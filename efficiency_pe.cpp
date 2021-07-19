/**
 * Example: Selection efficiency from CAFs.
 *
 * Author: A. Mastbaum <mastbaum@physics.rutgers.edu>, 2020/10/21
 * Contributor: P. Englezos <p.englezos@physics.rutgers.edu> 
 *
 * This example loops through events in a CAF file, with the added field for
 * hadronic energy deposited in the inner collar region around an inactive central
 * detector module (defined by the central module in the 5x7 array plus 30 cm around it in x and z).
 *
 * It creates two sets of plots of the acceptance, for selecting ND-LAr events
 * imposing the following sets of cuts:
 * 
 * 1) Muon neutrino charged-current interaction (numuCC) occured and 
 *     its vertex was inside the fiducial volume of the detector (when all modules are active).  
 *
 * 2) (1), muon was reconstructed and the hadronic activity was contained in the region inside from the outer walls of the 5x7 module array by 30 cm (outer collar region).
 *
 * 3) (1), numuCC's vertex was *outside* the anti-fiducial volume, the hadronic activity was contained between the two collar regions (inner and outer) and muon's track did not end inside the  *    inactive module. 
 *
 *
 * We then compute acceptances for (2)/(1), and for (3)/(2), the loss in acceptance when additionally fiducializing around the center module.
 *
 * Usage:
 *
 *   ./efficiency filelist.txt
 *
 */

#include <fstream>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"

void efficiency(std::string infiles) {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  std::vector<std::string> filenames;
  std::ifstream myfile(infiles);
  std::copy(std::istream_iterator<std::string>(myfile),
            std::istream_iterator<std::string>(),
            std::back_inserter(filenames));
  std::cout << "File Count: " << filenames.size() << std::endl;

  // Histograms
  TFile* fout = TFile::Open("hists.root", "recreate");
  
  TH1F* h_enu_orig = new TH1F("h_enu_orig", ";E_{#nu} (GeV);Efficiency", 400, 0, 10);
  h_enu_orig->SetDirectory(0);
  h_enu_orig->Sumw2();
  TH1F* h_enu_veto = (TH1F*) h_enu_orig->Clone("h_enu_veto");
  TH1F* h_enu_dead = (TH1F*) h_enu_orig->Clone("h_enu_dead");

  std::vector<TH2F*> h_q0q3_orig(10);
  std::vector<TH2F*> h_q0q3_veto(10);
  std::vector<TH2F*> h_q0q3_dead(10);
  for (int i=0; i<10;i++){
   char nameo[50];
   char namev[50];
   char named[50];
   snprintf(nameo, 50, "h_q0q3_orig_%lu", i);
   snprintf(namev, 50, "h_q0q3_veto_%lu", i);
   snprintf(named, 50, "h_q0q3_dead_%lu", i);
   h_q0q3_orig[i] = new TH2F(nameo, ";q^{3} (GeV/c);q^{0} (GeV);Efficiency", 200, 0, 10, 200, 0, 10);
   h_q0q3_orig[i] ->SetDirectory(0);
   h_q0q3_orig[i] ->Sumw2();
   h_q0q3_veto[i] = new TH2F(namev, ";q^{3} (GeV/c);q^{0} (GeV);Efficiency", 200, 0, 10, 200, 0, 10);
   h_q0q3_veto[i] ->SetDirectory(0);
   h_q0q3_veto[i] ->Sumw2();
   h_q0q3_dead[i] = new TH2F(named, ";q^{3} (GeV/c);q^{0} (GeV);Efficiency", 200, 0, 10, 200, 0, 10);
   h_q0q3_dead[i] ->SetDirectory(0);
   h_q0q3_dead[i] ->Sumw2(); 
 }

  std::vector<TH2F*> h_KE_orig(10);
  std::vector<TH2F*> h_KE_veto(10);
  std::vector<TH2F*> h_KE_dead(10);
  for (int i=0; i<10;i++){
    char nameKEo[50];
    char nameKEv[50];
    char nameKEd[50];
    snprintf(nameKEo, 50, "h_KE_orig_%lu", i);
    snprintf(nameKEv, 50, "h_KE_veto_%lu", i);
    snprintf(nameKEd, 50, "h_KE_dead_%lu", i);
    h_KE_orig[i] = new TH2F(nameKEo, ";cos(#theta) ;KE (GeV);Efficiency", 200, - 1, 1, 200, 0, 10);
    h_KE_orig[i]->SetDirectory(0);
    h_KE_orig[i]->Sumw2();
    h_KE_veto[i] = new TH2F(nameKEv, ";cos(#theta) ;KE (GeV);Efficiency", 200, - 1, 1, 200, 0, 10);
    h_KE_veto[i]->SetDirectory(0);
    h_KE_veto[i]->Sumw2(); 
    h_KE_dead[i] = new TH2F(nameKEd, ";cos(#theta) ;KE (GeV);Efficiency", 200, - 1, 1, 200, 0, 10);
    h_KE_dead[i]->SetDirectory(0);
    h_KE_dead[i]->Sumw2();
  } 


  TH2F* h_q0q3sqr_orig = new TH2F("h_q0q3sqr_orig", ";(q^{3})^{2} (GeV/c)^{2};E_{#nu} (GeV);Efficiency", 200, 0, 10, 200, 0, 10);
  h_q0q3sqr_orig->SetDirectory(0);
  h_q0q3sqr_orig->Sumw2();
  TH2F* h_q0q3sqr_veto = (TH2F*) h_q0q3sqr_orig->Clone("h_q0q3sqr_veto");
  TH2F* h_q0q3sqr_dead = (TH2F*) h_q0q3sqr_orig->Clone("h_q0q3sqr_dead");

  TH2F* h_lost = new TH2F("h_lost", ";E_{#nu} (GeV);E_{lost} (MeV);Probability/bin", 200, 0, 20, 200, 0, 500);
  h_lost->SetDirectory(0);

  std::vector<TH2F*> h_vtx_xz(10);
  std::vector<TH2F*> h_vtx_xz_veto(10);
  std::vector<TH2F*> h_vtx_xz_dead(10);
  std::vector<TH2F*> h_vtx_xz_lost(10);
  for (int i=0; i<10;i++){
   char namevo[50];
   char namevv[50];
   char namevd[50];
   char namevl[50];
   snprintf(namevo, 50, "h_vtx_xz_%lu", i);
   snprintf(namevv, 50, "h_vtx_xz_veto_%lu", i);
   snprintf(namevd, 50, "h_vtx_xz_dead_%lu", i);
   snprintf(namevl, 50, "h_vtx_xz_lost_%lu", i);
   h_vtx_xz[i] = new TH2F(namevo, "Vertex XZ", 200, -350, 350, 200, 0, 500);
   h_vtx_xz[i] ->SetDirectory(0);
   h_vtx_xz[i] ->Sumw2();
   h_vtx_xz_veto[i] = new TH2F(namevv, "Vertex XZ", 200, -350, 350, 200, 0, 500);
   h_vtx_xz_veto[i] ->SetDirectory(0);
   h_vtx_xz_veto[i] ->Sumw2();
   h_vtx_xz_dead[i] = new TH2F(namevd, "Vertex XZ", 200, -350, 350, 200, 0, 500);
   h_vtx_xz_dead[i] ->SetDirectory(0);
   h_vtx_xz_dead[i] ->Sumw2();
   h_vtx_xz_lost[i] = new TH2F(namevl, "Vertex XZ", 200, -350, 350, 200, 0, 500);
   h_vtx_xz_lost[i] ->SetDirectory(0);
   h_vtx_xz_lost[i] ->Sumw2(); 
}

  // Loop through files
  for (auto& infile : filenames) {
    std::cout << "File: " << infile << std::endl;

    // Set up the input tree
    TFile* f = TFile::Open(infile.c_str());
    if (!(f && f->IsOpen() && !f->IsZombie())) {
      std::cout << "Bad file! " << infile << std::endl;
      continue;
    }
    TTree* caf = (TTree*) f->Get("caf");
    if (!(caf && caf->GetEntries())) {
      std::cout << "Bad file! " << infile << std::endl;
      continue;
    }

    int isCC, nuPDG, muon_contained, muon_tracker, muon_ecal, muon_exit, geg=1;
    double Ev, pxv, pyv, pzv, El, pxl, pyl, pzl, vx, vy, vz;
    double Eveto, Edeadveto, Edead, lepKE;
    bool flagr;

    caf->SetBranchStatus("*", 0);
    caf->SetBranchStatus("nuPDG", 1);
    caf->SetBranchStatus("isCC", 1);
    caf->SetBranchStatus("muon_contained", 1);
    caf->SetBranchStatus("muon_tracker", 1);
    caf->SetBranchStatus("muon_ecal", 1);
    caf->SetBranchStatus("muon_exit", 1);
    caf->SetBranchStatus("Ev", 1);
    caf->SetBranchStatus("NuMomX",1); 
    caf->SetBranchStatus("NuMomY", 1);
    caf->SetBranchStatus("NuMomZ", 1);
    caf->SetBranchStatus("LepE", 1);
    caf->SetBranchStatus("LepMomX", 1);
    caf->SetBranchStatus("LepMomY", 1);
    caf->SetBranchStatus("LepMomZ", 1);
    caf->SetBranchStatus("vtx_x", 1);
    caf->SetBranchStatus("vtx_y", 1);
    caf->SetBranchStatus("vtx_z", 1);
    caf->SetBranchStatus("Ehad_veto", 1); 
    caf->SetBranchStatus("dead_Ehad_veto", 1);
    caf->SetBranchStatus("inner_dead_Ehad_veto", 1);
    caf->SetBranchStatus("flagr", 1);

    caf->SetBranchAddress("nuPDG", &nuPDG);
    caf->SetBranchAddress("isCC", &isCC);
    caf->SetBranchAddress("Ev", &Ev);
    caf->SetBranchAddress("muon_contained", &muon_contained);
    caf->SetBranchAddress("muon_tracker", &muon_tracker);
    caf->SetBranchAddress("muon_ecal", &muon_ecal);
    caf->SetBranchAddress("muon_exit", &muon_exit);
    caf->SetBranchAddress("NuMomX", &pxv);
    caf->SetBranchAddress("NuMomY", &pyv);
    caf->SetBranchAddress("NuMomZ", &pzv);
    caf->SetBranchAddress("LepE", &El);
    caf->SetBranchAddress("LepMomX", &pxl);
    caf->SetBranchAddress("LepMomY", &pyl);
    caf->SetBranchAddress("LepMomZ", &pzl);
    caf->SetBranchAddress("vtx_x", &vx);
    caf->SetBranchAddress("vtx_y", &vy);
    caf->SetBranchAddress("vtx_z", &vz);
    caf->SetBranchAddress("Ehad_veto", &Eveto);
    caf->SetBranchAddress("dead_Ehad_veto", &Edeadveto);
    caf->SetBranchAddress("inner_dead_Ehad_veto", &Edead);
    caf->SetBranchStatus("flagr", &flagr);

     std::cout << caf->GetEntries() << std::endl;
    // Event loop
    for (long i=0; i<caf->GetEntries(); i++) {
      caf->GetEntry(i);

      // numuCC inclusive
      if (!(isCC == 1 && nuPDG == 14)) {
        continue;
      }

      // Vertex in the active and fiducial volume
      bool vtx_active = (vx < 350 && vx > -350 &&
                         vy < 150 && vy > -150 &&
                         vz < 500 && vz >    0);

      bool vtx_fiducial = (vx < 300 && vx > -300 &&
                           vy < 100 && vy > -100 &&
                           vz < 350 && vz >   50);

      if (!(vtx_active && vtx_fiducial)) {
        continue;
      }

      bool muon_reco = muon_contained || muon_tracker || muon_ecal;

      bool collar_edep = (Eveto > 10);

      // Vertex within 50 cm (150 cm upstream) anti-fiducial volume around
      // the central module of the 5x7 array.
      //bool vtx_dead_antifid = (vx > -100 && vx < 100 &&
        //                       vy > -150 && vy < 150 &&
          //                     vz >   50 && vz < 350);
     //Update_To_the_left
      bool vtx_dead_antifid = (vx > -350 && vx < -200 &&
                                 vy > -150 && vy < 150 &&
                                 vz >   50 && vz < 350);

      bool dead_collar_edep = (Edeadveto > 10);

      // Vertex inside the central module.
      //bool vtx_dead = (vx >  -50 && vx <  50 &&
        //               vy > -150 && vy < 150 &&
          //             vz >  200 && vz < 300);
      
      //Update_To_the_left
      bool vtx_dead = (vx >  -350 && vx <  -250 &&
                       vy > -150 && vy < 150 &&
                       vz >  200 && vz < 300);

      // Kinematic variables
      float m = sqrt(El*El - pxl*pxl - pyl*pyl - pzl*pzl); 
      lepKE = El - m; 
      float q0 = Ev - El;
      float costheta = pzl/sqrt(pxl*pxl + pyl*pyl + pzl*pzl);
      float dpx = pxv - pxl;
      float dpy = pyv - pyl;
      float dpz = pzv - pzl;
      float q3 = sqrt(dpx*dpx + dpy*dpy + dpz*dpz);

      // Fill histograms for selected events
     if (Ev > 0 && Ev < 10) {
        //All events with FV vertex
        h_q0q3sqr_orig->Fill(q3*q3,Ev);
        h_enu_orig->Fill(Ev);
       
        // Hadronic collar veto
           if (!collar_edep && muon_reco) {
              h_enu_veto->Fill(Ev);
              h_q0q3sqr_veto->Fill(q3*q3,Ev);
           }
        
          // Dead module veto: fiducial region and collar around central module
                if (!(vtx_dead_antifid || dead_collar_edep || collar_edep) && muon_reco && !flagr) {
                    h_enu_dead->Fill(Ev);
                    h_q0q3sqr_dead->Fill(q3*q3,Ev);
               
                   // Accepted events depositing energy in the dead module
                          if (Edead > 0) {
                              h_lost->Fill(Ev, Edead);
                }  
            }
       } 

      for (long i=0; i<10; i++) {
        if (Ev >= i && Ev <= i+1) {
        // All events with FV vertex
        h_q0q3_orig[i]->Fill(q3, q0);
        h_vtx_xz[i]->Fill(vx, vz);
        h_KE_orig[i]->Fill(costheta, lepKE);

        // Hadronic collar veto
        if (!collar_edep && muon_reco){
          h_q0q3_veto[i]->Fill(q3, q0);
          h_vtx_xz_veto[i]->Fill(vx, vz);
          h_KE_veto[i]->Fill(costheta, lepKE);
        }

        // Dead module veto: fiducial region and collar around central module
        if (!(vtx_dead_antifid || dead_collar_edep || collar_edep) && muon_reco && !flagr) {
          h_q0q3_dead[i]->Fill(q3, q0);
          h_vtx_xz_dead[i]->Fill(vx, vz);
          h_KE_dead[i]->Fill(costheta, lepKE);
             
          // Accepted events depositing energy in the dead module
          if (Edead > 0) {
            h_vtx_xz_lost[i]->Fill(vx, vz);
          }
        }
      }
     }
    }
    f->Close();
  }

  //std::cout << "Events: "
    //        << h_enu_dead->GetEntries() << " / "
    //      << h_enu_veto->GetEntries() << " / "
   //    << h_enu_orig->GetEntries() << std::endl;

  //Update
  //veto
  assert(TEfficiency::CheckConsistency(*h_q0q3sqr_veto, *h_q0q3sqr_orig));
  TEfficiency* eff_veto_q0q3sqr = new TEfficiency(*h_q0q3sqr_veto, *h_q0q3sqr_orig);
  eff_veto_q0q3sqr->SetDirectory(0);
  eff_veto_q0q3sqr->SetName("eff_veto_q0q3sqr");
  TCanvas* c16 = new TCanvas();
  gStyle->SetPalette(60);
  eff_veto_q0q3sqr->Draw("colz");
  c16->SaveAs("eff_veto_q0q3sqr.pdf");

  std::vector<TEfficiency*> eff_veto_KE(10);
  for (long i=0; i<10; i++) {
  assert(TEfficiency::CheckConsistency(*h_KE_veto[i], *h_KE_orig[i]));
  char nameeffvKE[50];
  snprintf(nameeffvKE, 50, "eff_veto_KE_%lu", i);
  eff_veto_KE[i] = new TEfficiency(*h_KE_veto[i], *h_KE_orig[i]);
  eff_veto_KE[i]->SetDirectory(0);
  eff_veto_KE[i]->SetName(nameeffvKE);
  TCanvas* c20 = new TCanvas();
  gStyle->SetPalette(60);
  eff_veto_KE[i]->Draw("colz");
  c20->SaveAs(TString(eff_veto_KE[i]->GetName()) + ".pdf");
  }
  
  //dead
  assert(TEfficiency::CheckConsistency(*h_q0q3sqr_dead, *h_q0q3sqr_veto));
  TEfficiency* eff_dead_q0q3sqr = new TEfficiency(*h_q0q3sqr_dead, *h_q0q3sqr_veto);
  eff_dead_q0q3sqr->SetDirectory(0);
  eff_dead_q0q3sqr->SetName("eff_dead_q0q3sqr");
  TCanvas* c17 = new TCanvas();
  gStyle->SetPalette(60);
  eff_dead_q0q3sqr->Draw("colz");
  c17->SaveAs("eff_dead_q0q3sqr.pdf");

  std::vector<TEfficiency*> eff_dead_KE(10);
  for (long i=0; i<10; i++) {
  assert(TEfficiency::CheckConsistency(*h_KE_dead[i], *h_KE_veto[i]));
  char nameeffdKE[50];
  snprintf(nameeffdKE, 50, "eff_dead_KE_%lu", i);
  eff_dead_KE[i] = new TEfficiency(*h_KE_dead[i], *h_KE_veto[i]);
  eff_dead_KE[i]->SetDirectory(0);
  eff_dead_KE[i]->SetName(nameeffdKE);
  TCanvas* c21 = new TCanvas();
  gStyle->SetPalette(60);
  eff_dead_KE[i]->Draw("colz");
  c21->SaveAs(TString(eff_dead_KE[i]->GetName()) + ".pdf");
  }
  
  // Efficiency: (collar veto)/total

  assert(TEfficiency::CheckConsistency(*h_enu_veto, *h_enu_orig));
  TEfficiency* eff_veto = new TEfficiency(*h_enu_veto, *h_enu_orig);
  eff_veto->SetDirectory(0);
  eff_veto->SetName("eff_veto");
  TCanvas* c1 = new TCanvas("eff_veto");
  eff_veto->Draw("ap");
  c1->SaveAs("eff_veto.pdf");
  

  std::vector<TEfficiency*> eff_veto_q0q3(10);
  for (long i=0; i<10; i++) { 
  assert(TEfficiency::CheckConsistency(*h_q0q3_veto[i], *h_q0q3_orig[i]));
  char nameeffv[50];
  snprintf(nameeffv, 50, "eff_veto_q0q3_%lu", i);
  eff_veto_q0q3[i] = new TEfficiency(*h_q0q3_veto[i], *h_q0q3_orig[i]);
  eff_veto_q0q3[i]->SetDirectory(0);
  eff_veto_q0q3[i]->SetName(nameeffv);
  TCanvas* c3 = new TCanvas();
  gStyle->SetPalette(60);
  eff_veto_q0q3[i]->Draw("colz");
  c3->SaveAs(TString(eff_veto_q0q3[i]->GetName()) + ".pdf");
  }

  // Efficiency: (collar + dead module veto)/(collar veto only)
  assert(TEfficiency::CheckConsistency(*h_enu_dead, *h_enu_veto));
  TEfficiency* eff_dead = new TEfficiency(*h_enu_dead, *h_enu_veto);
  eff_dead->SetDirectory(0);
  eff_dead->SetName("eff_dead");
  TCanvas* c2 = new TCanvas();
  eff_dead->Draw("ap");
  c2->SaveAs("eff_dead.pdf");
  
    
  std::vector<TEfficiency*> eff_dead_q0q3(10);
  for (long i=0; i<10; i++) {
  assert(TEfficiency::CheckConsistency(*h_q0q3_dead[i], *h_q0q3_veto[i]));
  char nameeffd[50];
  snprintf(nameeffd, 50, "eff_dead_qoq3_%lu", i);
  eff_dead_q0q3[i] = new TEfficiency(*h_q0q3_dead[i], *h_q0q3_veto[i]);
  eff_dead_q0q3[i]->SetDirectory(0);
  eff_dead_q0q3[i]->SetName(nameeffd);
  TCanvas* c4 = new TCanvas();
  gStyle->SetPalette(60);
  eff_dead_q0q3[i]->Draw("colz");
  c4->SaveAs(TString(eff_dead_q0q3[i]->GetName()) + ".pdf");
  }
 
  TCanvas* c5 = new TCanvas();
  h_lost->Scale(1.0 / h_lost->Integral());
  gStyle->SetPalette(60);
  h_lost->Draw("colz");
  c5->SetLogz();
  c5->SetRightMargin(0.13);
  c5->SaveAs("h_lost.pdf");


 for (long i=0; i<10; i++) {
  TCanvas* c6 = new TCanvas();
  gStyle->SetPalette(60);
  h_vtx_xz[i] ->Draw("colz");
  c6->SaveAs(TString(h_vtx_xz[i]->GetName()) + ".pdf");

  TCanvas* c7 = new TCanvas();
  gStyle->SetPalette(60);
  h_vtx_xz_veto[i] -> Draw("colz");
  c7->SaveAs(TString(h_vtx_xz_veto[i]->GetName()) + ".pdf");

  TCanvas* c8 = new TCanvas();
  gStyle->SetPalette(60);
  h_vtx_xz_dead[i] ->Draw("colz");
  c8->SaveAs(TString(h_vtx_xz_dead[i]->GetName()) + ".pdf");

  TCanvas* c9 = new TCanvas();
  gStyle->SetPalette(60);
  h_vtx_xz_lost[i] ->Draw("colz");
  c9->SaveAs(TString(h_vtx_xz_lost[i]->GetName()) + ".pdf");
}

  // Write
  fout->cd();

 for (long i=0; i<10; i++) { 
 h_q0q3_orig[i]->Write();
 h_q0q3_veto[i]->Write();
 h_q0q3_dead[i]->Write();
 eff_veto_q0q3[i]->Write();
 eff_dead_q0q3[i]->Write();
 h_vtx_xz[i]->Write();
 h_vtx_xz_veto[i]->Write();
 h_vtx_xz_dead[i]->Write();
 h_vtx_xz_lost[i]->Write();
 h_KE_orig[i]->Write();
 h_KE_veto[i]->Write(); 
 h_KE_dead[i]->Write();
 eff_veto_KE[i]->Write();
 eff_dead_KE[i]->Write();
}
  h_enu_orig->Write();
  h_enu_veto->Write();
  h_enu_dead->Write();
  eff_veto->Write();
  eff_dead->Write();
  h_lost->Write(); 
  fout->Close();
}


int main(int argc, char* argv[]) {
   efficiency(argv[1]);
   return 0;
}


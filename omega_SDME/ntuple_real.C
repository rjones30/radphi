#define ntuple_real_cxx
// The class definition in ntuple_real.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("ntuple_real.C")
// root> T->Process("ntuple_real.C","some options")
// root> T->Process("ntuple_real.C+")
//

#define DUMP_KINE_TREE 0
#define MAKE_SKIM_TREE 1

#include <TH2.h>
#include <TStyle.h>

#include "ntuple_real.h"

//#include "sqr.h"
inline unsigned int sqr(unsigned int x) { return x*x; }
inline Int_t sqr(Int_t x) { return x*x; }
inline Float_t sqr(Float_t x) { return x*x; }
inline Double_t sqr(Double_t x) { return x*x; }


void ntuple_real::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void ntuple_real::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   int n=0;
   TString name;
   TString title;
   for (int run=7415; run < 8720; run += 30) {
      name.Form("m3g%d", n);
      title.Form("Elost vs mass3g, runs %d - %d", run, run + 29);
      h2mass3g[n] = new TH2D(name, title, 150, 0, 1.5, 30, -1.0, 2.0);
      h2mass3g[n]->Sumw2();
      GetOutputList()->Add(h2mass3g[n]);
      name.Form("dphi%d", n);
      title.Form("Elost vs dphi, runs %d - %d", run, run + 29);
      h2dphi[n] = new TH2D(name, title, 90, -180, 180, 30, -1.0, 2.0);
      h2dphi[n]->Sumw2();
      GetOutputList()->Add(h2dphi[n]);
      name.Form("E3g%d", n);
      title.Form("Emiss vs mass3g, runs %d - %d", run, run + 29);
      h2E3g[n] = new TH2D(name, title, 150, 0, 1.5, 30, -1.0, 2.0);
      h2E3g[n]->Sumw2();
      GetOutputList()->Add(h2E3g[n]);
      name.Form("heli%d", n);
      title.Form("helicity frame cos(theta) vs phi, runs %d - %d", 
                 run, run + 29);
      h2heli[n] = new TH2D(name, title, 100, -1.0, 1.0, 90, -3.1416, 3.1416);
      h2heli[n]->Sumw2();
      GetOutputList()->Add(h2heli[n]);
      n++;
   }
   tkine = 0;
   tskim = 0;
   dumpfile = 0;
   skimfile = 0;
}

Bool_t ntuple_real::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   // fReader.SetEntry(entry);
   if (fChain == 0) 
      return kFALSE;
   int nb = fChain->GetEntry(entry);

   // omega selection cuts
   if (nfrwd != 3 || nphot != 3)
      return kFALSE; // cuts 5, 7

   double mleast = 999;
   double ptot[4] = {0};
   for (int i=0; i < nphot; ++i) {
      if (pvect[i][0] < 0.050)
         return kFALSE; // cut 1
      for (int mu=0; mu < 4; ++mu) {
         ptot[mu] += pvect[i][mu];
      }
      double pperp = sqrt(sqr(pvect[i][1]) + sqr(pvect[i][2]));
      double theta = atan2(pperp, pvect[i][3]) * 180/3.1415926;
      double A = 0.13;
      double B = 1.033;
      if (pvect[i][0] < 0.5 && pvect[i][0] < B - A * theta)
         return kFALSE; // cut 6
      for (int j=i+1; j < nphot; ++j) {
         double p[4] = {0};
         for (int mu=0; mu < 4; ++mu) {
            p[mu] = pvect[i][mu] + pvect[j][mu];
         }
         double minv = sqrt(sqr(p[0]) - sqr(p[1]) - sqr(p[2]) - sqr(p[3]));
         if (minv < mleast)
            mleast = minv;
      }
   }

   if (ptot[0] < 3.0)
      return kFALSE; // cut 2

   if (mleast < 0.10 || mleast > 0.18)
      return kFALSE; // cut 8

   if (nrec != 1)
      return kFALSE; // cut 3

   int ncp = 0;
   for (int i=0; i < nhcpv; ++i) {
      int h0 = t1cpv[i];
      for (int h=0; h < ntcpv[i]; ++h) {
         if (fabs(le[h0 + h] - trec0) < 3.0)
            ++ncp;
      }
   }
   if (ncp > 0)
      return kFALSE; // cut 4

   double m3g = sqrt(sqr(ptot[0]) - sqr(ptot[1]) - sqr(ptot[2]) - sqr(ptot[3]));
   double dphi = (atan2(ptot[2], ptot[1]) - phirec[0]) * 180/3.1415926;
   dphi += 180;
   while (dphi < -180) 
      dphi += 360;
   while (dphi > +180) 
      dphi -= 360;

   double naccid = 0;
   double ntrues = 0;
   int nhist = (runno - 7415) / 30;
   for (int i=0; i < ncoin; ++i) {
#if DO_CPV_AGAINST_TAGGER
      int nocp = 1;
      for (int i=0; nocp && i < nhcpv; ++i) {
         int h0 = t1cpv[i];
         for (int h=0; nocp && h < ntcpv[i]; ++h) {
            if (fabs(le[h0 + h] - trec0) < 3.0)
               nocp = 0;
         }
      }
      if (!nocp)
         continue;
#endif
      double tagw = 0;
      if (fabs(trec0 - cotime[i] - 2.0) < 3.0) {
         ntrues += 1;
         tagw = 1;
      }
      else if (fabs(trec0 - cotime[i] - 8.0) < 3.0) {
         naccid += 1;
         tagw = -1;
      }
      double Elost = coenergy[i] - Efrwd;
      double Emiss = coenergy[i] - ptot[0];
      double tabs = sqr(ptot[1]) + sqr(ptot[2]) +
                    sqr(ptot[3] - coenergy[i]) -
                    sqr(ptot[0] - coenergy[i]);
      h2mass3g[nhist]->Fill(m3g, Elost, tagw);
      h2dphi[nhist]->Fill(dphi, Elost, tagw);
      h2E3g[nhist]->Fill(m3g, Emiss, tagw);
      TLorentzVector ppi0 = GetPpi0Helicity(coenergy[i]);
#if MAKE_SKIM_TREE
      if (tskim == 0) {
         if (skimfile == 0)
            skimfile = new TFile("/tmp/proof/jonesrt/ntuple_real_skim.root", "recreate");
         tskim = new TTree("skim", "omega->pi0,gamma event skim");
         tskim->Branch("runevent", &tskimrow.runevent, "runevent/D");
         tskim->Branch("Ebeam", &tskimrow.Ebeam, "Ebeam/D");
         tskim->Branch("Emiss", &tskimrow.Emiss, "Emiss/D");
         tskim->Branch("m3g", &tskimrow.m3g, "m3g/D");
         tskim->Branch("dphi", &tskimrow.dphi, "dphi/D");
         tskim->Branch("tabs", &tskimrow.tabs, "tabs/D");
         tskim->Branch("costhetahel", &tskimrow.costhetahel, "costhetahel/D");
         tskim->Branch("phihel", &tskimrow.phihel, "phihel/D");
         tskim->Branch("weight", &tskimrow.weight, "weight/D");
      }
      if (tagw != 0) {
         tskimrow.runevent = runno + eventno/1.0e7;
         tskimrow.Ebeam = coenergy[i];
         tskimrow.Emiss = Emiss;
         tskimrow.m3g = m3g;
         tskimrow.dphi = dphi;
         tskimrow.tabs = tabs;
         tskimrow.costhetahel = ppi0.CosTheta();
         tskimrow.phihel = ppi0.Phi();
         tskimrow.weight = tagw;
         tskim->Fill();
      }
#endif
      if (fabs(m3g - 0.8) < 0.1 &&
          Emiss < 0.5 && 
          fabs(dphi) < 50)
      {
         h2heli[nhist]->Fill(ppi0.CosTheta(), ppi0.Phi(), tagw);
      }
   }
   return kTRUE;
}

TLorentzVector ntuple_real::GetPpi0Helicity(double Ebeam)
{
   // Boost momenta into the helicty frame and there find the 4-momentum
   // of the pi0 assuming the 3gamma final state pi0(2gamma),gamma. To
   // To solve the kinematics of the reaction, I need to know the energy
   // of the beam Ebeam [GeV].

   struct p4set {
      TLorentzVector gam0;
      TLorentzVector gam1;
      TLorentzVector gam2;
      TLorentzVector fwd;
      TLorentzVector pi0;
      TLorentzVector bea;
      TLorentzVector rec;
      TLorentzVector ini;
      TLorentzVector tot;
   };
   const int p4set_size = 9;
   TLorentzVector *p4;
   int p4i;

   struct p4set p4lab;
   TLorentzVector *pgam = &p4lab.gam0;
   for (int n=0; n < 3; ++n) {
      pgam[n].SetPxPyPzE(pvect[n][1], pvect[n][2], pvect[n][3], pvect[n][0]);
   }

   // Find the pi0 pair, and its lab momentum.
   int nsolo;
   double mpi0 = 0.2;
   p4lab.fwd = pgam[0] + pgam[1] + pgam[2];
   for (int n=0; n < 3; ++n) {
      double m = (p4lab.fwd - pgam[n]).Mag();
      if (m < mpi0) {
         nsolo = n;
         mpi0 = m;
      }
   }
   p4lab.pi0 = p4lab.fwd - pgam[nsolo];

   // Solve for the final proton momentum, using energy conservation
   // and assuming that the entire energy transfer in the collision 
   // is absorbed on a single proton, which does not escapes with all
   // of this energy on the way out of the nucleus.
   double Mrecoil = 0.938;
   double nuclear_binding = 0.030; // average from MC
   double KErecoil = Ebeam - p4lab.fwd.E() - nuclear_binding;

   // There is an effective lower bound on the KE of the final proton
   // set by its range and the fact that it had to fire the trigger
   // scintillators, which I estimate to be 50 MeV. This depends on
   // angle, but the energy resolution does not justify trying to
   // be more accurate than this.
   KErecoil = (KErecoil > 0.05)? KErecoil : 0.05;

   double Erecoil = KErecoil + Mrecoil;
   double Precoil = sqrt(Erecoil*Erecoil - Mrecoil*Mrecoil);
   p4lab.rec.SetXYZM(Precoil * sin(therec[0]) * cos(phirec[0]),
                     Precoil * sin(therec[0]) * sin(phirec[0]),
                     Precoil * cos(therec[0]), 
                     Mrecoil);
 
   // Solve for the initial proton momentum, using momentum conservation
   p4lab.bea.SetXYZT(0, 0, Ebeam, Ebeam);
   p4lab.tot = p4lab.rec + p4lab.fwd;
   p4lab.ini = p4lab.tot - p4lab.bea;

   // Boost into the reaction rest frame
   struct p4set p4rea(p4lab);
   TVector3 CMbeta(p4lab.tot.BoostVector());
   for (p4 = &p4rea.gam0, p4i = 0; p4i < p4set_size; ++p4i, ++p4)
      p4->Boost(-CMbeta);

   // Put the beam back along the z axis in the reaction CM frame
   double phi0 = atan2(p4rea.bea.Py(), p4rea.bea.Px());
   for (p4 = &p4rea.gam0, p4i = 0; p4i < p4set_size; ++p4i, ++p4)
      p4->RotateZ(-phi0);
   double theta0 = atan2(p4rea.bea.Px(), p4rea.bea.Pz());
   for (p4 = &p4rea.gam0, p4i = 0; p4i < p4set_size; ++p4i, ++p4)
      p4->RotateY(-theta0);

   // Put the forward system momentum in the xz plane
   double phi1 = atan2(p4rea.fwd.Py(), p4rea.fwd.Px());
   for (p4 = &p4rea.gam0, p4i = 0; p4i < p4set_size; ++p4i, ++p4)
      p4->RotateZ(-phi1);

   // Rotate and boost into the s-channel helicity frame
   struct p4set p4hel(p4rea);
   double theta1 = atan2(p4rea.fwd.Px(), p4rea.fwd.Pz());
   for (p4 = &p4hel.gam0, p4i = 0; p4i < p4set_size; ++p4i, ++p4)
      p4->RotateY(-theta1);
   TVector3 helibeta(p4hel.fwd.BoostVector());
   for (p4 = &p4hel.gam0, p4i = 0; p4i < p4set_size; ++p4i, ++p4)
      p4->Boost(-helibeta);

#if DUMP_KINE_TREE
   // Save results in a tree for checkout offline
   if (tkine == 0) {
      if (dumpfile == 0)
         dumpfile = new TFile("ntuple_real_dump.root", "recreate");
      tkine = new TTree("kine", "3gamma event kinematics");
      tkine->Branch("Ebeam", &tkinerow.Ebeam, "Ebeam/D");
      tkine->Branch("plabfwd", tkinerow.plabfwd, "plabfwd[4]/D");
      tkine->Branch("plabpi0", tkinerow.plabpi0, "plabpi0[4]/D");
      tkine->Branch("plabbac", tkinerow.plabbac, "plabbac[4]/D");
      tkine->Branch("plabrec", tkinerow.plabrec, "plabrec[4]/D");
      tkine->Branch("plabini", tkinerow.plabini, "plabini[4]/D");
      tkine->Branch("plabtot", tkinerow.plabtot, "plabtot[4]/D");
      tkine->Branch("preafwd", tkinerow.preafwd, "preafwd[4]/D");
      tkine->Branch("preapi0", tkinerow.preapi0, "preapi0[4]/D");
      tkine->Branch("preabac", tkinerow.preabac, "preabac[4]/D");
      tkine->Branch("prearec", tkinerow.prearec, "prearec[4]/D");
      tkine->Branch("preaini", tkinerow.preaini, "preaini[4]/D");
      tkine->Branch("preatot", tkinerow.preatot, "preatot[4]/D");
      tkine->Branch("preabea", tkinerow.preabea, "preabea[4]/D");
      tkine->Branch("phelfwd", tkinerow.phelfwd, "phelfwd[4]/D");
      tkine->Branch("phelpi0", tkinerow.phelpi0, "phelpi0[4]/D");
      tkine->Branch("phelbac", tkinerow.phelbac, "phelbac[4]/D");
      tkine->Branch("phelrec", tkinerow.phelrec, "phelrec[4]/D");
      tkine->Branch("phelini", tkinerow.phelini, "phelini[4]/D");
      tkine->Branch("pheltot", tkinerow.pheltot, "pheltot[4]/D");
      tkine->Branch("phelbea", tkinerow.phelbea, "phelbea[4]/D");
      tkine->Branch("theta0", &tkinerow.theta0, "theta0/D");
      tkine->Branch("theta1", &tkinerow.theta1, "theta1/D");
      tkine->Branch("costhetahel", &tkinerow.costhetahel, "costhetahel/D");
      tkine->Branch("phihel", &tkinerow.phihel, "phihel/D");
   }
   tkinerow.Ebeam = Ebeam;
   p4lab.fwd.GetXYZT(tkinerow.plabfwd);
   p4lab.pi0.GetXYZT(tkinerow.plabpi0);
   (p4 = &p4lab.gam0)[nsolo].GetXYZT(tkinerow.plabbac);
   p4lab.rec.GetXYZT(tkinerow.plabrec);
   p4lab.ini.GetXYZT(tkinerow.plabini);
   p4lab.tot.GetXYZT(tkinerow.plabtot);
   p4rea.fwd.GetXYZT(tkinerow.preafwd);
   p4rea.pi0.GetXYZT(tkinerow.preapi0);
   (p4 = &p4rea.gam0)[nsolo].GetXYZT(tkinerow.preabac);
   p4rea.rec.GetXYZT(tkinerow.prearec);
   p4rea.ini.GetXYZT(tkinerow.preaini);
   p4rea.tot.GetXYZT(tkinerow.preatot);
   p4rea.bea.GetXYZT(tkinerow.preabea);
   p4hel.fwd.GetXYZT(tkinerow.phelfwd);
   p4hel.pi0.GetXYZT(tkinerow.phelpi0);
   (p4 = &p4hel.gam0)[nsolo].GetXYZT(tkinerow.phelbac);
   p4hel.rec.GetXYZT(tkinerow.phelrec);
   p4hel.ini.GetXYZT(tkinerow.phelini);
   p4hel.tot.GetXYZT(tkinerow.pheltot);
   p4hel.bea.GetXYZT(tkinerow.phelbea);
   tkinerow.theta0 = theta0;
   tkinerow.theta1 = theta1;
   tkinerow.costhetahel = p4hel.pi0.CosTheta();
   tkinerow.phihel = p4hel.pi0.Phi();
   tkine->Fill();
#endif

   return p4hel.pi0;
}

void ntuple_real::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

   if (dumpfile) {
      dumpfile->cd();
      tkine->Write();
      dumpfile->Close();
   }
   if (skimfile) {
      skimfile->cd();
      tskim->Write();
      skimfile->Close();
   }
}

void ntuple_real::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   TFile fout("ntuple_real.root", "recreate");
   TIter next(GetOutputList());
   TObject *obj;
   while ((obj = next()))
      obj->Write();
}

//////////////////////////////////////////////////////////
//
// mlfitter - A helper class for doing maximum likelihood
//            unbinned fit of experimental data describing
//            2-body decays from a vector to pi0, gamma.
//
// author: richard.t.jones at uconn.edu
// version: first written november 14, 2017
//
// usage: (a typical session, yours may be different)
// root [1] .L mlfitter.C++
// root [2] .L xyplots.C++
// root [3] mlfitter fit;
// root [4] fit.LoadAccept();
// root [5] fit.LoadNormInt();
// root [6] fit.LoadRealData();
// root [7] fit.Do("migrad");
//   ... repeat until happy with the fit ...
// root [n] xyplots(fit);
// root [n] hdatax->Draw();
// root [n] hfitx->Draw("hist same c");
// root [n] c1->Print("fitcheck1.png");
// root [n] .q
//
//////////////////////////////////////////////////////////

#include "mlfitter.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <limits>

TString mcsim_templ("/home/radphi/omega_SDME/mctuple_sim_skim-%d.root");
TString mcgen_templ("/home/radphi/omega_SDME/mctuple_gen_skim-%d.root");
TString expdata_templ("/home/radphi/omega_SDME/ntuple_real_skim-%2.2d.root");

const int nbins_cost = 50;
const int nbins_phi = 50;

// For multi-threaded operation this needs to be thread_local, but
// it causes a crash in Root 6.08 so until this is fixed, comment it.
/*thread_local*/ mlfitter *mlfitter::fFitter;

mlfitter::mlfitter()
{
   // set limits to default values

   SetLimits_Ebeam(0., 100.0);
   SetLimits_Emiss(-0.2, 0.5);
   SetLimits_tabs(-100., 0.5);
   SetLimits_dphi(-60., 60., -60., 60.);
   SetLimits_mass(0.70, 0.90, 0.60, 1.0);
   fRandom = new TRandom3(0);
}

mlfitter::~mlfitter()
{}

void mlfitter::LoadAccept()
{
   // Form an acceptance histogram from simulated events.
 
   TString haccept;
   if (fHaccept == 0) {
      for (int n=0; n < 99999; ++n) {
         haccept.Form("haccept%d", n);
         if (gROOT->FindObject(haccept) == 0)
            break;
      }
      fHaccept = new TH2D(haccept, "MC acceptance phi vs cos(theta)",
                          nbins_cost, -1, 1, nbins_phi, -M_PI, M_PI);
      fHaccept->SetStats(0);
      fHaccept->Sumw2();
   }
   fHaccept->Reset();

   // sum reconstructed MC for these cuts
   fChain = ChainMCsim(mcsim_templ);
   fReader.SetTree(fChain);
   while (fReader.Next()) {
      if ( *Ebeam > fEbeamlim[0] && *Ebeam < fEbeamlim[1] &&
           *Emiss > fEmisslim[0] && *Emiss < fEmisslim[1] &&
            *tabs > ftabslim[0]  &&  *tabs < ftabslim[1] )
      {
         double wgt_dphi, wgt_m3g;
         if (*dphi > fdphilim[0]  &&  *dphi < fdphilim[1])
            wgt_dphi = 1;
         else if (*dphi > fdphiban[0]  &&  *dphi < fdphiban[1])
            wgt_dphi = -fdphibanwgt;
         else
            wgt_dphi = 0;
         if (*m3g > fmasslim[0]  &&   *m3g < fmasslim[1])
            wgt_m3g = 1;
         else if (*m3g > fmassban[0]  &&   *m3g < fmassban[1])
            wgt_m3g = -fmassbanwgt;
         else
            wgt_m3g = 0;
         fHaccept->Fill(*costhetahel, *phihel, *weight * wgt_dphi * wgt_m3g);
      }
   }

   // normalize to the number generated
   fChain = ChainMCgen(mcgen_templ);
   fReader.SetTree(fChain);
   fNgenerated = 0;
   double nchecked = 0;
   while (fReader.Next()) {
      if ( *Ebeam > fEbeamlim[0] && *Ebeam < fEbeamlim[1] &&
           *Emiss > fEmisslim[0] && *Emiss < fEmisslim[1] &&
            *tabs > ftabslim[0]  &&  *tabs < ftabslim[1] )
      {
         double wgt_dphi, wgt_m3g;
         if (*dphi > fdphilim[0]  &&  *dphi < fdphilim[1])
            wgt_dphi = 1;
         else if (*dphi > fdphiban[0]  &&  *dphi < fdphiban[1])
            wgt_dphi = -fdphibanwgt;
         else
            wgt_dphi = 0;
         if (*m3g > fmasslim[0]  &&   *m3g < fmasslim[1])
            wgt_m3g = 1;
         else if (*m3g > fmassban[0]  &&   *m3g < fmassban[1])
            wgt_m3g = -fmassbanwgt;
         else
            wgt_m3g = 0;
         fNgenerated += *weight * wgt_dphi * wgt_m3g;
         // There are excessive statistics in this chain, no need
         // to count beyond 1M to know the normalization to 0.1%.
         if (fNgenerated > 1000000) {
            fNgenerated *= fChain->GetEntries() / nchecked;
            break;
         }
      }
      nchecked += 1;
   }
}

void mlfitter::LoadNormInt()
{
   // Evaluate the normalization integrals from input trees.

   if (fRfunction == 0)
      ClearRfunction();
   if (fNgenerated == 0)
      LoadAccept();

   // compute the normalization sums
   fNormInt11 = 0;
   fNormInt00 = 0;
   fNormInt1m = 0;
   fNormInt10 = 0;
   fChain = ChainMCsim(mcsim_templ);
   fReader.SetTree(fChain);
   double nsimulated = 0;
   double nRsimulated = 0;
   while (fReader.Next()) {
      if ( *Ebeam > fEbeamlim[0] && *Ebeam < fEbeamlim[1] &&
           *Emiss > fEmisslim[0] && *Emiss < fEmisslim[1] &&
            *tabs > ftabslim[0]  &&  *tabs < ftabslim[1] )
      {
         double wgt_dphi, wgt_m3g;
         if (*dphi > fdphilim[0]  &&  *dphi < fdphilim[1])
            wgt_dphi = 1;
         else if (*dphi > fdphiban[0]  &&  *dphi < fdphiban[1])
            wgt_dphi = -fdphibanwgt;
         else
            wgt_dphi = 0;
         if (*m3g > fmasslim[0]  &&   *m3g < fmasslim[1])
            wgt_m3g = 1;
         else if (*m3g > fmassban[0]  &&   *m3g < fmassban[1])
            wgt_m3g = -fmassbanwgt;
         else
            wgt_m3g = 0;
         int ibin = fRfunction->FindBin(*costhetahel, *phihel);
         double Rvalue = fRfunction->GetBinContent(ibin);
         double wR = *weight * wgt_dphi * wgt_m3g * Rvalue;
         fNormInt11 += W11(*costhetahel, *phihel) * wR / fNgenerated;
         fNormInt00 += W00(*costhetahel, *phihel) * wR / fNgenerated;
         fNormInt1m += W1m(*costhetahel, *phihel) * wR / fNgenerated;
         fNormInt10 += W10(*costhetahel, *phihel) * wR / fNgenerated;
         nsimulated += *weight;
         nRsimulated += wR;
      }
   }
   std::cout << "average acceptance is " 
             << nsimulated / fNgenerated
             << ", with effective statistics of "
             << rint(nRsimulated) << " events."
             << std::endl;
}

void mlfitter::ApplyCosthetaCut(double costlim0, double costlim1)
{
   // Set the value of the R function to zero in any bins on the
   // costhetahel,phihel plane outside the range (costlim0, costlim1).
   // After the R function is changed, LoadNormInt() and LoadRealData()
   // need to be invoked again.

   if (fRfunction == 0)
      ClearRfunction();

   // Zero bins in the R function outside the allowed region.
   for (int i=1; i <= nbins_cost; ++i) {
      double cost = fRfunction->GetXaxis()->GetBinCenter(i);
      if (cost < costlim0 || cost > costlim1)
         for (int j=1; j <= nbins_phi; ++j)
            fRfunction->SetBinContent(i, j, 0);
   }
}

void mlfitter::ApplyAcceptanceThreshold(double acclim0)
{
   // Set the value of the R function to zero in any bins on the
   // costhetahel,phihel plane where the acceptance drops below
   // some threshold value acclim0. After the R function is changed
   // LoadNormInt() and LoadRealData() need to be invoked again.

   if (fRfunction == 0)
      ClearRfunction();
   if (fHaccept == 0)
      LoadAccept();

   // Zero bins in the R function where acceptance is too low
   double accept_max = fHaccept->GetMaximum();
   for (int i=1; i <= nbins_cost; ++i) {
      for (int j=1; j <= nbins_phi; ++j) {
         double acc = fHaccept->GetBinContent(i, j);
         if (acc < acclim0 * accept_max)
            fRfunction->SetBinContent(i, j, 0);
      }
   }
}

void mlfitter::LoadRealData(int nevents, int nstart)
{
   // Load the real decay kinematics from an input ROOT tree
   // and load them in memory for looping over in GetW. Arguments
   // nevents, nstart select a subset of events from the dataset,
   // default values nevents=0, nstart=0 taking the whole thing.

   // check the R function
   if (fRfunction == 0) {
      std::cerr << "Error in LoadRealData: you must load Monte Carlo data"
                << " with LoadNormInt() before loading the real data,"
                << " without changing the R function in between!"
                << std::endl;
      return;
   }

   // create weight histograms if not already existing
   if (fHavewgt == 0) {
      TString hname;
      int n;
      for (n=0; n < 99999; ++n) {
         hname.Form("hwgt1_%d", n);
         if (gROOT->FindObject(hname) == 0)
            break;
      }
      fHavewgt = new TH2D(hname, "average weight on phi vs cos(theta)",
                          nbins_cost, -1, 1, nbins_phi, -M_PI, M_PI);
      fHavewgt->SetStats(0);
      fHavewgt->Sumw2();
      hname.Form("hwgt2_%d", n);
      fHavewgt2 = new TH2D(hname, "average weight**2 on phi vs cos(theta)",
                           nbins_cost, -1, 1, nbins_phi, -M_PI, M_PI);
      fHavewgt2->SetStats(0);
      fHavewgt2->Sumw2();
      hname.Form("hm3g_%d", n);
      fHm3g = new TH1D(hname, "3gamma mass spectrum", 150, 0, 1.5);
      fHm3g->Sumw2();
      hname.Form("hdphi_%d", n);
      fHdphi = new TH1D(hname, "delta phi", 90, -180., 180.);
      fHdphi->Sumw2();
   }
   fHavewgt->Reset();
   fHavewgt2->Reset();
   fHm3g->Reset();
   fHdphi->Reset();

   // now load the real data
   fChain = ChainRealData(expdata_templ);
   fSample.clear();
   fReader.SetTree(fChain);
   double sum_weights = 0;
   long int nrows = 0;
   while (fReader.Next()) {
      if (nrows++ < nstart)
         continue;
      else if (nevents > 0 && nrows > nstart + nevents)
         break;
      int ibin = fRfunction->FindBin(*costhetahel, *phihel);
      double Rvalue = fRfunction->GetBinContent(ibin);
      if ( *Ebeam > fEbeamlim[0] && *Ebeam < fEbeamlim[1] &&
           *Emiss > fEmisslim[0] && *Emiss < fEmisslim[1] &&
            *tabs > ftabslim[0]  &&  *tabs < ftabslim[1] )
      {
         double wgt_dphi, wgt_m3g;
         if (*dphi > fdphilim[0]  &&  *dphi < fdphilim[1])
            wgt_dphi = 1;
         else if (*dphi > fdphiban[0]  &&  *dphi < fdphiban[1])
            wgt_dphi = -fdphibanwgt;
         else
            wgt_dphi = 0;
         if (*m3g > fmasslim[0]  &&   *m3g < fmasslim[1])
            wgt_m3g = 1;
         else if (*m3g > fmassban[0]  &&   *m3g < fmassban[1])
            wgt_m3g = -fmassbanwgt;
         else
            wgt_m3g = 0;
         struct event_t datarow;
         datarow.W11 = W11(*costhetahel, *phihel);
         datarow.W00 = W00(*costhetahel, *phihel);
         datarow.W1m = W1m(*costhetahel, *phihel);
         datarow.W10 = W10(*costhetahel, *phihel);
         datarow.wgt = *weight * wgt_dphi * wgt_m3g;
         datarow.R = Rvalue;
         fSample.push_back(datarow);
         sum_weights += datarow.wgt * datarow.R;
         fHavewgt->Fill(*costhetahel, *phihel, datarow.wgt);
         fHavewgt2->Fill(*costhetahel, *phihel, pow(datarow.wgt, 2));
         if (*dphi > fdphilim[0]  &&  *dphi < fdphilim[1])
            fHm3g->Fill(*m3g, datarow.wgt * datarow.R);
         if (*m3g > fmasslim[0]  &&   *m3g < fmasslim[1])
            fHdphi->Fill(*dphi, datarow.wgt * datarow.R);
      }
   }
   std::cout << "total events loaded is " << fSample.size()
             << " out of " << nrows << " events scanned,"
             << " weighted sum is " << sum_weights
             << std::endl;
}

void mlfitter::LoadRfunction()
{
   // Compute the R function from Havewgt and Havewgt2
   // assuming that both of these exist, otherwise print
   // an error and exit. After the R function is changed
   // LoadNormInt() and LoadRealData() need to be invoked again.

   if (fHavewgt == 0 || fHavewgt2 == 0 || fRfunction == 0) {
      std::cerr << "Error in LoadRfunction: you must load real data"
                << " with LoadRealData() before loading the R function."
                << std::endl;
      return;
   }

   for (int i=1; i <= nbins_cost; ++i) {
      for (int j=1; j <= nbins_phi; ++j) {
         double wgt = fHavewgt->GetBinContent(i,j);
         double wgt2 = fHavewgt2->GetBinContent(i,j) + 1e-99;
         if (wgt2 > 0)
            fRfunction->SetBinContent(i, j, wgt / wgt2);
         else
            fRfunction->SetBinContent(i, j, 0);
      }
   }
}

void mlfitter::ClearRfunction()
{
   // If Rfunction histogram does not exist, create a new one;
   // either way, set the contents of all bins to unity. After
   // the R function is changed, LoadNormInt() and LoadRealData()
   // need to be invoked again.

   if (fRfunction == 0) {
      TString hrfunc;
      for (int n=0; n < 99999; ++n) {
         hrfunc.Form("hRfunc%d", n);
         if (gROOT->FindObject(hrfunc) == 0)
            break;
      }
      fRfunction = new TH2D(hrfunc, "R function on phi vs cos(theta)",
                            nbins_cost, -1, 1, nbins_phi, -M_PI, M_PI);
      fRfunction->SetStats(0);
   }
   fRfunction->Reset();

   for (int i=1; i <= nbins_cost; ++i)
      for (int j=1; j <= nbins_phi; ++j)
         fRfunction->SetBinContent(i, j, 1);
}

TMinuit *mlfitter::GetMinuit()
{
   // Initialize the Minuit search engine and attach it to this object,
   // or return the existing one, if it already exists.

   if (fMinuit)
      return fMinuit;

   fMinuit = new TMinuit(4);
   fMinuit->SetFCN(FCN);

   // define the parameters and starting values
   for (int i=0; i < 4; ++i) {
      fLastParam[i] = fRandom->Rndm();
   }
   int ierflg = 0;
   double stepsize = 1.e-4;
   fMinuit->mnparm(0, "rho11", fLastParam[0], stepsize, 0, 0, ierflg);
   fMinuit->mnparm(1, "rho00", fLastParam[1], stepsize, 0, 0, ierflg);
   fMinuit->mnparm(2, "rho1M", fLastParam[2], stepsize, 0, 0, ierflg);
   fMinuit->mnparm(3, "rho10", fLastParam[3], stepsize, 0, 0, ierflg);

   // set default values for other parameters
   Do("SET PRINTOUT", 1);
   Do("SET ERRORDEF", 0.5);
   Do("SET WARNINGS");
   Do("SET GRADIENT", 0);
   Do("SET STRATEGY", 2);
   return fMinuit;
}

TChain *mlfitter::ChainMCsim(const TString templ)
{
   // Construct and return a TChain containing all of the ROOT trees
   // that constitute the MC simulation sample for this analysis.

   TChain *ch = new TChain("skim");
   int n = 0;
   while (++n) {
      TString fpath;
      fpath.Form(templ, n);
      std::ifstream froot(fpath);
      if (froot)
         ch->Add(fpath);
      else
         break;
   }
   return ch;
}

TChain *mlfitter::ChainMCgen(const TString templ)
{
   // Construct and return a TChain containing all of the ROOT trees
   // that constitute the MC generated sample for this analysis.

   TChain *ch = new TChain("skimmc");
   int n = 0;
   while (++n) {
      TString fpath;
      fpath.Form(templ, n);
      std::ifstream froot(fpath);
      if (froot)
         ch->Add(fpath);
      else
         break;
   }
   return ch;
}

TChain *mlfitter::ChainRealData(const TString templ)
{
   // Construct and return a TChain containing all of the ROOT trees
   // that constitute the real data sample for this analysis.

   TChain *ch = new TChain("skim");
   int n = 0;
   while (++n) {
      TString fpath;
      fpath.Form(templ, n);
      std::ifstream froot(fpath);
      if (froot)
         ch->Add(fpath);
      else
         break;
   }
   return ch;
}

double mlfitter::W11(double costheta, double phi) const
{
   // The decay angular distribution that goes with rho11.

   return 3 / (8 * M_PI) * (1 + costheta * costheta);
}

double mlfitter::W00(double costheta, double phi) const
{
   // The decay angular distribution that goes with rho00.
 
   return 3 / (8 * M_PI) * (1 - costheta * costheta);
}

double mlfitter::W1m(double costheta, double phi) const
{
   // The decay angular distribution that goes with rho1M.

   return 3 / (8 * M_PI) * (1 - costheta * costheta) * cos(2 * phi);
}

double mlfitter::W10(double costheta, double phi) const
{
   // The decay angular distribution that goes with rho10.

   return 3 / (8 * M_PI) * costheta * 
          sqrt(8 * (1 - costheta * costheta)) * cos(phi);
}

double mlfitter::GetFCN(double rho11, double rho00, double rho1M, double rho10)
{
   // This is the core method of the mlfitter class. It computes
   // the likelihood of the data for the given model parameters
   // and Monte Carlo acceptance. It also computes the gradient,
   // and the second derivative matrix as a by-product.

   fLastFCN = 0;
   fLastParam[0] = rho11;
   fLastParam[1] = rho00;
   fLastParam[2] = rho1M;
   fLastParam[3] = rho10;
   for (int i=0; i < 4; i++) {
      fLastGradient[i] = 0;
      for (int j=0; j < 4; j++) {
         fLastHesse[i][j] = 0;
      }
   }

   // Do the log sum over real events
   long int nrow = 0;
   std::vector<struct event_t>::const_iterator iter;
   for (iter = fSample.begin(); iter != fSample.end(); ++iter, ++nrow) {
      if (fParts > 1 && ((fPart > 0 && nrow % fParts != fPart - 1) ||
                         (fPart < 0 && nrow % fParts == -fPart - 1)))
      {
         continue;
      }
      
      double a = rho11 * iter->W11 + rho00 * iter->W00 +
                 rho1M * iter->W1m + rho10 * iter->W10 + 1e-99;
      // compress negative values for a into the range [0,fWlowcut]
      double dlogada = 1 / a;
      double d2logada2 = -1 / (a*a);
      if (fWlowcut > 0 && a < fWlowcut) {
         a = fWlowcut * exp(a - fWlowcut) + 1e-99;
         dlogada = 1;
         d2logada2 = 0;
      }
      double wR = iter->wgt * iter->R;
      fLastFCN -= wR * log(a);
      fLastGradient[0] -= wR * iter->W11 * dlogada;
      fLastGradient[1] -= wR * iter->W00 * dlogada;
      fLastGradient[2] -= wR * iter->W1m * dlogada;
      fLastGradient[3] -= wR * iter->W10 * dlogada;
      fLastHesse[0][0] -= wR * iter->W11 * iter->W11 * d2logada2;
      fLastHesse[0][1] -= wR * iter->W11 * iter->W00 * d2logada2;
      fLastHesse[0][2] -= wR * iter->W11 * iter->W1m * d2logada2;
      fLastHesse[0][3] -= wR * iter->W11 * iter->W10 * d2logada2;
      fLastHesse[1][1] -= wR * iter->W00 * iter->W00 * d2logada2;
      fLastHesse[1][2] -= wR * iter->W00 * iter->W1m * d2logada2;
      fLastHesse[1][3] -= wR * iter->W00 * iter->W10 * d2logada2;
      fLastHesse[2][2] -= wR * iter->W1m * iter->W1m * d2logada2;
      fLastHesse[2][3] -= wR * iter->W1m * iter->W10 * d2logada2;
      fLastHesse[3][3] -= wR * iter->W10 * iter->W10 * d2logada2;
   }

   // Apply the normalization
   fLastFCN += rho11 * fNormInt11 + rho00 * fNormInt00 +
               rho1M * fNormInt1m + rho10 * fNormInt10;
   fLastGradient[0] += fNormInt11;
   fLastGradient[1] += fNormInt00;
   fLastGradient[2] += fNormInt1m;
   fLastGradient[3] += fNormInt10;
 
   // Add a barrier factor for violating physical bounds on parameters
   // by checking that all rulecheck variables are non-negative.
   double rulecheck1 = rho11;
   double rulecheck2 = rho00;
   double rulecheck3 = rho11 * rho11 - rho1M * rho1M;
   double rulecheck4 = rho11 * rho00 - rho10 * rho10;
   if (fRule1 > 0) {
      double barrier = exp(-fRule1 * rulecheck1);
      fLastFCN += barrier;
      fLastGradient[0] += -fRule1 * barrier;
      fLastHesse[0][0] += fRule1 * fRule1 * barrier;
   }
   if (fRule2 > 0) {
      double barrier = exp(-fRule2 * rulecheck2);
      fLastFCN += barrier;
      fLastGradient[1] += -fRule2 * barrier;
      fLastHesse[1][1] += fRule2 * fRule2 * barrier;
   }
   if (fRule3 > 0) {
      double barrier = exp(-fRule3 * rulecheck3);
      fLastFCN += barrier;
      fLastGradient[0] += -fRule3 * barrier * (+2 * rho11);
      fLastGradient[2] += -fRule3 * barrier * (-2 * rho1M);
      fLastHesse[0][0] += fRule3 * fRule3 * barrier * (4 * rho11 * rho11) -
                          fRule3 * barrier * 2;
      fLastHesse[0][2] += fRule3 * fRule3 * barrier * (-4 * rho11 * rho1M);
      fLastHesse[2][2] += fRule3 * fRule3 * barrier * (4 * rho1M * rho1M) +
                          fRule3 * barrier * 2;
   }
   if (fRule4 > 0) {
      double barrier = exp(-fRule4 * rulecheck4);
      fLastFCN += barrier;
      fLastGradient[0] += -fRule4 * barrier * (rho00);
      fLastGradient[1] += -fRule4 * barrier * (rho11);
      fLastGradient[3] += -fRule4 * barrier * (-2 * rho10);
      fLastHesse[0][0] += fRule4 * fRule4 * barrier * (rho00 * rho00);
      fLastHesse[0][1] += fRule4 * fRule4 * barrier * (rho00 * rho11) -
                          fRule4 * barrier;
      fLastHesse[0][3] += fRule4 * fRule4 * barrier * (-2 * rho00 * rho10);
      fLastHesse[1][1] += fRule4 * fRule4 * barrier * (rho11 * rho11);
      fLastHesse[1][3] += fRule4 * fRule4 * barrier * (-2 * rho11 * rho10);
      fLastHesse[3][3] += fRule4 * fRule4 * barrier * (4 * rho10 * rho10) +
                          fRule4 * barrier * 2;
   }

   // Fill out the lower triangular part of the Hesse matrix
   fLastHesse[1][0] = fLastHesse[0][1];
   fLastHesse[2][0] = fLastHesse[0][2];
   fLastHesse[3][0] = fLastHesse[0][3];
   fLastHesse[2][1] = fLastHesse[1][2];
   fLastHesse[3][1] = fLastHesse[1][3];
   fLastHesse[3][2] = fLastHesse[2][3];
   return fLastFCN;
}

void mlfitter::FCN(int &npar, double *gin, double &f, double *par, int iflag)
{
   // This is a static function, so it does not have direct
   // access to any of the member variables of the class.
   // This gets called directly by the Minuit fit engine
   // once per step in the minimization search algorithm.

   if (iflag == 1) {                     // initialization
      return;
   }
   else if (iflag == 3) {                // termination
      return;
   }

   f = fFitter->GetFCN(par[0], par[1], par[2], par[3]);

   if (iflag == 2) {                // compute FCN with derivatives
      gin[0] = fFitter->fLastGradient[0];
      gin[1] = fFitter->fLastGradient[1];
      gin[2] = fFitter->fLastGradient[2];
      gin[3] = fFitter->fLastGradient[3];
   }
}

TH2D *mlfitter::GenerateDataPlot(const TString &hname)
{
   // Construct and return a 2D histogram of the decay
   // angular distribution of decays from real data.

   TH2D *h2 = (TH2D*)gROOT->FindObject(hname);
   if (h2 == 0) {
      h2 = new TH2D(hname, "measured omega decay phi vs cos(theta)",
                    nbins_cost, -1, 1, nbins_phi, -M_PI, M_PI);
      h2->Sumw2();
   }
   h2->Reset();
   fChain = ChainRealData(expdata_templ);
   fReader.SetTree(fChain);
   while (fReader.Next()) {
      if ( *Ebeam > fEbeamlim[0] && *Ebeam < fEbeamlim[1] &&
           *Emiss > fEmisslim[0] && *Emiss < fEmisslim[1] &&
            *tabs > ftabslim[0]  &&  *tabs < ftabslim[1] )
      {
         double wgt_dphi, wgt_m3g;
         if (*dphi > fdphilim[0]  &&  *dphi < fdphilim[1])
            wgt_dphi = 1;
         else if (*dphi > fdphiban[0]  &&  *dphi < fdphiban[1])
            wgt_dphi = -fdphibanwgt;
         else
            wgt_dphi = 0;
         if (*m3g > fmasslim[0]  &&   *m3g < fmasslim[1])
            wgt_m3g = 1;
         else if (*m3g > fmassban[0]  &&   *m3g < fmassban[1])
            wgt_m3g = -fmassbanwgt;
         else
            wgt_m3g = 0;
         h2->Fill(*costhetahel, *phihel, *weight * wgt_dphi * wgt_m3g);
      }
   }
   return h2;
}

TH2D *mlfitter::GenerateFitModel(const TString &hname)
{
   // Construct and return a 2D histogram of the decay
   // angular distribution of decays from the fit model,
   // including acceptance from Monte Carlo. The model
   // parameter values are taken from fLastParam.

   TH2D *h2 = (TH2D*)gROOT->FindObject(hname);
   if (h2 == 0) {
      h2 = new TH2D(hname, "model omega decay phi vs cos(theta)",
                    nbins_cost, -1, 1, nbins_phi, -M_PI, M_PI);
      h2->Sumw2();
   }
   h2->Reset();

   fChain = ChainMCsim(mcsim_templ);
   fReader.SetTree(fChain);
   while (fReader.Next()) {
      if ( *Ebeam > fEbeamlim[0] && *Ebeam < fEbeamlim[1] &&
           *Emiss > fEmisslim[0] && *Emiss < fEmisslim[1] &&
            *tabs > ftabslim[0]  &&  *tabs < ftabslim[1] )
      {
         double wgt_dphi, wgt_m3g;
         if (*dphi > fdphilim[0]  &&  *dphi < fdphilim[1])
            wgt_dphi = 1;
         else if (*dphi > fdphiban[0]  &&  *dphi < fdphiban[1])
            wgt_dphi = -fdphibanwgt;
         else
            wgt_dphi = 0;
         if (*m3g > fmasslim[0]  &&   *m3g < fmasslim[1])
            wgt_m3g = 1;
         else if (*m3g > fmassban[0]  &&   *m3g < fmassban[1])
            wgt_m3g = -fmassbanwgt;
         else
            wgt_m3g = 0;
         double wgt_prod = *weight * wgt_dphi * wgt_m3g;
         int ibin = fRfunction->FindBin(*costhetahel, *phihel);
         double W = GetW(*costhetahel, *phihel);
         h2->Fill(*costhetahel, *phihel, wgt_prod * W / fNgenerated);
      }
   }
   return h2;
}

void mlfitter::PartitionEvents(int part, int parts)
{
   // Tell the fitter which subset of events in the real dataset
   // should be included in the fits. The events are divided into
   // "parts" equal subsets. If part is in the range [1,parts]
   // then events in subset "part" are included. If part is negative
   // then all events are included EXCEPT those in partition |part|.
   // Values of |part| less than 1 or greater than parts are illegal.

   if (parts == 0 or abs(part) == 0 or abs(part) > parts) {
      std::cerr << "Error in PartitionEvents - combination"
                << " part=" << part << " parts=" << parts
                << " is not supported, command ignored!"
                << std::endl;
      return;
   }
   fParts = parts;
   fPart = part;
}

void mlfitter::RandomizeParams(double scale)
{
   // Generate random values for the model parameters
   // and call GetFCN to compute the model distribution
   // for these random parameters. The scale for the
   // parameter size is given by the last value for that
   // parameter, or unity if the last value is zero.

   double par[4];
   for (int i=0; i < 4; ++i) {
      par[i] = fRandom->Rndm();
      if (scale > 0) {
         if (i < 2)
            par[i] = par[i] * scale;
         else
            par[i] = (par[i] - 0.5) * scale/10.;
      }
      Do("SET PAR", i+1 , par[i]);
   }
   GetFCN(par[0], par[1], par[2], par[3]);
   Do("SHOW FCN");
}

int mlfitter::Fit(int redundancy, int maxtries)
{
   // Find the maximum likelihood values for the model parameters,
   // starting from the last parameter values seen by the fitter,
   // by seeking the minimum of the NLL function using the migrad
   // algorithm of Minuit. Migrad commands are issued one after
   // another until either (a) <redundancy> attempts are successful
   // and terminate at the same minimum, or (b) the number of tries
   // exceeds maxtries. Return value is 0 for case (a), and nonzero
   // for case (b). Upon success, the Minuit object returned by
   // GetMinuit can be queried for the final parameter values and
   // errors, the estimated error matrix, correlations, EDM, etc.

   double finalFCN(std::numeric_limits<double>::max());
   double tolerFCN = 0.01;
   int goodfinds = 0;
   for (int tries=0; tries < maxtries; ++tries) {
      RandomizeParams();
      int res = Do("migrad");
      std::cout << "============================"
                << "============================"
                << "============================"
                << " " << res << std::endl;
      if (res == 0) {
         Do("call 4");
         if (fLastFCN < finalFCN + tolerFCN) {
            if (++goodfinds == redundancy) {
               return 0;
            }
            else if (fLastFCN < finalFCN) {
               finalFCN = fLastFCN;
            }
         }
         std::cout << "goodfinds=" << goodfinds
                   << std::setprecision(15)
                   << ", finalFCN=" << finalFCN
                   << ", lastFCN=" << fLastFCN
                   << std::endl;
      }
   }
   return 1;
}

void mlfitter::Print() const
{
   // Print a summary of the fitter state to the screen

   std::cout << "lastFCN: " << fLastFCN << std::endl
             << "lastParam: " << fLastParam[0] << " "
                               << fLastParam[1] << " "
                               << fLastParam[2] << " "
                               << fLastParam[3] << std::endl
             << "lastGradient: " << fLastGradient[0] << " "
                               << fLastGradient[1] << " "
                               << fLastGradient[2] << " "
                               << fLastGradient[3] << std::endl
             << "lastHesse: " << fLastHesse[0][0] << " "
                               << fLastHesse[0][1] << " "
                               << fLastHesse[0][2] << " "
                               << fLastHesse[0][3] << std::endl
             << "            " << fLastHesse[1][0] << " "
                               << fLastHesse[1][1] << " "
                               << fLastHesse[1][2] << " "
                               << fLastHesse[1][3] << std::endl
             << "            " << fLastHesse[2][0] << " "
                               << fLastHesse[2][1] << " "
                               << fLastHesse[2][2] << " "
                               << fLastHesse[2][3] << std::endl
             << "            " << fLastHesse[3][0] << " "
                               << fLastHesse[3][1] << " "
                               << fLastHesse[3][2] << " "
                               << fLastHesse[3][3] << std::endl
             << "Ebeam_limits: " << fEbeamlim[0] << " "
                               << fEbeamlim[1] << std::endl
             << "Emiss_limits: " << fEmisslim[0] << " "
                               << fEmisslim[1] << std::endl
             << "tabs_limits: " << ftabslim[0] << " "
                               << ftabslim[1] << std::endl
             << "mass_limits: " << fmasslim[0] << " "
                               << fmasslim[1] << std::endl
             << "mass_sidebands: " << fmassban[0] << " "
                               << fmassban[1] << " "
                               << fmassbanwgt << std::endl
             << "dphi_limits: " << fdphilim[0] << " "
                               << fdphilim[1] << std::endl
             << "dphi_sidebands: " << fdphiban[0] << " "
                               << fdphiban[1] << " "
                               << fdphibanwgt << std::endl
             << std::endl;

   double sum_weights = 0;
   std::vector<struct event_t>::const_iterator iter;
   for (iter = fSample.begin(); iter != fSample.end(); ++iter) {
      sum_weights += iter->wgt * iter->R;
   } 

   std::cout << "real events loaded: " << fSample.size()
             << ", " <<  sum_weights << " weighted" << std::endl;
   std::cout << "selected partition: " << fPart
             << " out of " << fParts << std::endl;
   std::cout << "normalization integrals: " << fNormInt11 << " "
             << fNormInt00 << " " << fNormInt1m << " "
             << fNormInt10 << std::endl;
   std::cout << "normalization sum: " << fNormInt11 * fLastParam[0] +
                                         fNormInt00 * fLastParam[1] +
                                         fNormInt1m * fLastParam[2] +
                                         fNormInt10 * fLastParam[3]
             << std::endl;
   std::cout << "acceptance histogram: "
             << ((fHaccept)? fHaccept->GetName() : "") << std::endl;
   std::cout << "average wgt histogram: "
             << ((fHavewgt)? fHavewgt->GetName() : "") << std::endl;
   std::cout << "average wgt**2 histogram: "
             << ((fHavewgt2)? fHavewgt2->GetName() : "") << std::endl;
   std::cout << "R function histogram: "
             << ((fRfunction)? fRfunction->GetName() : "") << std::endl;
}

void mlfitter::SaveState(const char* filename) const
{
   // Save fitter state variables to a file,
   // complementary format to RestoreState.

   std::ofstream fsave(filename);
   if (fsave) {
      fsave << "lastFCN: " << fLastFCN << std::endl
            << "lastParam: " << fLastParam[0] << " "
                              << fLastParam[1] << " "
                              << fLastParam[2] << " "
                              << fLastParam[3] << std::endl
            << "lastGradient: " << fLastGradient[0] << " "
                              << fLastGradient[1] << " "
                              << fLastGradient[2] << " "
                              << fLastGradient[3] << std::endl
            << "lastHesse: " << fLastHesse[0][0] << " "
                              << fLastHesse[0][1] << " "
                              << fLastHesse[0][2] << " "
                              << fLastHesse[0][3] << std::endl
            << "            " << fLastHesse[1][0] << " "
                              << fLastHesse[1][1] << " "
                              << fLastHesse[1][2] << " "
                              << fLastHesse[1][3] << std::endl
            << "            " << fLastHesse[2][0] << " "
                              << fLastHesse[2][1] << " "
                              << fLastHesse[2][2] << " "
                              << fLastHesse[2][3] << std::endl
            << "            " << fLastHesse[3][0] << " "
                              << fLastHesse[3][1] << " "
                              << fLastHesse[3][2] << " "
                              << fLastHesse[3][3] << std::endl
            << "Ebeam_limits: " << fEbeamlim[0] << " "
                              << fEbeamlim[1] << std::endl
            << "Emiss_limits: " << fEmisslim[0] << " "
                              << fEmisslim[1] << std::endl
            << "tabs_limits: " << ftabslim[0] << " "
                              << ftabslim[1] << std::endl
            << "mass_limits: " << fmasslim[0] << " "
                              << fmasslim[1] << std::endl
            << "mass_sidebands: " << fmassban[0] << " "
                              << fmassban[1] << " "
                              << fmassbanwgt << std::endl
            << "dphi_limits: " << fdphilim[0] << " "
                              << fdphilim[1] << std::endl
            << "dphi_sidebands: " << fdphiban[0] << " "
                              << fdphiban[1] << " "
                              << fdphibanwgt << std::endl;
   }
   else {
      std::cerr << "Error in SaveState: could not open file "
                << filename << " for output."
                << std::endl;
   }
}

void mlfitter::RestoreState(const char* filename)
{
   // Load fitter state variables from a file
   // that was previously written by SaveState.

   std::ifstream fsave(filename);
   if (fsave) {
      std::string key;
      while (fsave >> key) {
         if (key == "lastFCN:")
            fsave >> fLastFCN;
         else if (key == "lastParam:")
            fsave >> fLastParam[0] >> fLastParam[1]
                  >> fLastParam[2] >> fLastParam[3];
         else if (key == "lastGradient:")
            fsave >> fLastGradient[0] >> fLastGradient[1]
                  >> fLastGradient[2] >> fLastGradient[3];
         else if (key == "lastHesse:")
            fsave >> fLastHesse[0][0] >> fLastHesse[0][1]
                  >> fLastHesse[0][2] >> fLastHesse[0][3]
                  >> fLastHesse[1][0] >> fLastHesse[1][1]
                  >> fLastHesse[1][2] >> fLastHesse[1][3]
                  >> fLastHesse[2][0] >> fLastHesse[2][1]
                  >> fLastHesse[2][2] >> fLastHesse[2][3]
                  >> fLastHesse[3][0] >> fLastHesse[3][1]
                  >> fLastHesse[3][2] >> fLastHesse[3][3];
         else if (key == "Ebeam_limits:")
            fsave >> fEbeamlim[0] >> fEbeamlim[1];
         else if (key == "Emiss_limits:")
            fsave >> fEmisslim[0] >> fEmisslim[1];
         else if (key == "tabs_limits:")
            fsave >> ftabslim[0] >> ftabslim[1];
         else if (key == "mass_limits:")
            fsave >> fmasslim[0] >> fmasslim[1];
         else if (key == "mass_sidebands:")
            fsave >> fmassban[0] >> fmassban[1] >> fmassbanwgt;
         else if (key == "dphi_limits:")
            fsave >> fdphilim[0] >> fdphilim[1];
         else if (key == "dphi_sidebands:")
            fsave >> fdphiban[0] >> fdphiban[1] >> fdphibanwgt;
      }
   }
   else {
      std::cerr << "Error in RestoreState: could not open file "
                << filename << " for input."
                << std::endl;
   }
}

/////////////////////////////////////////////////////////
//
// scans.C - Action routines for SDME analysis of the
//           decays of vector mesons with Radphi.
//
// author: richard.t.jones at uconn.edu
// version: first written november 20, 2017
//
//////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TH2D.h>
#include <TCanvas.h>

#include "mlfitter.h"
#include "xyplots.h"

void scanRcut(mlfitter &fit, 
              int nbins=15, 
              double cutlim0=0.0585,
              double cutlim1=0.0735,
              double cutcost0=-1.0,
              double cutcost1=0.8)
{
   // Do a scan over the selection threshold for the R function,
   // and plot the SDME values from the fit as a function of the
   // threshold. The desired kinematic cuts for the scan must be
   // already set in the mlfitter fit, but the data do not need
   // to be loaded.

   TH1D *hrhos[4];
   hrhos[0] = new TH1D("hrho11s", "rho11 vs hdcut selection threshold",
                       nbins, cutlim0, cutlim1);
   hrhos[1] = new TH1D("hrho00s", "rho00 vs hdcut selection threshold",
                       nbins, cutlim0, cutlim1);
   hrhos[2] = new TH1D("hrho1Ms", "rho1M vs hdcut selection threshold",
                       nbins, cutlim0, cutlim1);
   hrhos[3] = new TH1D("hrho10s", "rho10 vs hdcut selection threshold",
                       nbins, cutlim0, cutlim1);
   hrhos[0]->GetYaxis()->SetTitle("\\rho 11");
   hrhos[1]->GetYaxis()->SetTitle("\\rho 00");
   hrhos[2]->GetYaxis()->SetTitle("\\rho 1-1");
   hrhos[3]->GetYaxis()->SetTitle("\\rho 10");
   for (int n=0; n < 4; ++n) {
      hrhos[n]->GetXaxis()->SetTitle("ethresh");
      hrhos[n]->SetStats(0);
   }

   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1");

   for (int bin=0; bin < nbins; ++bin) {
      fit.SetLimits_mass(0.74, 0.86);
      fit.LoadAccept();
      fit.ClearRfunction();
      TH2D *hR = (TH2D*)gROOT->FindObject("hRfunc0");
      double ethresh = hrhos[0]->GetXaxis()->GetBinCenter(bin+1);
      TString title;
      title.Form("R function on phi vs cos(theta), cut value %f", ethresh);
      hR->SetTitle(title);
      if (fabs(ethresh) > 1e-6) {
         xyplots(fit);
         TH2D *hdcor = select2D(ethresh);
         mask2D(hR, hdcor);
         holefill2D(hR);
         fit.ApplyCosthetaCut(cutcost0, cutcost1);
      }
      hR->Draw("colz");
      c1->Update();
      fit.SetLimits_mass(0.70, 0.90, 0.58, 1.0);
      fit.LoadAccept();
      fit.LoadNormInt();
      fit.LoadRealData();
      if (fit.Fit() == 0) {
         for (int par=0; par < 4; ++par) {
            int npar;
            double parval;
            double parerr;
            fit.GetMinuit()->GetParameter(par, parval, parerr);
            hrhos[par]->SetBinContent(bin+1, parval);
            hrhos[par]->SetBinError(bin+1, parerr);
         }
      }
   }

   TH1D *hrhon[4];
   hrhon[0] = copy1D(hrhos[0], "hrho11n");
   hrhon[1] = copy1D(hrhos[1], "hrho00n");
   hrhon[2] = copy1D(hrhos[2], "hrho1Mn");
   hrhon[3] = copy1D(hrhos[3], "hrho10n");
   for (int n=1; n <= nbins; ++n) {
      double rho11 = hrhon[0]->GetBinContent(n);
      double rho00 = hrhon[1]->GetBinContent(n);
      double rho1M = hrhon[2]->GetBinContent(n);
      double rho10 = hrhon[3]->GetBinContent(n);
      double rho11e = hrhon[0]->GetBinError(n);
      double rho00e = hrhon[1]->GetBinError(n);
      double rho1Me = hrhon[2]->GetBinError(n);
      double rho10e = hrhon[3]->GetBinError(n);
      double trace = 2 * rho11 + rho00;
      hrhon[0]->SetBinContent(n, rho11 /trace);
      hrhon[1]->SetBinContent(n, rho00 /trace);
      hrhon[2]->SetBinContent(n, rho1M /trace);
      hrhon[3]->SetBinContent(n, rho10 /trace);
      hrhon[0]->SetBinError(n, rho11e /trace);
      hrhon[1]->SetBinError(n, rho00e /trace);
      hrhon[2]->SetBinError(n, rho1Me /trace);
      hrhon[3]->SetBinError(n, rho10e /trace);
   }

   for (int n=0; n < 4; ++n) {
      hrhos[n]->SetDirectory(0);
      hrhon[n]->SetDirectory(0);
   }
   TFile fout("scans.root", "update");
   for (int n=0; n < 4; ++n) {
      hrhos[n]->Write();
      hrhon[n]->Write();
   }
}

void scanMcut(mlfitter &fit, 
              int nbins=15, 
              double cutlim0=0.955,
              double cutlim1=1.555)
{
   // Do a scan over the mass subtraction weight, and plot
   // the SDME values from the fit as a function of the cut.
   // All of the other kinematic cuts for the scan must be
   // already set in the mlfitter fit, but the data do not
   // need to be loaded.

   TH1D *hrhos[4];
   hrhos[0] = new TH1D("hrho11s", "rho11 vs mass subtraction weight",
                       nbins, cutlim0, cutlim1);
   hrhos[1] = new TH1D("hrho00s", "rho00 vs hdcut subtraction weight",
                       nbins, cutlim0, cutlim1);
   hrhos[2] = new TH1D("hrho1Ms", "rho1M vs hdcut subtraction weight",
                       nbins, cutlim0, cutlim1);
   hrhos[3] = new TH1D("hrho10s", "rho10 vs hdcut subtraction weight",
                       nbins, cutlim0, cutlim1);
   hrhos[0]->GetYaxis()->SetTitle("\\rho 11");
   hrhos[1]->GetYaxis()->SetTitle("\\rho 00");
   hrhos[2]->GetYaxis()->SetTitle("\\rho 1-1");
   hrhos[3]->GetYaxis()->SetTitle("\\rho 10");
   for (int n=0; n < 4; ++n) {
      hrhos[n]->GetXaxis()->SetTitle("wgt");
      hrhos[n]->SetStats(0);
   }

   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1");

   for (int bin=0; bin < nbins; ++bin) {
      double subw = hrhos[0]->GetXaxis()->GetBinCenter(bin+1);
      fit.SetLimits_mass(0.70, 0.90, 0.6, 1.0, subw);
      fit.LoadAccept();
      fit.LoadNormInt();
      fit.LoadRealData();
      if (fit.Fit() == 0) {
         for (int par=0; par < 4; ++par) {
            int npar;
            double parval;
            double parerr;
            fit.GetMinuit()->GetParameter(par, parval, parerr);
            hrhos[par]->SetBinContent(bin+1, parval);
            hrhos[par]->SetBinError(bin+1, parerr);
         }
      }
   }

   TH1D *hrhon[4];
   hrhon[0] = copy1D(hrhos[0], "hrho11n");
   hrhon[1] = copy1D(hrhos[1], "hrho00n");
   hrhon[2] = copy1D(hrhos[2], "hrho1Mn");
   hrhon[3] = copy1D(hrhos[3], "hrho10n");
   for (int n=1; n <= nbins; ++n) {
      double rho11 = hrhon[0]->GetBinContent(n);
      double rho00 = hrhon[1]->GetBinContent(n);
      double rho1M = hrhon[2]->GetBinContent(n);
      double rho10 = hrhon[3]->GetBinContent(n);
      double rho11e = hrhon[0]->GetBinError(n);
      double rho00e = hrhon[1]->GetBinError(n);
      double rho1Me = hrhon[2]->GetBinError(n);
      double rho10e = hrhon[3]->GetBinError(n);
      double trace = 2 * rho11 + rho00;
      hrhon[0]->SetBinContent(n, rho11 /trace);
      hrhon[1]->SetBinContent(n, rho00 /trace);
      hrhon[2]->SetBinContent(n, rho1M /trace);
      hrhon[3]->SetBinContent(n, rho10 /trace);
      hrhon[0]->SetBinError(n, rho11e /trace);
      hrhon[1]->SetBinError(n, rho00e /trace);
      hrhon[2]->SetBinError(n, rho1Me /trace);
      hrhon[3]->SetBinError(n, rho10e /trace);
   }

   for (int n=0; n < 4; ++n) {
      hrhos[n]->SetDirectory(0);
      hrhon[n]->SetDirectory(0);
   }
   TFile fout("scans.root", "update");
   for (int n=0; n < 4; ++n) {
      hrhos[n]->Write();
      hrhon[n]->Write();
   }
}

void jackknife(mlfitter &fit)
{
   // Use a jack-knife technique to compute the statistical error
   // on the SDME fit results.

   int parts = 10;
   double fpar[4][10] = {0};
   double epar[4][10] = {0};

   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1");

   fit.ClearRfunction();
   TH2D *hR = (TH2D*)gROOT->FindObject("hRfunc0");
   double ethresh = 0.068;
   TString title;
   title.Form("R function on phi vs cos(theta), cut value %f", ethresh);
   hR->SetTitle(title);
   mask2D(hR, select2D(ethresh));
   hR->Draw("colz");
   c1->Update();
   fit.LoadAccept();
   fit.LoadNormInt();
   fit.LoadRealData();
   xyplots(fit);
   double sums[4][4] = {0};
   for (int part=1; part <= parts; ++part) {
      fit.PartitionEvents(-part, parts);
      if (fit.Fit() == 0) {
         for (int n=0; n < 4; ++n) {
            int npar;
            double parval;
            double parerr;
            fit.GetMinuit()->GetParameter(n, parval, parerr);
            fpar[n][part-1] = parval;
            epar[n][part-1] = parerr;
            sums[n][0] += 1;
            sums[n][1] += parval;
            sums[n][2] += parval * parval;
            sums[n][3] += parerr;
         }
      }
   }
   std::cout << setprecision(9);
   std::cout << "Results from jack knife analysis:" << std::endl;
   std::cout << "  rho11   " << "  rho00   " 
             << "  rho1M   " << "  rho10   " << std::endl;
   for (int part=1; part <= parts; ++part) {
      for (int n=0; n < 4; ++n) {
         std::cout << " " << epar[n][part-1];
      }
      std::cout << std::endl;
   }
   std::cout << "mean errors:";
   for (int n=0; n < 4; ++n) {
      sums[n][0] += 1e-99;
      double emean = sums[n][3] / sums[n][0];
      std::cout << " " << emean;
   }
   std::cout << std::endl;
   for (int part=1; part <= parts; ++part) {
      for (int n=0; n < 4; ++n) {
         std::cout << " " << fpar[n][part-1];
      }
      std::cout << std::endl;
   }
   std::cout << "jackknife errors:";
   for (int n=0; n < 4; ++n) {
      double sigma = sqrt(sums[n][2] - sums[n][1]*sums[n][1]/sums[n][0]);
      std::cout << " " << sigma * parts / (parts - 1.);
   }
   std::cout << std::endl;
}

void scan3Rcut(mlfitter &fit,
               int nbins=15, 
               double cutcost0=-1.0,
               double cutcost1=0.8)
{
   // Do three separate scans over the selection thresholds for
   // the R function, and plot the SDME values from the fit as a
   // function of the threshold for 3 separate bins in |t|.

   fit.SetLimits_tabs(0.8, 100);
   scanRcut(fit, nbins, 0.155, 0.305, cutcost0, cutcost1);
   fit.SetLimits_tabs(0.3, 0.8);
   scanRcut(fit, nbins, 0.0655, 0.0805, cutcost0, cutcost1);
   fit.SetLimits_tabs(-100, 0.3);
   scanRcut(fit, nbins, 0.0555, 0.0705, cutcost0, cutcost1);
}

void scan1Rcut(mlfitter &fit,
               int nbins=15, 
               double cutlim0=0.0455,
               double cutlim1=0.0605,
               double cutcost0=-1.0,
               double cutcost1=0.8)
{
   // Do a scan over the selection threshold for the R function,
   // and plot the SDME values from the fit as a function of the
   // threshold for 3 separate bins in |t|. The desired kinematic
   // cuts for the scan must be already set in the mlfitter fit,
   // apart from the limits on |t|.

   TH1D *hrhos[4][3];
   TString tsuffix[3] = {"a", "b", "c"};
   TString ttitle[3] = {"|t|<0.3", "0.3<|t|<0.8", "|t|>0.8"};
   TString pname[4] = {"rho11", "rho00", "rho1M", "rho10"};
   TString ptitle[4] = {"\\rho 11", "\\rho 00", "\\rho 1M", "\\rho 10"};
   for (int tbin=0; tbin < 3; ++tbin) {
      for (int par=0; par < 4; ++par) {
         TString hname;
         TString htitle;
         hname.Form("h%s%s", pname[par].Data(), tsuffix[tbin].Data());
         htitle.Form("%s vs hdcut selection threshold, %s",
                     pname[par].Data(), ttitle[tbin].Data());
         hrhos[par][tbin] = new TH1D(hname, htitle, nbins, cutlim0, cutlim1);
         hrhos[par][tbin]->GetYaxis()->SetTitle(ptitle[par]);
         hrhos[par][tbin]->GetXaxis()->SetTitle("ethresh");
         hrhos[par][tbin]->SetStats(0);
      }
   }

   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1");

   TH2D *hRn[nbins];
   for (int bin=0; bin < nbins; ++bin) {
      fit.SetLimits_tabs(-100., 100.);
      fit.LoadAccept();
      fit.ClearRfunction();
      TString name("hRfunc0");
      TH2D *hR = (TH2D*)gROOT->FindObject(name);
      double ethresh = hrhos[0][0]->GetXaxis()->GetBinCenter(bin+1);
      TString title;
      title.Form("R function on phi vs cos(theta), cut value %f", ethresh);
      hR->SetTitle(title);
      if (fabs(ethresh) > 1e-6) {
         xyplots(fit);
         mask2D(hR, select2D(ethresh));
         fit.ApplyCosthetaCut(cutcost0, cutcost1);
      }
      hR->Draw("colz");
      c1->Update();
      name.Form("hRfunc0_%d", bin+1);
      hRn[bin] = (TH2D*)hR->Clone(name);
      for (int tbin=0; tbin < 3; ++tbin) {
         if (tbin == 0)
            fit.SetLimits_tabs(-100, 0.3);
         else if (tbin == 1)
            fit.SetLimits_tabs(0.3, 0.8);
         else
            fit.SetLimits_tabs(0.8, 100);
         fit.LoadAccept();
         fit.LoadNormInt();
         fit.LoadRealData();
         if (fit.Fit() == 0) {
            for (int par=0; par < 4; ++par) {
               int npar;
               double parval;
               double parerr;
               fit.GetMinuit()->GetParameter(par, parval, parerr);
               hrhos[par][tbin]->SetBinContent(bin+1, parval);
               hrhos[par][tbin]->SetBinError(bin+1, parerr);
            }
         }
      }
   }

   TH1D *hrhon[4][3];
   for (int tbin=0; tbin < 3; ++tbin) {
      for (int par=0; par < 4; ++par) {
         TString hname(hrhos[par][tbin]->GetName());
         hrhon[par][tbin] = copy1D(hrhos[par][tbin], hname + "n");
      }
      for (int n=1; n <= nbins; ++n) {
         double rho11 = hrhon[0][tbin]->GetBinContent(n);
         double rho00 = hrhon[1][tbin]->GetBinContent(n);
         double rho1M = hrhon[2][tbin]->GetBinContent(n);
         double rho10 = hrhon[3][tbin]->GetBinContent(n);
         double rho11e = hrhon[0][tbin]->GetBinError(n);
         double rho00e = hrhon[1][tbin]->GetBinError(n);
         double rho1Me = hrhon[2][tbin]->GetBinError(n);
         double rho10e = hrhon[3][tbin]->GetBinError(n);
         double trace = 2 * rho11 + rho00;
         hrhon[0][tbin]->SetBinContent(n, rho11 / trace);
         hrhon[1][tbin]->SetBinContent(n, rho00 / trace);
         hrhon[2][tbin]->SetBinContent(n, rho1M / trace);
         hrhon[3][tbin]->SetBinContent(n, rho10 / trace);
         hrhon[0][tbin]->SetBinError(n, rho11e / trace);
         hrhon[1][tbin]->SetBinError(n, rho00e / trace);
         hrhon[2][tbin]->SetBinError(n, rho1Me / trace);
         hrhon[3][tbin]->SetBinError(n, rho10e / trace);
      }
      for (int par=0; par < 4; ++par) {
         hrhos[par][tbin]->SetDirectory(0);
         hrhon[par][tbin]->SetDirectory(0);
      }
   }

   TFile fout("scans.root", "update");
   for (int tbin=0; tbin < 3; ++tbin) {
      for (int par=0; par < 4; ++par) {
         hrhos[par][tbin]->Write();
         hrhon[par][tbin]->Write();
      }
   }
   for (int bin=0; bin < nbins; ++bin) {
      hRn[bin]->Write();
   }
}

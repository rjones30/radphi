//////////////////////////////////////////////////////////
//
// xyplots.C - Utility functions for SDME analysis of the
//             decays of vector mesons with Radphi.
//
// author: richard.t.jones at uconn.edu
// version: first written november 14, 2017
//
//////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TCanvas.h>

#include "xyplots.h"

TH2D *triclip_h2 = 0;
int triclip_nclick = 0;
double triclip_x[3];
double triclip_y[3];
double triclip_zero;

void xyplots(mlfitter &fit)
{
   // Generate phiH vs costhetaH 2D histograms of the selected
   // real data (hdata) and fit function (hfit), decorating
   // them with axis labels and readying them for display,
   // together with their 1D projections.

   TH2D *hdata = fit.GenerateDataPlot("hdata");
   hdata->GetXaxis()->SetTitle("\\cos\\theta_H");
   hdata->GetYaxis()->SetTitle("\\phi_H");
   TH1D *hdatax = hdata->ProjectionX("hdatax");
   TH1D *hdatay = hdata->ProjectionY("hdatay");
   hdata->SetStats(0);
   hdatax->SetStats(0);
   hdatay->SetStats(0);
   TH2D *hfit = fit.GenerateFitModel("hfit");
   hfit->GetXaxis()->SetTitle("\\cos\\theta_H");
   hfit->GetYaxis()->SetTitle("\\phi_H");
   hfit->SetLineColor(kRed);
   TH1D *hfitx = hfit->ProjectionX("hfitx");
   TH1D *hfity = hfit->ProjectionY("hfity");
   hfit->SetStats(0);
   hfitx->SetStats(0);
   hfity->SetStats(0);
}

void truncate2D(TH2D *h2, double hmax)
{
   // Find all bins in histogram h2 with contents greater than hmax
   // and replace the contents with hmax. The input histogram is
   // modified in-place.

   int nx = h2->GetNbinsX();
   int ny = h2->GetNbinsY();
   for (int i=1; i <= nx; ++i)
      for (int j=1; j <= ny; ++j)
         if (h2->GetBinContent(i,j) > hmax)
            h2->SetBinContent(i,j,hmax);
}

void zero2D(TH2D *h2, double hmin)
{
   // Find all bins in histogram h2 with contents less than hmin
   // and replace the contents with zero. The input histogram is
   // modified in-place.

   int nx = h2->GetNbinsX();
   int ny = h2->GetNbinsY();
   for (int i=1; i <= nx; ++i)
      for (int j=1; j <= ny; ++j)
         if (h2->GetBinContent(i,j) < hmin)
            h2->SetBinContent(i,j,0);
}

void mask2D(TH2D *h2, const TH2D *hmask)
{
   // Set all of the bins in histograph h2 to zero that have zero in the
   // corresponding bin in hmask. The two histograms must have the same
   // number of bins in x and y.

   int nx = h2->GetNbinsX();
   int ny = h2->GetNbinsY();
   if (hmask->GetNbinsX() != nx || hmask->GetNbinsY() != ny) {
      std::cerr << "Error in mask: two histograms provided as arguments"
                << " do not have the same binning, cannot continue."
                << std::endl;
      return;
   }

   for (int i=1; i <= nx; ++i)
      for (int j=1; j <= ny; ++j)
         if (hmask->GetBinContent(i,j) == 0)
            h2->SetBinContent(i,j,0);
}

void holefill2D(TH2D *h2)
{
   // Find all zero-value pixels in histograph h2 that are surrounded
   // on all 4 sides with non-zero values, and overwrite them with the
   // mean value of their 4 neighbors.

   int nx = h2->GetNbinsX();
   int ny = h2->GetNbinsY();

   for (int i=1; i <= nx; ++i) {
      for (int j=1; j <= ny; ++j) {
         if (h2->GetBinContent(i,j) == 0) {
            double vno = h2->GetBinContent(i,j+1);
            double vso = h2->GetBinContent(i,j-1);
            double vea = h2->GetBinContent(i+1,j);
            double vwe = h2->GetBinContent(i-1,j);
            if (vno > 0 && vso > 0 && vea > 0 && vwe > 0)
               h2->SetBinContent(i,j,(vno+vso+vea+vwe)/4);
         }
      }
   }
}

TH2D *Wplot(mlfitter &fit)
{
   // Create a new histogram called wplot and fill it with the
   // function W(costhetahel,phihel) parameterized from the
   // latest fit, sampled at the centers of each bin. If wplot
   // already exists, the former histogram is overwritten.
 
   TH2D *wplot = (TH2D*)gROOT->FindObject("wplot");
   if (wplot == 0) {
      wplot = (TH2D*)gROOT->FindObject("haccept0")->Clone("wplot");
      wplot->SetTitle("W function");
      wplot->GetXaxis()->SetTitle("\\cos\\theta_H");
      wplot->GetYaxis()->SetTitle("\\phi_H");
      wplot->SetStats(0);
   }
   int nx = wplot->GetNbinsX();
   int ny = wplot->GetNbinsY();
   for (int i=1; i <= nx; ++i) {
      double cost = wplot->GetXaxis()->GetBinCenter(i);
      for (int j=1; j <= ny; ++j) {
         double phi = wplot->GetYaxis()->GetBinCenter(j);
         wplot->SetBinContent(i,j,fit.GetW(cost,phi));
      }
   }
   return wplot;
}

void etch2D(TH2D *h2, double ethresh)
{
   // It is assumed that input histogram h2 contains an inner simply connected
   // region with non-zero contents, which tapers off to zero around the edges.
   // Starting from the outside, this algorithm searches for the edges of this
   // region, "etching" the pixels outside this region to zero and stopping at
   // the edge. The edge is characterized by threshold value ethresh. The
   // histogram is modified in-place.

   int netched=1;
   while (netched > 0) {
      netched = 0;
      for (int j=1; j <= 51; ++j) {
         for (int i=1; i <= 51; ++i) {
            double e = h2->GetBinContent(i,j);
            if (e > 0 && e < ethresh) {
               if ((h2->GetBinContent(i,j+1) < ethresh * 1.2 &&
                    h2->GetBinContent(i,j-1) == 0) ||
                   (h2->GetBinContent(i,j-1) < ethresh * 1.2 &&
                    h2->GetBinContent(i,j+1) == 0))
               {
                  h2->SetBinContent(i,j,0);
                  ++netched;
               }
               else if (h2->GetBinContent(i,j+1) < ethresh &&
                        h2->GetBinContent(i,j-1) < ethresh &&
                        h2->GetBinContent(i-1,j-1) == 0 &&
                        h2->GetBinContent(i-1,j-1) == 0)
               {
                  h2->SetBinContent(i,j,0);
                  ++netched;
               }   
               else
                  break;
            }
            else if (e > 0)
               break;
         }
         for (int i=50; i >= 1; --i) {
            double e = h2->GetBinContent(i,j);
            if (e > 0 && e < ethresh) {
               if ((h2->GetBinContent(i,j+1) < ethresh * 1.2 &&
                    h2->GetBinContent(i,j-1) == 0) ||
                   (h2->GetBinContent(i,j-1) < ethresh * 1.2 &&
                    h2->GetBinContent(i,j+1) == 0))
               {
                  h2->SetBinContent(i,j,0);
                  ++netched;
               }
               else if (h2->GetBinContent(i,j+1) < ethresh &&
                        h2->GetBinContent(i,j-1) < ethresh &&
                        h2->GetBinContent(i+1,j-1) == 0 &&
                        h2->GetBinContent(i+1,j-1) == 0)
               {
                  h2->SetBinContent(i,j,0);
                  ++netched;
               }
               else
                  break;
            }
            else if (e > 0)
               break;
         }
      }
   }
}

TH2D *select2D(double ethresh, double athresh)
{
   // Apply an acceptance correction to the real data decay angles
   // distribution hdata (see xyplots above) and then etch the
   // boundaries of this plot up to threshold value ethresh. If
   // no ethresh is given, display an interactive loop to allow
   // the user to chose the best threshold.

   TH2D *haccept0 = (TH2D*)gROOT->FindObject("haccept0");
   TH2D *hdata = (TH2D*)gROOT->FindObject("hdata");
   if (haccept0 == 0 || hdata == 0) {
      std::cerr << "Error in select2D - must call xyplots"
                << " before invoking select2D."
                << std::endl;
      return 0;
   }
   TH2D *hacor = copy2D(haccept0, "hacor");
   zero2D(hacor, athresh * haccept0->GetMaximum());
   TH2D *hdcor = copy2D(hdata, "hdcor");
   zero2D(hdcor, 1e-9);
   hdcor->Divide(hacor);
   TH2D *hdsel = copy2D(hdata, "hdsel");

   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1", "c1");
   if (ethresh == 0) {
      TString ans("no");
      while (true) {
         etch2D(hdcor, ethresh);
         mask2D(hdsel, hdcor);
         hdcor->Draw("colz");
         c1->Update();
         std::cout << "ethresh=" << ethresh 
                << ", event fraction: "
                << hdsel->Integral() / hdata->Integral()
                << ", ethresh, +/-delta, or q to quit: ";
         std::cin >> ans;
         if (ans(0) == '+' || ans(0) == '-')
            ethresh += ans.Atof();
         else if (ans.IsFloat())
            ethresh = ans.Atof();
         else if (ans == "q")
            break;
         hacor = copy2D(haccept0, "hacor");
         zero2D(hacor, athresh * haccept0->GetMaximum());
         hdcor = copy2D(hdata, "hdcor");
         zero2D(hdcor, 1e-9);
         hdcor->Divide(hacor);
         hdsel = copy2D(hdata, "hdsel");
      }
   }
   else {
      etch2D(hdcor, ethresh);
      mask2D(hdsel, hdcor);
      hdcor->Draw("colz");
   }
   return hdcor;
}

TH2D *copy2D(const TH2D *h2, const TString &hname, const TString *title)
{
   // Make a copy of 2D histogram *h2, name the new copy hname,
   // and delete any previous objects with that name.

   TObject *hold;
   while ((hold = gROOT->FindObject(hname)))
      hold->Delete();
   TH2D *hnew  = (TH2D*)h2->Clone(hname);
   if (title)
      hnew->SetTitle(*title);
   return hnew;
}

TH1D *copy1D(const TH1D *h1, const TString &hname, const TString *title)
{
   // Make a copy of 1D histogram *h2, name the new copy hname,
   // and delete any previous objects with that name.

   TObject *hold;
   while ((hold = gROOT->FindObject(hname)))
      hold->Delete();
   TH1D *hnew  = (TH1D*)h1->Clone(hname);
   if (title)
      hnew->SetTitle(*title);
   return hnew;
}

TH2D *triclip2D(TH2D *h2, double zero)
{
   // Display 2D histogram h2 and wait for the user to left-click
   // three points inside the displayed box. These three points are
   // the three corners of a triangle whose contents are replaced
   // with the value zero, after which the histogram display is
   // updated. This continues until another plot is displayed.
 
   h2->Draw("colz");
   triclip_h2 = h2;
   triclip_nclick = 0;
   triclip_zero = zero;
   gPad->AddExec("triclip_callback","triclip_callback()");
   std::cout << "To complete this action, click three corners"
             << " of a triangle on the plot." 
             << std::endl;
   return h2;
}

void triclip_callback()
{
   // Helper action for the triclip2D function.

   int event = gPad->GetEvent();
   if (event == 1) { // left mouse button click
      int px = gPad->GetEventX();
      int py = gPad->GetEventY();
      double xd = gPad->AbsPixeltoX(px);
      double yd = gPad->AbsPixeltoY(py);
      triclip_x[triclip_nclick] = gPad->PadtoX(xd);
      triclip_y[triclip_nclick] = gPad->PadtoY(yd);
      ++triclip_nclick;
   }
   if (triclip_nclick == 3) {
      TH2D *h2 = triclip_h2;
      double vx[3];
      double vy[3];
      for (int n=0; n < 3; ++n) {
         int ix = triclip_h2->GetXaxis()->FindBin(triclip_x[n]);
         int iy = triclip_h2->GetYaxis()->FindBin(triclip_y[n]);
         vx[n] = triclip_h2->GetXaxis()->GetBinCenter(ix);
         vy[n] = triclip_h2->GetYaxis()->GetBinCenter(iy);
      }
      double phi[2][2];
      double deltaphi[2];
      double phi0, phi1, dphi;
      phi0 = atan2(vy[1] - vy[0], vx[1] - vx[0]);
      phi1 = atan2(vy[2] - vy[0], vx[2] - vx[0]);
      dphi = phi1 - phi0;
      if (dphi >= 0 && dphi < M_PI) {
         phi[0][0] = phi0;
         phi[0][1] = phi1;
         deltaphi[0] = dphi;
      }
      else if (dphi >= M_PI) {
         phi[0][0] = phi1;
         phi[0][1] = phi0;
         deltaphi[0] = 2*M_PI - dphi;
      }
      else if (dphi < 0 && dphi > -M_PI) {
         phi[0][0] = phi1;
         phi[0][1] = phi0;
         deltaphi[0] = -dphi;
      }
      else {
         phi[0][0] = phi0;
         phi[0][1] = phi1;
         deltaphi[0] = 2*M_PI + dphi;
      }
      phi0 = atan2(vy[2] - vy[1], vx[2] - vx[1]);
      phi1 = atan2(vy[0] - vy[1], vx[0] - vx[1]);
      dphi = phi1 - phi0;
      if (dphi >= 0 && dphi < M_PI) {
         phi[1][0] = phi0;
         phi[1][1] = phi1;
         deltaphi[1] = dphi;
      }
      else if (dphi >= M_PI) {
         phi[1][0] = phi1;
         phi[1][1] = phi0;
         deltaphi[1] = 2*M_PI - dphi;
      }
      else if (dphi < 0 && dphi > -M_PI) {
         phi[1][0] = phi1;
         phi[1][1] = phi0;
         deltaphi[1] = -dphi;
      }
      else {
         phi[1][0] = phi0;
         phi[1][1] = phi1;
         deltaphi[1] = 2*M_PI + dphi;
      }
      for (int ix=1; ix <= h2->GetNbinsX(); ++ix) {
         for (int iy=1; iy <= h2->GetNbinsY(); ++iy) {
            double x = h2->GetXaxis()->GetBinCenter(ix);
            double y = h2->GetYaxis()->GetBinCenter(iy);
            if ( (x == vx[0] && y == vy[0]) ||
                 (x == vx[1] && y == vy[1]) ||
                 (x == vx[2] && y == vy[2]) )
            {
               h2->SetBinContent(ix, iy, triclip_zero);
               continue;
            }
            phi0 = atan2(y - vy[0], x - vx[0]);
            phi1 = atan2(y - vy[1], x - vx[1]);
            if (((phi0 >= phi[0][0] && phi0 <= phi[0][0] + deltaphi[0]) ||
                 (phi0 <= phi[0][1] && phi0 >= phi[0][1] - deltaphi[0])) && 
                ((phi1 >= phi[1][0] && phi1 <= phi[1][0] + deltaphi[1]) ||
                 (phi1 <= phi[1][1] && phi1 >= phi[1][1] - deltaphi[1])))
            {
               h2->SetBinContent(ix, iy, triclip_zero);
            }
         }
      }
      triclip_nclick = 0;
      gPad->DeleteExec("triclip_callback");
      TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
      h2->Draw("colz");
      c1->Update();
   }
}

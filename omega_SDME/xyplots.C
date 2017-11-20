#include <TROOT.h>
#include <TH2D.h>
#include <TCanvas.h>

#include "mlfitter.h"

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
               if ((h2->GetBinContent(i,j+1) < ethresh &&
                    h2->GetBinContent(i,j-1) == 0) ||
                   (h2->GetBinContent(i,j-1) < ethresh &&
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
               if ((h2->GetBinContent(i,j+1) < ethresh &&
                    h2->GetBinContent(i,j-1) == 0) ||
                   (h2->GetBinContent(i,j-1) < ethresh &&
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

TH2D *select2D(double ethresh=0)
{
   // Apply an acceptance correction to the real data decay angles
   // distribution hdata (see xyplots above) and then etch the
   // boundaries of this plot up to threshold value ethresh. If
   // no ethresh is given, display an interactive loop to allow
   // the user to chose the best threshold.

   TH2D *haccept0 = (TH2D*)gROOT->FindObject("haccept0");
   TH2D *hdata = (TH2D*)gROOT->FindObject("hdata");
   if (haccept0 == 0 || hdata == 0) {
      std::cerr << "Error in genR2D - must call mlfitter::LoadAccept"
                << " and LoadRealData() before invoking genR2D."
                << std::endl;
      return 0;
   }
   TH2D *hacor = (TH2D*)haccept0->Clone("hacor");
   zero2D(hacor, 0.1 * haccept0->GetMaximum());
   TH2D *hdcor = (TH2D*)hdata->Clone("hdcor");
   zero2D(hdcor, 1e-9);
   hdcor->Divide(hacor);
   TH2D *hdsel = (TH2D*)hdata->Clone("hdsel");

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
         hacor->Delete();
         hacor = (TH2D*)haccept0->Clone("hacor");
         zero2D(hacor, 0.1 * haccept0->GetMaximum());
         hdcor->Delete();
         hdcor = (TH2D*)hdata->Clone("hdcor");
         zero2D(hdcor, 1e-9);
         hdcor->Divide(hacor);
         hdsel->Delete();
         hdsel = (TH2D*)hdata->Clone("hdsel");
      }
   }
   else {
      etch2D(hdcor, ethresh);
      mask2D(hdsel, hdcor);
      hdcor->Draw("colz");
   }
   return hdcor;
}

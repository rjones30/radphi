//////////////////////////////////////////////////////////
//
// xyplots.C - Utility functions for SDME analysis of the
//             decays of vector mesons with Radphi.
//
// author: richard.t.jones at uconn.edu
// version: first written november 14, 2017
//
//////////////////////////////////////////////////////////

#include <TH2D.h>

#include "mlfitter.h"

void xyplots(mlfitter &fit);
void truncate2D(TH2D *h2, double hmax);
void zero2D(TH2D *h2, double hmin);
void mask2D(TH2D *h2, const TH2D *hmask);
void holefill2D(TH2D *h2);
TH2D *Wplot(mlfitter &fit);
void etch2D(TH2D *h2, double ethresh);
TH2D *select2D(double ethresh=0, double athresh=0.1);
TH2D *copy2D(const TH2D *h2, const TString &hname, const TString *title=0);
TH1D *copy1D(const TH1D *h1, const TString &hname, const TString *title=0);
TH2D *triclip2D(TH2D *h2, double zero=0);

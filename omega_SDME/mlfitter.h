//////////////////////////////////////////////////////////
//
// mlfitter - A helper class for doing maximum likelihood
//            unbinned fit of experimental data describing
//            2-body decays from a vector to pi0, gamma.
//
// author: richard.t.jones at uconn.edu
// version: first written november 14, 2017
//
//////////////////////////////////////////////////////////

#ifndef mlfitter_h
#define mlfitter_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TMinuit.h>
#include <TRandom3.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>

#include <vector>


class mlfitter {
 public :
   mlfitter();
   ~mlfitter();

   void LoadAccept();
   void LoadNormInt();
   void LoadRealData(int nevents=0, int nstart=0);

   TMinuit *GetMinuit();
   int Do(const char *cmd);
   int Do(const char *cmd, double arg1);
   int Do(const char *cmd, double arg1, double arg2);
   int Do(const char *cmd, double arg1, double arg2, double arg3);
   int Do(const char *cmd, double arg1, double arg2, double arg3, double arg4);
   int Fit(int redundancy=2, int maxtries=9);

   double GetW(double costheta, double phi) const;
   double W11(double costheta, double phi) const;
   double W00(double costheta, double phi) const;
   double W1m(double costheta, double phi) const;
   double W10(double costheta, double phi) const;

   TChain *ChainMCsim(const TString templ);
   TChain *ChainMCgen(const TString templ);
   TChain *ChainRealData(const TString templ);

   double GetFCN(double rho11, double rho00, double rh1m, double rho10);

   double GetLastFCN() const;
   double GetLastParam(unsigned int i) const;
   const double *GetLastParam() const;
   double GetLastGradient(unsigned int i) const;
   const double *GetLastGradient() const;
   double GetLastHesse(unsigned int i, unsigned int j) const;
   const double *GetLastHesse() const;

   void SetLimits_Ebeam(double Ebeamlim0, double Ebeamlim1);
   void SetLimits_Emiss(double Emisslim0, double Emisslim1);
   void SetLimits_tabs(double tlim0, double tlim1);
   void SetLimits_mass(double masslim0, double masslim1);
   void SetLimits_dphi(double dphilim0, double dphilim1);
   
   double GetR(double costheta, double phi) const;
   double GetA(double costheta, double phi) const;

   TH2D *GenerateDataPlot(const TString &hname);
   TH2D *GenerateFitModel(const TString &hname);

   void RandomizeParams();
   void LoadRfunction();
   void ClearRfunction();
   void ApplyCosthetaCut(double costlim0, double costlim1);
   void ApplyAcceptanceThreshold(double acclim0);

   void SaveState(const char *filename) const;
   void RestoreState(const char *filename);
   void Print() const;

   struct event_t {
      float W11;
      float W00;
      float W1m;
      float W10;
      float wgt;
      float R;
   };

   double fRule1 = 0;
   double fRule2 = 0;
   double fRule3 = 0;
   double fRule4 = 0;
   double fWlowcut = 1e-20;

 private:
   mlfitter(const mlfitter &src) {}

 protected:
   TMinuit *fMinuit = 0;
   TRandom3 *fRandom = 0;

   static void FCN(int &npar, double *gin, double &f, double *par, int iflag);
   // For multi-threaded operation this needs to be thread_local, but
   // it causes a crash in Root 6.08 so until this is fixed, comment it.
   static /*thread_local*/ mlfitter *fFitter;

   // Fit parameters are ordered: rho11, rho00, rho1m, rho10
   double fLastFCN = 0;
   double fLastParam[4] = {0};
   double fLastGradient[4] = {0};
   double fLastHesse[4][4] = {{0}};
   double fEbeamlim[2];
   double fEmisslim[2];
   double ftabslim[2];
   double fmasslim[2];
   double fdphilim[2];

   std::vector<struct event_t> fSample;
   double fNormInt11 = 0;
   double fNormInt00 = 0;
   double fNormInt1m = 0;
   double fNormInt10 = 0;
   double fNgenerated = 0;
   TH2D *fHaccept = 0;
   TH2D *fHavewgt = 0;
   TH2D *fHavewgt2 = 0;
   TH2D *fRfunction = 0;

   TTree *fChain = 0;
   TTreeReader fReader;
   TTreeReaderValue<Double_t> runevent = {fReader, "runevent"};
   TTreeReaderValue<Double_t> Ebeam = {fReader, "Ebeam"};
   TTreeReaderValue<Double_t> Emiss = {fReader, "Emiss"};
   TTreeReaderValue<Double_t> m3g = {fReader, "m3g"};
   TTreeReaderValue<Double_t> dphi = {fReader, "dphi"};
   TTreeReaderValue<Double_t> tabs = {fReader, "tabs"};
   TTreeReaderValue<Double_t> costhetahel = {fReader, "costhetahel"};
   TTreeReaderValue<Double_t> phihel = {fReader, "phihel"};
   TTreeReaderValue<Double_t> weight = {fReader, "weight"};
};

inline double mlfitter::GetLastFCN() const {
   return fLastFCN;
}

inline double mlfitter::GetLastParam(unsigned int i) const {
   if (i < 4)
      return fLastParam[i];
   else
      return 0;
}

inline double mlfitter::GetLastGradient(unsigned int i) const {
   if (i < 4)
      return fLastGradient[i];
   else
      return 0;
}

inline const double *mlfitter::GetLastParam() const {
   return fLastParam;
}

inline const double *mlfitter::GetLastGradient() const {
   return fLastGradient;
}

inline double mlfitter::GetLastHesse(unsigned int i, unsigned int j) const {
   if (i < 4 && j < 4)
      return fLastHesse[i][j];
   else
      return 0;
}

inline const double *mlfitter::GetLastHesse() const {
   return fLastHesse[0];
}

inline void mlfitter::SetLimits_tabs(double tlim0, double tlim1) {
   ftabslim[0] = tlim0;
   ftabslim[1] = tlim1;
}

inline void mlfitter::SetLimits_Ebeam(double Ebeamlim0, double Ebeamlim1) {
   fEbeamlim[0] = Ebeamlim0;
   fEbeamlim[1] = Ebeamlim1;
}

inline void mlfitter::SetLimits_Emiss(double Emisslim0, double Emisslim1) {
   fEmisslim[0] = Emisslim0;
   fEmisslim[1] = Emisslim1;
}

inline void mlfitter::SetLimits_mass(double masslim0, double masslim1) {
   fmasslim[0] = masslim0;
   fmasslim[1] = masslim1;
}

inline void mlfitter::SetLimits_dphi(double dphilim0, double dphilim1) {
   fdphilim[0] = dphilim0;
   fdphilim[1] = dphilim1;
}

inline double mlfitter::GetR(double costheta, double phi) const {
   if (fRfunction)
      return fRfunction->Interpolate(costheta, phi);
   else
      return 0;
}

inline double mlfitter::GetA(double costheta, double phi) const {
   if (fHaccept)
      return fHaccept->Interpolate(costheta, phi);
   else
      return 0;
}

inline double mlfitter::GetW(double costheta, double phi) const {
   return fLastParam[0] * W11(costheta, phi) +
          fLastParam[1] * W00(costheta, phi) +
          fLastParam[2] * W1m(costheta, phi) +
          fLastParam[3] * W10(costheta, phi);
}

inline int mlfitter::Do(const char *cmd) {
   if (fMinuit == 0)
      fMinuit = GetMinuit();
   fFitter = this;
   int icondn;
   fMinuit->mncomd(cmd, icondn);
   return icondn;
}

inline int mlfitter::Do(const char *cmd, double arg1) {
   if (fMinuit == 0)
      fMinuit = GetMinuit();
   fFitter = this;
   int icondn;
   double arglist[] = {arg1};
   fMinuit->mnexcm(cmd, arglist, 1, icondn);
   return icondn;
}

inline int mlfitter::Do(const char *cmd, double arg1, double arg2) {
   if (fMinuit == 0)
      fMinuit = GetMinuit();
   fFitter = this;
   int icondn;
   double arglist[] = {arg1, arg2};
   fMinuit->mnexcm(cmd, arglist, 2, icondn);
   return icondn;
}

inline int mlfitter::Do(const char *cmd, double arg1, double arg2, double arg3) {
   if (fMinuit == 0)
      fMinuit = GetMinuit();
   fFitter = this;
   int icondn;
   double arglist[] = {arg1, arg2, arg3};
   fMinuit->mnexcm(cmd, arglist, 3, icondn);
   return icondn;
}

inline int mlfitter::Do(const char *cmd, double arg1, double arg2, double arg3, double arg4) {
   if (fMinuit == 0)
      fMinuit = GetMinuit();
   fFitter = this;
   int icondn;
   double arglist[] = {arg1, arg2, arg3, arg4};
   fMinuit->mnexcm(cmd, arglist, 4, icondn);
   return icondn;
}

#endif

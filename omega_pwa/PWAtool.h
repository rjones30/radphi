//
// class PWAtool - provides the quantity to be minimized in the
//                 extended maximum likelihood method for the
//                 partial-wave analysis of omega decay angular
//                 distributions from the Radphi experiment
//
// authors: richard.t.jones at uconn.edu, fridah.mokaya at uconn.edu
// version: april 5, 2016
//
// For implementation details, see the source file PWAtool.cc
//

#include <vector>
#include <TChain.h>
#include <TMinuit.h>
#include <TMatrixD.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TRandom3.h>

#ifndef PWATOOL_H
#define PWATOOL_H 1

class PWAtool : public TObject
{
 public:
   PWAtool();
   ~PWAtool();

   struct angles_t {
      double costheta; // range [-1, 1]
      double phi;      // range [-pi, pi]
   };
   typedef struct angles_t angles_T;
   struct event_t {
      double costheta; // range [-1, 1]
      double phi;      // range [-pi, pi]
      double weight;
   };
   typedef struct event_t event_T;
   struct treerow_t {
      double tabs;
      double Ebeam;
      double omegaphi;
      double helcosthe;
      double helthe;
      double helphi;
      double gjcosthe;
      double gjthe;
      double gjphi;
      double wgt;
      double wgt_mass;
      double wgt_delphi;
      double Ewgt;
      double Ewgt_mass;
      double Ewgt_delphi;
   };
   typedef struct treerow_t treerow_T;
   struct mctreerow_t {
      double tabs;
      double Ebeam;
      double omegaphi;
      double helcosthe;
      double helthe;
      double helphi;
      double gjcosthe;
      double gjthe;
      double gjphi;
      double wgt;
      double wgt_mass;
      double wgt_delphi;
      double Ewgt;
      double Ewgt_mass;
      double Ewgt_delphi;
   };
   typedef struct mctreerow_t mctreerow_T;

   double W(const angles_T &Omega);

   double get_sfact() const;
   double get_Re_rho(int m1, int m2) const;
   double get_Im_rho(int m1, int m2) const;
   void set_sfact(double s=0);
   void set_rho(double sdm11, double sdm00, double sdmM1, 
                double Re_sdm10, double Im_sdm10=0);
   void set_parameters13(double param[4]);
   void set_parameters04(double param[4]);

   void print_rho() const;
   void print_last_param();
   const TMatrixD &get_PartialD() { return PartialD; }

   double get_phi_cut0() const { return phi_cut[0]; }
   double get_phi_cut1() const { return phi_cut[1]; }
   double get_costheta_cut0() const { return costheta_cut[0]; }
   double get_costheta_cut1() const { return costheta_cut[1]; }
   void set_phi_cuts(double cut0=-10, double cut1=+10) {
      phi_cut[0] = cut0;
      phi_cut[1] = cut1;
   }
   void set_costheta_cuts(double cut0=-10, double cut1=+10) {
      costheta_cut[0] = cut0;
      costheta_cut[1] = cut1;
   }
   unsigned int get_slice() const { return Nslice; }
   unsigned int get_slices() const { return Nslices; }
   void set_slicing(unsigned int slice=0, unsigned int slices=1) {
      Nslice = slice;
      Nslices = slices;
   }
   unsigned int get_stride_step() const { return Nstride_step; }
   unsigned int get_stride_start() const { return Nstride_start; }
   void set_stride(int start=0, int step=1) {
      Nstride_start = start;
      Nstride_step = step;
   }

   double get_real_Nevents(int powgt=0) const;
   double get_mc_Nevents(int powgt=0) const;
   void gen_real_events(int nevents);
   void gen_mc_events(int nevents);
   const event_T &get_real_event(unsigned int i) const;
   const event_T &get_mc_event(unsigned int i) const;

   TH1D *hist_real_costheta(const char *name, double phimin=-M_PI,
                                              double phimax=M_PI);
   TH1D *hist_real_phi(const char *name, double costmax=-1,
                                         double costmin=1);
   TH1D *hist_mc_costheta(const char *name, int weighted=0, double phimin=-M_PI,
                                                            double phimax=M_PI);
   TH1D *hist_mc_phi(const char *name, int weighted=0, double costmin=-1,
                                                       double costmax=1);
   TH2D *hist_real_weight(const char *name, int pow);
   TH2D *hist_mc_weight(const char *name, int pow);

   const TH2D *get_weight_2dhist() const { return weight_2dhist; }
   void set_weight_2dhist(TH2D *h2d) { weight_2dhist = h2d; }
   
   double tf1_W(double *var, double *par);

   // configure class functor to serve as TMinuit fcn //
   double operator()(int &npar, double grad[4], int iflag);
   static void FCN(int &npar,double grad[4],  double &nll,
                   double param[4], int iflag);
   void SetFCN(TMinuit &fitter);

   int load_real_events(TChain *data,
                        double tabsmin=0, double tabsmax=99,
                        double Emin=-1, double Emax=99);
   int load_mc_events(TChain *mc,
                      double tabsmin=0, double tabsmax=99,
                      double Emin=-1, double Emax=99);

 protected:
   double sfact;
   double rho11;
   double rho00;
   double rhoM1;
   double Re_rho10;
   double Im_rho10;
   double rhoMM() const { return rho11; }
   double rho1M() const { return rhoM1; }
   double Re_rho01() const { return +Re_rho10; }
   double Im_rho01() const { return -Im_rho10; }
   double Re_rho0M() const { return -Re_rho10; }
   double Im_rho0M() const { return +Im_rho10; }
   double Re_rhoM0() const { return -Re_rho10; }
   double Im_rhoM0() const { return -Im_rho10; }
   double a[3];
   double b[3];
   double c[3];
   double Wprime[4];
   TMatrixD PartialD;

   double phi_cut[2];
   double costheta_cut[2];
   unsigned int Nslice;
   unsigned int Nslices;
   unsigned int Nstride_step;
   unsigned int Nstride_start;

   TH2D *weight_2dhist;
   std::vector<event_T> real_events;
   std::vector<event_T> mc_events;

   ClassDef(PWAtool, 1);

 private:
   PWAtool(PWAtool& src); // copy constructor not supported
   static PWAtool* owner;
   static TRandom3* randoms;
};

#endif

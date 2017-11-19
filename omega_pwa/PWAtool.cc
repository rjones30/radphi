//
// class PWAtool - provides the quantity to be minimized in the
//                 extended maximum likelihood method for the
//                 partial-wave analysis of omega decay angular
//                 distributions from the Radphi experiment.
//
// authors: richard.t.jones at uconn.edu, fridah.mokaya at uconn.edu
// version: april 5, 2016
//
// implementation notes:
// 1) To use this class, the following steps are required:
//     * construct a new PWAtool instance and call methods
//       load_real_events and load_mc_events to prepare for
//       subsequent calls to the functor.
//     * construct a new TMinuit instance and call its SetFCN
//       method with a pointer to the PWAtool functor created
//       in the previous step.
//     * tell your TMinuit to do the fit, and examine the results.
// 2) The specific reaction that is coded into this class implementation
//    is gamma(unpolarized),p -> omega(pi0,gamma),p where the direction
//    of the pi0 in the omega decay rest frame is given by the pair
//    (costheta,phi). The decay channel-specific code is all local to
//    the W() method, and so it would be straight-forward to extend
//    it to cover other decay channels or to the decays of other 
//    vector mesons.
// 3) For carrying out the full PWA project, it is expected that the
//    user will want to write utility functions to handle the loading
//    of the data from various sources, execution of the fits, plotting
//    of results, and estimation of errors. This is most conveniently
//    done in python -- see example code in PWAomega.py for this.

#define VERBOSE 1

#include "PWAtool.h"
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>

//#define FIT_PARAMETERIZATION_13 1

PWAtool *PWAtool::owner = 0;
TRandom3 *PWAtool::randoms = 0;

PWAtool::PWAtool() 
 : PartialD(4, 4),
   weight_2dhist(0),
   sfact(1)
{
   // takes ownership of FCN

   owner = this;
   if (randoms == 0)
      randoms = new TRandom3();

   // initialize the angle cut windows to include everything
   set_stride();
   set_slicing();
   set_phi_cuts();
   set_costheta_cuts();

   // random initial values for fit parameters
   double pars[] = {0.5788, -0.3415, -0.8993, 1.0082};
   set_parameters04(pars);
}

PWAtool::~PWAtool()
{ 
   // empty destructor

   if (owner == this)
      owner = 0;
}

double PWAtool::W(const angles_T &Omega)
{
   // compute W = sum_{m,m'} [A*_m(Omega) rho(m,m') A_m'(Omega)]
   // where A_m(Omega) is the quantum amplitude for decay from
   // omega state m into the final two-body direction described
   // by angles_t object Omega. Note that the condition
   // trace{rho} = 1 is NOT assumed in this model.
 
   double cost = Omega.costheta;
   double sint = sqrt(1 - cost*cost);
   double phi = Omega.phi;
   Wprime[0] = (1 + cost*cost);
   Wprime[1] = (2 - Wprime[0]);
   Wprime[2] = (2 - Wprime[0]) * cos(2*phi);
   Wprime[3] = 2*sqrt(2.) * sint*cost * cos(phi);
   
   double  w_ang = Wprime[0] * rho11    +
                   Wprime[1] * rho00    + 
                   Wprime[2] * rhoM1    +
                   Wprime[3] * Re_rho10;
   return sfact * w_ang;
}

void PWAtool::set_sfact(double s)

{
   // directly load the sfact constant that normalizes W,
   // or else set sfact = Tr{rho} and renomarlize rho
   // to have unit trace, leaving W invariant.

   if (s > 0) {
      sfact = s;
   }
   else {
      double trace_rho = 2 * rho11 + rho00;
      rho11 /= trace_rho;
      rho00 /= trace_rho;
      rhoM1 /= trace_rho;
      Re_rho10 /= trace_rho;
      sfact *= trace_rho;
   }
} 

void PWAtool::set_rho(double sdm11, double sdm00, double sdmM1, 
                      double Re_sdm10, double Im_sdm10)
{
   // load from arguments directly into the density matrix,
   // no check is made that the values actually satisfy the
   // positivity requirements for W.

   rho11 = sdm11;
   rho00 = sdm00;
   rhoM1 = sdmM1;
   Re_rho10 = Re_sdm10;
   Im_rho10 = Im_sdm10;
}

void PWAtool::set_parameters04(double param[4])
{
   // load from arguments directly into the density matrix,
   // using a parameterization that guarantees the positivity
   // and parity conservation constraints for rho in the decays
   // of a vector meson. The parameterization is as follows.
   //
   //                _3_ 
   //                \     / a_i \*
   //    rho[m,m'] = /__   | c_i |  ( a_i  c_i  b_i )
   //                i=1   \ b_i /
   //
   // where the spinors ( a_i  c_i  b_i ) are a complete set
   // of mutually orthogonal basis vectors whose normalization
   // is not fixed in the parameterization. This construction
   // guarantees the positivity of the W pdf for any set of
   // complex coefficients {a_i, b_i, c_i}. There are 18 free
   // parameters in this set, but this set gets reduced to 9
   // through the following conventions:
   //   1) arbitrary phase on each basis vector (-3)
   //   2) orthogonality conditions (-6)
   // Defining c_i as pure real takes care of (1) and (2) must
   // be enforced by constraining some of the coefficients from
   // the others.
   //
   // Parity conservation further requires that
   //      rho[-m,-m'] = (-1)^(m-m') rho[m,m']
   // which is satisfied by either one of two conditions that
   // must hold for each i:
   //      (a) a_i = -b_i, for any c_i, or
   //      (b) a_i = +b_i, for c_i = 0
   // The parity constraints reduce the number of free parameters
   // from 9 to 5.
   //
   // We select a_1, c_1, a_2, and c_3 as the parameters: a_1 has
   // two degrees of freedom because it is complex, whereas a_2 has
   // a spurious phase and so it is taken to be real. In the special
   // case of an unpolarized beam, it turns out that only the real
   // elements of rho contribute to the differential cross section,
   // so we furthermore restrict a_1 to be real. Thus we end up with
   // a complete parameterization of rho which simultaneously ensures
   // constraints from positivity and parity conservation, plus the
   // unpolarized beam condition, for arbitrary values of the inputs.
   //
   // WARNING: This parameterization has a built-in dualism: two
   // completely different parameter vectors lead to exactly the same
   // rho matrix. One of the two choices typically leads to very poor
   // convergence because the scales of one of the parameters is very
   // far from the others. The equivalence transformation is as follows.
   //    a'_1 = c_1 c_3, c'_1 = -2 a_1 c_3, c'_3 = -1 / (2 c_3)
   // If the scale of the c_3 parameter gets far from the scales of the
   // a_1 and c_1, switching from one choice to the other can improve
   // the convergence rate of the fit.

   a[0] = param[0];
   c[0] = param[1];
   b[0] = -param[0];
   a[1] = param[2];
   c[1] = 0;
   b[1] = param[2];
   a[2] = param[1]*param[3];
   c[2] = -2*param[0]*param[3];
   b[2] = -param[1]*param[3];

   // compute the 4 independent rho matrix elements from the parameters
   rho11          = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
   rho00          = c[0]*c[0]             + c[2]*c[2];
   rhoM1          = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
   Re_rho10       = b[0]*c[0]             + b[2]*c[2];
   Im_rho10       = 0;

   // evaluate the partial derivatives of the rho matrix elements
   // with respect to the parameters, ordered as PartialD[row_i][param_j] 
   //     row_i = 0 (rho11), 1 (rho00), 2 (rhoM1), or 3 (Re_rho10); and
   PartialD[0][0] = 2*a[0];
   PartialD[0][1] = 2*a[2]*param[3];
   PartialD[0][2] = 2*a[1];
   PartialD[0][3] = 2*a[2]*param[1];
   PartialD[1][0] = -4*c[2]*param[3];
   PartialD[1][1] = 2*c[0];
   PartialD[1][2] = 0;
   PartialD[1][3] = -4*c[2]*param[0];
   PartialD[2][0] = -2*a[0];
   PartialD[2][1] = -2*a[2]*param[3];
   PartialD[2][2] = 2*a[1];
   PartialD[2][3] = -2*a[2]*param[1];
   PartialD[3][0] = -c[0] - 2*b[2]*param[3];
   PartialD[3][1] = b[0] - c[2]*param[3];
   PartialD[3][2] = 0;
   PartialD[3][3] = 4*a[2]*param[0];
}

void PWAtool::set_parameters13(double param[4])
{
   // same as set_parameters04, except that the meaning of the param
   // array argument is changed to separate out one scale factor
   // param[0] from the other three which determine the shape of the
   // spin density matrix. The mapping from these param values to
   // the ones needed for the equivalent call to set_parameters04
   // is given in the first few lines below.

   double par04[4];
   par04[0] = param[0];
   par04[1] = param[1] * param[0];
   par04[2] = param[2] * param[0];
   par04[3] = param[3] * param[0];
   set_parameters04(par04);

   // propagate the derivative matrix to the new parameterization
   PartialD[0][0] += PartialD[0][1] * param[1] +
                     PartialD[0][2] * param[2] +
                     PartialD[0][3] * param[3];
   PartialD[1][0] += PartialD[1][1] * param[1] +
                     PartialD[1][2] * param[2] +
                     PartialD[1][3] * param[3];
   PartialD[2][0] += PartialD[2][1] * param[1] +
                     PartialD[2][2] * param[2] +
                     PartialD[2][3] * param[3];
   PartialD[3][0] += PartialD[3][1] * param[1] +
                     PartialD[3][2] * param[2] +
                     PartialD[3][3] * param[3];
   for (int i=0; i < 4; ++i)
      for (int j=1; j < 4; ++j)
         PartialD[i][j] *= param[0];
}

void PWAtool::print_last_param()
{
   printf("%12.7f  %12.7f  %12.7f\n", a[0], a[1], a[2]);
   printf("%12.7f  %12.7f  %12.7f\n", c[0], c[1], c[2]);
   printf("%12.7f  %12.7f  %12.7f\n", b[0], b[1], b[2]);
}

double PWAtool::get_sfact() const
{
   return sfact;
}
double PWAtool::get_Re_rho(int m1, int m2) const
{
   if (m1 == 0 && m2 == 0)
      return rho00;
   else if (m1 == 0 && m2 == 1)
      return Re_rho01();
   else if (m1 == 0 && m2 == -1)
      return Re_rho0M();
   else if (m1 == 1 && m2 == 0)
      return Re_rho10;
   else if (m1 == 1 && m2 == 1)
      return rho11;
   else if (m1 == 1 && m2 == -1)
      return rho1M();
   else if (m1 == -1 && m2 == 0)
      return Re_rhoM0();
   else if (m1 == -1 && m2 == 1)
      return rhoM1;
   else if (m1 == -1 && m2 == -1)
      return rhoMM();
   else
      throw std::range_error("wave index out of range");
}

double PWAtool::get_Im_rho(int m1, int m2) const
{
   if (m1 == 0 && m2 == 0)
      return 0;
   else if (m1 == 0 && m2 == 1)
      return Im_rho01();
   else if (m1 == 0 && m2 == -1)
      return Im_rho0M();
   else if (m1 == 1 && m2 == 0)
      return Im_rho10;
   else if (m1 == 1 && m2 == 1)
      return 0;
   else if (m1 == 1 && m2 == -1)
      return 0;
   else if (m1 == -1 && m2 == 0)
      return Im_rhoM0();
   else if (m1 == -1 && m2 == 1)
      return 0;
   else if (m1 == -1 && m2 == -1)
      return 0;
   else
      throw std::range_error("wave index out of range");
}

void PWAtool::print_rho() const
{
   std::cout << "sdm11=" << rho11
             << ", sdm00=" << rho00
             << ", sdm10=" << Re_rho10
             << ", sdmM1=" << rhoM1
             << std::endl;
}

double PWAtool::get_real_Nevents(int powgt) const
{
   // return the number of real events in the real data buffer,
   // weighted by the factor pow(weight[i], powgt)

   if (powgt == 0)
      return real_events.size();
   long long int nreal = real_events.size();
   long long int nreal0 = nreal * Nslice / Nslices;
   long long int nreal1 = nreal * (Nslice + 1) / Nslices;
   nreal0 = (nreal0 < 0)? 0 : (nreal0 > nreal)? nreal : nreal0;
   nreal1 = (nreal1 < 0)? 0 : (nreal1 > nreal)? nreal : nreal1;
   nreal0 += Nstride_start;
   double sum=0;
   for (long long int i=nreal0; i < nreal1; i += Nstride_step) {
      double weight = real_events[i].weight;
      if (weight != 0) 
         sum += pow(weight, powgt);
   }
   return sum;
}

const PWAtool::event_T &PWAtool::get_real_event(unsigned int i) const
{
   if (i < real_events.size())
      return real_events[i];
   else
      throw std::range_error("event index out of range");
}

double PWAtool::get_mc_Nevents(int powgt) const
{
   // return the number of MC events in the MC data buffer,
   // weighted by the factor pow(weight[i], powgt)

   if (powgt == 0)
      return mc_events.size();
   long long int nmc = mc_events.size();
   long long int nmc0 = nmc * Nslice / Nslices;
   long long int nmc1 = nmc * (Nslice + 1) / Nslices;
   nmc0 = (nmc0 < 0)? 0 : (nmc0 > nmc)? nmc : nmc0;
   nmc1 = (nmc1 < 0)? 0 : (nmc1 > nmc)? nmc : nmc1;
   nmc0 += Nstride_start;
   double sum=0;
   //for (long long int i=nmc0; i < nmc1; i += Nstride_step) {
   for (long long int i=0; i < nmc; i += 1) {
      double weight = mc_events[i].weight;
      if (weight != 0) 
         sum += pow(weight, powgt);
   }
   return sum;
}

const PWAtool::event_T &PWAtool::get_mc_event(unsigned int i) const
{
   if (i < mc_events.size())
      return mc_events[i];
   else
      throw std::range_error("event index out of range");
}

int PWAtool::load_real_events(TChain *data,
                              double tabsmin, double tabsmax,
                              double Emin, double Emax)
{
   // read in real event information from a root tree
   // and store it in memory in the real_events vector.
 
   treerow_T row;
   data->SetBranchAddress("pwa_ntuple", (double*)&row);
   int nevents = 0;
   long long int nrows = data->GetEntries();
   for (long long int r = 0; r < nrows; ++r) {
      event_T event;
      data->GetEntry(r);
      if (row.tabs < tabsmin || row.tabs > tabsmax ||
          row.Ebeam < Emin || row.Ebeam > Emax)
      {  
          continue;
      }
      event.costheta = row.helcosthe;
      event.phi = row.helphi;
      event.weight = row.wgt;
      if (event.weight != 0) {
         real_events.push_back(event);
         ++nevents;
      }
   }
   return nevents;
}

int PWAtool::load_mc_events(TChain *mc,
                            double tabsmin, double tabsmax,
                            double Emin, double Emax)
{
   // read in mc event information from a root tree
   // and store it in memory in the mc_events vector.

   mctreerow_T row;
   mc->SetBranchAddress("pwa_ntuple", (double*)&row);
   int nevents = 0;
   long long int nrows = mc->GetEntries();
   for (long long int r = 0; r < nrows; ++r) {
      event_T event;
      mc->GetEntry(r);
      if (row.tabs < tabsmin || row.tabs > tabsmax ||
          row.Ebeam < Emin || row.Ebeam > Emax)
      {
         continue;
      }
      event.costheta = row.helcosthe;
      event.phi = row.helphi;
      event.weight = row.wgt;
      if (event.weight != 0) {
         mc_events.push_back(event);
         ++nevents;
      }
   }
   return nevents;
}

void PWAtool::gen_real_events(int nevents)
{
   // generate random events and store them in the real_events vector

   for (int i=0; i < nevents; ++i) {
      event_T event;
      event.costheta = randoms->Uniform(-1, 1);
      event.phi = randoms->Uniform(-M_PI, M_PI);
      event.weight = 1;
      real_events.push_back(event);
   }
}

void PWAtool::gen_mc_events(int nevents)
{
   // generate random events and store them in the mc_events vector

   for (int i=0; i < nevents; ++i) {
      event_T event;
      event.costheta = randoms->Uniform(-1, 1);
      event.phi = randoms->Uniform(-M_PI, M_PI);
      event.weight = 1;
      mc_events.push_back(event);
   }
}

// configure class functor to serve as TMinuit fcn //
double PWAtool::operator() (int &npar, double grad[4], int iflag)
{
   // evaluate the negative-log-likelihood function on the loaded data set,
   // where flag=1 : do initialization to prepare for a fit
   //       flag=2 : return the gradient in parameter space
   //       flag=3 : do any necessary wrap-up steps after a fit
   //       flag=4 : evaluate nll only (also implied for 1-3)

   double nll = 0;
   double gradnll[4] = {};
   long long int nreal = real_events.size();
   long long int nreal0 = nreal * Nslice / Nslices;
   long long int nreal1 = nreal * (Nslice + 1) / Nslices;
   nreal0 = (nreal0 < 0)? 0 : (nreal0 > nreal)? nreal : nreal0;
   nreal1 = (nreal1 < 0)? 0 : (nreal1 > nreal)? nreal : nreal1;
   const int Nmin=5;
   double Wmin[Nmin];
   const int Nmax=5;
   double Wmax[Nmax];
   for (int n=0; n < Nmin; ++n)
      Wmin[n] = 1e99;
   for (int n=0; n < Nmax; ++n)
      Wmax[n] = 0;
   nreal0 += Nstride_start;
   for (long long int i=nreal0; i < nreal1; i += Nstride_step) {
      angles_T Omega = {real_events[i].costheta,
                        real_events[i].phi};
      if (Omega.costheta < costheta_cut[0] ||
          Omega.costheta > costheta_cut[1] ||
          Omega.phi < phi_cut[0] ||
          Omega.phi > phi_cut[1])
      {
         continue;
      }
      double weight = real_events[i].weight;
      if (weight_2dhist) {
         int bin = weight_2dhist->FindBin(Omega.costheta, Omega.phi);
         weight *= weight_2dhist->GetBinContent(bin);
         if (weight == 0)
            continue;
      }
      double w = W(Omega);
      w = sqrt(w*w + 1e-6); // smooth out the minimum
      if (w < 0) {
         printf("error - real event %lld has negative W %f,"
                " invalidating result!\n", i, w);
         nll = sqrt(-1);
         break;
      }
      else {
         for (int n=0; n < Nmin; ++n) {
            if (w < Wmin[n]) {
               Wmin[n] = w;
               break;
            }
         }
         for (int n=0; n < Nmax; ++n) {
            if (w > Wmax[n]) {
               Wmax[n] = w;
               break;
            }
         }
      }
      nll -= log(w) * weight;
      gradnll[0] -= (Wprime[0] / w) * weight;
      gradnll[1] -= (Wprime[1] / w) * weight;
      gradnll[2] -= (Wprime[2] / w) * weight;
      gradnll[3] -= (Wprime[3] / w) * weight;
   }
   if (VERBOSE > 3) {
      printf(" sum over real events gives %f\n", nll);
   }
   if (VERBOSE > 5) {
      printf("%d lowest real W values are", Nmin);
      for (int n=0; n < Nmin; ++n)
         printf(" %e", Wmin[n]);
      printf("\n");
      printf("%d highest real W values are", Nmax);
      for (int n=0; n < Nmax; ++n)
         printf(" %e", Wmax[n]);
      printf("\n");
   }

   long long int nmc = mc_events.size();
   long long int nmc0 = nmc * Nslice / Nslices;
   long long int nmc1 = nmc * (Nslice + 1) / Nslices;
   nmc0 = (nmc0 < 0)? 0 : (nmc0 > nmc)? nmc : nmc0;
   nmc1 = (nmc1 < 0)? 0 : (nmc1 > nmc)? nmc : nmc1;
   for (int n=0; n < Nmin; ++n)
      Wmin[n] = 1e99;
   for (int n=0; n < Nmax; ++n)
      Wmax[n] = 0;
   nmc0 += Nstride_start;
   //for (long long int i=nmc0; i < nmc1; i += Nstride_step) {
   for (long long int i=0; i < nmc; i += 1) {
      angles_T Omega = {mc_events[i].costheta,
                        mc_events[i].phi};
      if (Omega.costheta < costheta_cut[0] ||
          Omega.costheta > costheta_cut[1] ||
          Omega.phi < phi_cut[0] ||
          Omega.phi > phi_cut[1])
      {
         continue;
      }
      double weight = mc_events[i].weight;
      if (weight_2dhist) {
         int bin = weight_2dhist->FindBin(Omega.costheta, Omega.phi);
         weight *= weight_2dhist->GetBinContent(bin);
         if (weight == 0)
            continue;
      }
      double w = W(Omega);
      if (w < 0) {
         printf("error - mc event %lld has negative W %f,"
               " invalidating result!\n", i,w);
         nll = sqrt(-1);
         break;
      }
      else {
         for (int n=0; n < Nmin; ++n) {
            if (w < Wmin[n]) {
               Wmin[n] = w;
               break;
            }
         }
         for (int n=0; n < Nmax; ++n) {
            if (w > Wmax[n]) {
               Wmax[n] = w;
               break;
            }
         }
      }
      nll += w * weight;
      gradnll[0] += Wprime[0] * weight;
      gradnll[1] += Wprime[1] * weight;
      gradnll[2] += Wprime[2] * weight;
      gradnll[3] += Wprime[3] * weight;
   }
   if (VERBOSE > 2) {
      printf(" sum over all terms gives %f\n", nll);
   }
   if (VERBOSE > 5) {
      printf("%d lowest mc W values are", Nmin);
      for (int n=0; n < Nmin; ++n)
         printf(" %e", Wmin[n]);
      printf("\n");
      printf("%d highest mc W values are", Nmin);
      for (int n=0; n < Nmax; ++n)
         printf(" %e", Wmax[n]);
      printf("\n");
   }

   if (iflag == 2 || VERBOSE > 0) {
      for (int par=0; par < 4; ++par) {
         grad[par] = 0;
         for (int comp=0; comp < 4; ++comp) {
            grad[par] += gradnll[comp] * PartialD[comp][par];
         }
         grad[par] *= sfact;
      }
   }
   return nll;
}

void PWAtool::FCN(int &npar,double grad[4], double &nll,
                  double param[4], int iflag)
{
    //play the role of the FCN function for the Minuit fit package

   if (VERBOSE > 0) {
      printf(" param is %12.9f,%12.9f,%12.9f, %12.9f\n", 
             param[0],param[1],param[2],param[3]);
      printf(" gradient is  %f,%f,%f,%f\n", 
             grad[0],grad[1],grad[2],grad[3]);
   }
#if FIT_PARAMETERIZATION_13
   owner->set_parameters13(param);
#else
   owner->set_parameters04(param);
#endif
   nll = (*owner)(npar, grad, iflag);
   static int nentries = 0;
   printf("fcn entry %d, iflag = %d, nll = %f\n", ++nentries, iflag, nll);
   if (VERBOSE > 0) {
       printf("rho11, rho00, rhoM1, Re_rho10, sfact =  %f, %f, %f, %f, %f\n",
              owner->rho11, owner->rho00, owner->rhoM1, owner->Re_rho10,
              owner->sfact);
   }
}

void PWAtool::SetFCN(TMinuit &fitter)
{
   // call TMinuit::SetFCN to associate this functor with the fitter

   fitter.SetFCN(&PWAtool::FCN);
}

TH1D *PWAtool::hist_real_costheta(const char *name, double phimin,
                                                    double phimax)
{
   // create histogram of the real data in costheta

   TH1D *h = new TH1D(name, "real costheta distribution", 100, -1.0, 1.0);
   //TH1D *h = new TH1D(name, "real costheta distribution", 100, -0.7, 0.7);
   long long int nreal = real_events.size();
   long long int nreal0 = nreal * Nslice / Nslices;
   long long int nreal1 = nreal * (Nslice + 1) / Nslices;
   nreal0 = (nreal0 < 0)? 0 : (nreal0 > nreal)? nreal : nreal0;
   nreal1 = (nreal1 < 0)? 0 : (nreal1 > nreal)? nreal : nreal1;
   nreal0 += Nstride_start;
   for (long long int i=nreal0; i < nreal1; i += Nstride_step) {
      if (real_events[i].phi >= phimin && real_events[i].phi <= phimax)
         h->Fill(real_events[i].costheta, real_events[i].weight);
   }
   return h;
}

TH1D *PWAtool::hist_real_phi(const char *name, double costmin,
                                               double costmax)
{
   // create histogram of the real data in phi

   TH1D *h = new TH1D(name, "real phi distribution", 100, -M_PI, M_PI);
   long long int nreal = real_events.size();
   long long int nreal0 = nreal * Nslice / Nslices;
   long long int nreal1 = nreal * (Nslice + 1) / Nslices;
   nreal0 = (nreal0 < 0)? 0 : (nreal0 > nreal)? nreal : nreal0;
   nreal1 = (nreal1 < 0)? 0 : (nreal1 > nreal)? nreal : nreal1;
   nreal0 += Nstride_start;
   for (long long int i=nreal0; i < nreal1; i += Nstride_step) {
      if (real_events[i].costheta >= costmin && 
          real_events[i].costheta <= costmax)
      {
         h->Fill(real_events[i].phi, real_events[i].weight);
      }
   }
   return h;
}

TH1D *PWAtool::hist_mc_costheta(const char *name, int weighted, double phimin,
                                                                double phimax)
{
   // create histogram of the mc data in costheta
   // and weight it with the fitted W pdf if weighted > 0.

   int print1=0;
   TH1D *h = new TH1D(name, "mc costheta distribution", 100, -1.0, 1.0);
   long long int nmc = mc_events.size();
   long long int nmc0 = nmc * Nslice / Nslices;
   long long int nmc1 = nmc * (Nslice + 1) / Nslices;
   nmc0 = (nmc0 < 0)? 0 : (nmc0 > nmc)? nmc : nmc0;
   nmc1 = (nmc1 < 0)? 0 : (nmc1 > nmc)? nmc : nmc1;
   nmc0 += Nstride_start;
   //for (long long int i=nmc0; i < nmc1; i += Nstride_step) {
   for (long long int i=0; i < nmc; i += 1) {
      if (mc_events[i].phi >= phimin && mc_events[i].phi <= phimax) {
         double weight = mc_events[i].weight;
         if (weighted == 1) {
            angles_T Omega = {mc_events[i].costheta, mc_events[i].phi};
            weight *= W(Omega);
         }
         h->Fill(mc_events[i].costheta, weight);
         if (print1++ == 0)
            print_rho();
      }
   }
   return h;
}

TH1D *PWAtool::hist_mc_phi(const char *name, int weighted, double costmin,
                                                           double costmax)
{
   // create histogram of the mc data in phi
   // and weight it with the fitted W pdf if weighted > 0.

   int print1=0;
   TH1D *h = new TH1D(name, "mc phi distribution", 100, -M_PI, M_PI);
   long long int nmc = mc_events.size();
   long long int nmc0 = nmc * Nslice / Nslices;
   long long int nmc1 = nmc * (Nslice + 1) / Nslices;
   nmc0 = (nmc0 < 0)? 0 : (nmc0 > nmc)? nmc : nmc0;
   nmc1 = (nmc1 < 0)? 0 : (nmc1 > nmc)? nmc : nmc1;
   nmc0 += Nstride_start;
   //for (long long int i=nmc0; i < nmc1; i += Nstride_step) {
   for (long long int i=0; i < nmc; i += 1) {
      if (mc_events[i].costheta >= costmin && mc_events[i].costheta <= costmax) 
      {
         double weight = mc_events[i].weight;
         if (weighted == 1 ) {
            angles_T Omega = {mc_events[i].costheta, mc_events[i].phi};
            weight *= W(Omega);
         }
         h->Fill(mc_events[i].phi, weight);
         if (print1++ == 0)
            print_rho();
      }
   }
   return h;
}

TH2D *PWAtool::hist_real_weight(const char *name, int power)
{
   // create histogram of the event weights ** power
   // and fill a 2D plot in costheta, phi for real data

   std::stringstream title;
   title << "real event weights ** " << power;
   TH2D *h = new TH2D(name, title.str().c_str(), 30, -1, 1, 30, -M_PI, M_PI);
   long long int nreal = real_events.size();
   long long int nreal0 = nreal * Nslice / Nslices;
   long long int nreal1 = nreal * (Nslice + 1) / Nslices;
   nreal0 = (nreal0 < 0)? 0 : (nreal0 > nreal)? nreal : nreal0;
   nreal1 = (nreal1 < 0)? 0 : (nreal1 > nreal)? nreal : nreal1;
   nreal0 += Nstride_start;
   for (long long int i=nreal0; i < nreal1; i += Nstride_step) {
      double wgt = real_events[i].weight;
      h->Fill(real_events[i].costheta, real_events[i].phi, 
              (power != 0)? pow(wgt, power) : 1);
   }
   h->SetStats(0);
   return h;
}

TH2D *PWAtool::hist_mc_weight(const char *name, int power)
{
   // create histogram of the event weights ** power
   // and fill a 2D plot in costheta, phi for monte carlo

   std::stringstream title;
   title << "mc event weights ** " << power;
   TH2D *h = new TH2D(name, title.str().c_str(), 30, -1, 1, 30, -M_PI, M_PI);
   long long int nmc = mc_events.size();
   long long int nmc0 = nmc * Nslice / Nslices;
   long long int nmc1 = nmc * (Nslice + 1) / Nslices;
   nmc0 = (nmc0 < 0)? 0 : (nmc0 > nmc)? nmc : nmc0;
   nmc1 = (nmc1 < 0)? 0 : (nmc1 > nmc)? nmc : nmc1;
   nmc0 += Nstride_start;
   for (long long int i=nmc0; i < nmc1; i += Nstride_step) {
      double wgt = mc_events[i].weight;
      h->Fill(mc_events[i].costheta, mc_events[i].phi, 
              (power != 0)? pow(wgt, power) : 1);
   }
   h->SetStats(0);
   return h;
}

double PWAtool::tf1_W(double *var, double *par)
{  
   rho11 = par[0];
   rho00 = par[1];
   rhoM1 = par[2];
   Re_rho10 = par[3];
   angles_T Omega = {var[0],var[1]};
   return  W(Omega);
}

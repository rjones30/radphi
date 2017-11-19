#!/usr/bin/python
#
# PWAomega.py - wrapper script about the c++ class PWAtool
#               for doing PWA fits to the decay angular distributions
#               of the omega->gamma,pi0 reaction in the Radphi data.
#
# authors: richard.t.jones at uconn.edu, fridah.mokaya at uconn.edu
# version: april 5, 2016

import os
import glob
import sys
import numpy
import math
import random
import re

#import mc
#import data

from ROOT import *
gSystem.Load("PWAtool_cc.so")

initial_value = random.random()
#initial_value = 0.3
initial_stepsize = 0.03

fitter = 0

def init(check=0, tlim=[0, 99],Elim=[0,99]):
   """
   Load real events and Monte Carlo events from their respective
   root trees, and define the fit parameters in preparation for
   a call to fit(). Setting the first argument to a non-zero value
   forces Minuit to check the computation of the local gradient of
   the negative-log-likelihood function against its approximation by
   finite differences. The tlim argument provides a way for the user
   to request a restricted range in kinematic t in the selection of
   events for the fit. Normally this function would be called before
   each invocation of fit().
   """
   realdata_dir = "/pnfs4/phys.uconn.edu/data/Gluex/radphi/data_pwa-7-2016/"
   realdata = TChain("pwa_tree", "real omega ntuple")
   for rfile in range(0,61):
       realdata.Add(realdata_dir + "pwa_tree_%d.root" %(rfile))
   mcdata_dir = "/pnfs4/phys.uconn.edu/data/Gluex/radphi/mc_pwa-7-2016/"
   mcdata = TChain("pwa_tree", "mc omega ntuple")
   for mfile in range(0,61):
       mcdata.Add(mcdata_dir + "pwa_tree_%d.root" %(mfile))

   print realdata.GetEntries(), "real events found"
   print mcdata.GetEntries(), "mc events found"
   print initial_value  , "initial value"
   global tool
   tool = PWAtool()
   tool.load_real_events(realdata, tlim[0], tlim[1],Elim[0],Elim[1])
   print "real events kept:", tool.get_real_Nevents()
   tool.load_mc_events(mcdata, tlim[0], tlim[1],Elim[0],Elim[1])
   print "mc events kept:", tool.get_mc_Nevents()
   
   global fitter
   if not fitter:
      fitter = TMinuit(4)
      tool.SetFCN(fitter)
      err1 = Long()
      fitter.mnparm(0, "a0", initial_value, initial_stepsize, 0, 0, err1)
      err2 = Long()
      fitter.mnparm(1, "c0", initial_value, initial_stepsize, 0, 0, err2)
      err3 = Long()
      fitter.mnparm(2, "a1", initial_value, initial_stepsize, 0, 0, err3)
      err4 = Long()
      fitter.mnparm(3, "c2", initial_value, initial_stepsize, 0, 0, err4)
      if err1 > 0 or err2 > 0 or err3 > 0 or err4 > 0:
         print "fitter.mnparm choked, cannot continue!"
         sys.exit(1)
      if check:
         fitter.Command("set gradient")
      else:
         fitter.Command("set gradient 1")
      fitter.Command("set errdef 0.5")

def init_canvas():
   """
   Several of the PWAomega methods require an open graphics window for
   plotting results. This function creates a TCanvas for this purpose
   with a set of dimensions that seem appropriate for viewing the
   results generated by the PWAomega pacakge.
   """
   global c1
   try:
      if c1:
         return
   except NameError:
      c1 = TCanvas("c1", "c1", 5, 10, 550, 500)
      c1.SetLeftMargin(0.12)
      c1.SetRightMargin(0.95)

def init_randoms():
   """
   Initialize the uniform random number generator for this pacakge
   """
   global rgen
   try:
      if rgen:
         pass
   except NameError:
      rgen = TRandom3();

def gen_weight_2dhist(smoothing=3, acceptance_threshold=0.1):
   """
   The weighted maximum likelihood method for unbinned fits to large
   multidimensional data sets requires a multidimensional function to
   be defined that characterizes the average behavior of the weight
   factors over the domain of the fit function. In the case of PWAomega
   fits, the fit function is defined on the 2d axes of (costheta,phi) 
   and the weight_2dhist function is the value of <w>/<w^2> over this
   finite 2d space, where w is the weight of real events that fall in
   the vicinity of the given point in (costheta,phi). If the weights
   are properly defined, this function should be posiitve-semidefinite
   over this full domain. This method uses statistical sampling over
   the set of real events to estimate this function as the contents
   of a TH2D histogram, which is passed back as the return value.
   The smoothing argument gives the degree of smoothing that should
   be applied to the histogram of <w> before it is divided by <w^2>.
   The acceptance_threshold is the minimum value of the Monte Carlo
   acceptance, as a fraction of the maximum value over the fit domain,
   below which the value of <w> is overwritten with zero.
   """
   global realw_1
   realw_1 = tool.hist_real_weight("realw_1", 1)
   for n in range(0, smoothing):
      realw_1.Smooth()

   global realw_2
   realw_2 = tool.hist_real_weight("realw_2", 2)
   realw_2.Smooth()

   global mcw_1
   mcw_1 = tool.hist_mc_weight("mc2_1", 1)
   mcw_1.Smooth()

   global weight_2dhist
   weight_2dhist = TH2D(realw_1)
   weight_2dhist.Divide(realw_2)
   min_acceptance = mcw_1.GetMaximum() * acceptance_threshold
   for ix in range(0, mcw_1.GetNbinsX()):
      for iy in range(0, mcw_1.GetNbinsY()):
         if weight_2dhist.GetBinContent(ix+1, iy+1) < 0 or \
                    mcw_1.GetBinContent(ix+1, iy+1) < min_acceptance:
            weight_2dhist.SetBinContent(ix+1, iy+1, 0)
   return weight_2dhist

def fit(maxcalls=1000, interactive=1):
   """
   Perform a PWA fit using the MIGRAD facility of Minuit. At the end
   of the first pass the user is given a chance to enter follow-up
   commands to refine the search, or simply hit enter to accept it.
   To suppress the interactive step, set argument interactive to 0.
   """
   global fitter
   try:
      if fitter:
         pass
   except NameError:
      print "PWAomega.fit error - must call init first!"
      return
   fitter.Command("set strategy  2")
   fitter.Command("set errordef  0.5")
   fitter.Command("migrad " + str(maxcalls))
   """
   fitter.Command("hesse " + str(maxcalls))
   fitter.Command("minos " + str(maxcalls))
   """
   while interactive:
      resp = raw_input("Enter a Minuit command, s to swap parameters," +
                       "or return to accept: ")
      if resp:
         if resp == 's':
            swap_parameters04()
         else:
            fitter.Command(resp)
      else:
         break

def multifit(nfits, maxcalls=1000):
   """
   Perform a series of PWA fits using the MIGRAD facility of Minuit,
   starting from a new set of random initial parameters for each fit.
   Each fit consists of a sequence of "migrad" steps which repeat
   until the fit converges. If convergence is not reached within 5
   "migrad" commands, an "improve" is issued, followed by a new
   cycle of "migrad"s.
   """
   global fitter
   try:
      if fitter:
         pass
   except NameError:
      print "PWAomega.multifit error - must call init first!"
      return

   init_randoms()
   fitter.Command("set strategy  2")
   fitter.Command("set errordef  0.5")
   for fit in range(0, nfits):
      print ">>>>> Starting fit {0} <<<<<".format(fit)
      a0 = rgen.Uniform(2.0) - 1.0
      c0 = rgen.Uniform(2.0) - 1.0
      a1 = rgen.Uniform(2.0) - 1.0
      c2 = rgen.Uniform(2,0) - 1.0
      fitter.Command("set parameter 1 {0:20.12f}".format(a0))
      fitter.Command("set parameter 2 {0:20.12f}".format(c0))
      fitter.Command("set parameter 3 {0:20.12f}".format(a1))
      fitter.Command("set parameter 4 {0:20.12f}".format(c2))
      fitter.Command("migrad " + str(maxcalls))
      count = 1
      while fitter.fCstatu != "CONVERGED ":
         if count == 5:
            fitter.Command("improve")
            count = 0
         fitter.Command("migrad " + str(maxcalls))
         count += 1
      print ">>>>> Completed fit {0} <<<<<".format(fit)
   return nfits

def fit30t(maxcalls=1000):
   """
   Perform a series of fits over the 30 separate divisions of the data
   according to |t|. The |t| divisions are hard-wired in the code below.
   """
   for sample in range(0,30):
      print "Now fitting |t| range", str(tbins[sample][0]) + ",", \
                                     str(tbins[sample][1]) + " GeV^2"
      global initial_stepsize
      """
      initial_stepsize = raw_input()
      resp = raw_input("starting step size? [" + str (initial_stepsize) + "] ")
      try:
         resp = float(resp)
      except ValueError:
         resp = initial_stepsize
      initial_stepsize = resp
      """
      global fitter
      fitter = 0
      init(0, tbins[sample])
      fit(maxcalls, 0)
      report_fit()

def report_fit(sample):
   """
   Print a report showing the best-fit parameters for the SDM and
   total acceptance-corrected yield, together with their errors.
   The sample argument should be one of the key values from the
   lookup tables tbins, Ebins, Nmcgen, etc. which it uses to look
   up the weight factors for renormalizing the errors.
   """
   a1 = Double()
   c1 = Double()
   a2 = Double()
   c3 = Double()
   a1err = Double()
   c1err = Double()
   a2err = Double()
   c3err = Double()
   fitter.GetParameter(0, a1, a1err)
   fitter.GetParameter(1, c1, c1err)
   fitter.GetParameter(2, a2, a2err)
   fitter.GetParameter(3, c3, c3err)
   param = numpy.array([a1, c1, a2, c3])
   grad = numpy.array([0., 0., 0., 0.])
   nll_min = Double()
   npars = Long(4)
   fitter.Eval(npars, grad, nll_min, param, 4)
   rho11 = tool.get_Re_rho(1,1)
   rho00 = tool.get_Re_rho(0,0)
   rhoM1 = tool.get_Re_rho(-1,1)
   rho10 = tool.get_Re_rho(1,0)
   traceRho = 2*rho11 + rho00
   rho11 /= traceRho
   rho00 /= traceRho
   rhoM1 /= traceRho
   rho10 /= traceRho
   Wintegral = 8 * math.pi / 3
   dt = tbins[sample][1] - tbins[sample][0]
   dsdt = traceRho * Wintegral / dt
   dsdt *= Nmcgen[sample] / (4 * math.pi * integrated_luminosity)
   dsdt /= decay_branching_ratio
   dsdt /= protons_per_atom

   # Now we need to propagate the error matrix from the (a0,c1,a2,c3) basis
   # to the (rho11,rho00,rhoM1,rho10) basis, and from there to the basis of
   # dsigma/dt and (nrho11, nrho00, nrhoM1, nrho10) normalized to Tr{nrho}=1.
   emat = numpy.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], dtype=numpy.float64)
   fitter.mnemat(emat, 4)
   rhocovar = TMatrixD(4, 4)
   n = 0
   for row in range(0,4):
      for col in range(0,4):
         rhocovar[row][col] = emat[n]
         n += 1
   D = tool.get_PartialD()
   DT = TMatrixD(4,4)
   DT.Transpose(D)
   Drhocovar = TMatrixD(4,4)
   Drhocovar.Mult(D, rhocovar)
   rhocovar.Mult(Drhocovar, DT)
   N = TMatrixD(4,4)
   N[0][0] = 2
   N[0][1] = 1
   N[0][2] = 0
   N[0][3] = 0
   N[1][0] = -2 * rho00 / traceRho
   N[1][1] = 2 * rho11 / traceRho
   N[1][2] = 0
   N[1][3] = 0
   N[2][0] = -2 * rhoM1 / traceRho
   N[2][1] = -rhoM1 / traceRho
   N[2][2] = 1 / traceRho
   N[2][3] = 0
   N[3][0] = -2 * rho10 / traceRho
   N[3][1] = -rho10 / traceRho
   N[3][2] = 0
   N[3][3] = 1 / traceRho
   NT = TMatrixD(4,4)
   NT.Transpose(N)
   Nrhocovar = TMatrixD(4,4)
   Nrhocovar.Mult(N, rhocovar)
   rhocovar.Mult(Nrhocovar, NT)
   traceVar = rhocar[0][0]
   dsdt_error = math.sqrt(traceVar) * Wintegral / dt
   dsdt_error *= Nmcgen[sample] / (4 * math.pi * integrated_luminosity)
   dsdt_error /= decay_branching_ratio
   dsdt_error /= protons_per_atom
   print "total acceptance-corrected differential cross section =", \
         "{0:8.5f} +/- {1:8.5f} ubarns/GeV^2".format(dsdt, dsdt_error)
   print "   rho11 = {0:8.5f} +/- {1:8.5f}".format(rho11, 
                                      0.5 * math.sqrt(rhocovar[1][1]))
   print "   rho00 = {0:8.5f} +/- {1:8.5f}".format(rho00, 
                                            math.sqrt(rhocovar[1][1]))
   print "   rhoM1 = {0:8.5f} +/- {1:8.5f}".format(rhoM1,
                                            math.sqrt(rhocovar[2][2]))
   print "Re_rho10 = {0:8.5f} +/- {1:8.5f}".format(rho10,
                                            math.sqrt(rhocovar[3][3]))
   print
   print "Final fit status was", fitter.fCstatu
   print "Final nll_min was", nll_min
   print "Final correlation matrix is:"
   for i in range(0,4):
      print "{0:8.3f}  {1:8.3f}  {2:8.3f}  {3:8.3f}".format(
            rhocovar[i][0] / math.sqrt(rhocovar[i][i] * rhocovar[0][0]),
            rhocovar[i][1] / math.sqrt(rhocovar[i][i] * rhocovar[1][1]),
            rhocovar[i][2] / math.sqrt(rhocovar[i][i] * rhocovar[2][2]),
            rhocovar[i][3] / math.sqrt(rhocovar[i][i] * rhocovar[3][3]))

def plot_fit_results(logfile, parname):
   """
   Results from a sequence of fits to subsets of a given event
   sample have been saved to logfile, in the form of output from 
   report_fit(). Now scan this logfile and plot the values on a
   single graph and fit them to a constant level to see if they are
   all consistent with each other.
   """
   parvalue = []
   parerror = []
   for line in open(logfile):
      match = re.match(r"^ *" + parname + " *= *([-.0-9]+) \+/- *([.0-9]+)",
                       line)
      if match:
         parvalue.append(float(match.group(1)))
         parerror.append(float(match.group(2)))
   nvals = len(parvalue)
   graph = TH1D(parname, parname, nvals, 0, nvals)
   for n in range(0, nvals):
      graph.SetBinContent(n+1, parvalue[n])
      graph.SetBinError(n+1, parerror[n])
   graph.Fit("pol0")
   graph.SetStats(0)
   gStyle.SetOptFit(1)
   init_canvas()
   graph.Draw()
   return graph

def swap_parameters04():
   """
   There is a dualism in the rho parameterization that is used in this 
   PWA, such that two very different parameter sets give rise to the
   same rho matrix and hence the same NLL value. Sometimes one of these
   two sets is easier to work with than the other because all of the
   parameters are of the same scale. This function fetches the current
   parameters from the fitter, computes and uploads the dual parameter
   set to the fitter, and checks that the NLL is the unchanged.
   """
   a1 = Double()
   c1 = Double()
   a2 = Double()
   c3 = Double()
   err = Double()
   fitter.GetParameter(0, a1, err)
   fitter.GetParameter(1, c1, err)
   fitter.GetParameter(2, a2, err)
   fitter.GetParameter(3, c3, err)
   param = numpy.array([a1, c1, a2, c3])
   grad = numpy.array([0., 0., 0., 0.])
   nll = Double()
   npars = Long(4)
   fitter.Eval(npars, grad, nll, param, 4)
   a_1 = c1 * c3
   c_1 = -2 * a1 * c3
   a_2 = a2
   c_3 = -1 / (2 * c3)
   fitter.Command("set parameter 1 {0:20.12f}".format(a_1))
   fitter.Command("set parameter 2 {0:20.12f}".format(c_1))
   fitter.Command("set parameter 3 {0:20.12f}".format(a_2))
   fitter.Command("set parameter 4 {0:20.12f}".format(c_3))
   nllnew = Double()
   fitter.Eval(npars, grad, nllnew, param, 4)
   print ""
   print "swapping param = {0:12.8f}, {1:12.8f}, {2:12.8f}, {3:12.8f}".format(
                              a1,        c1,        a2,        c3   )
   print "  for dual set = {0:12.8f}, {1:12.8f}, {2:12.8f}, {3:12.8f}".format(
                              a_1,       c_1,       a_2,       c_3  )
   print "nll old, new = {0:12.7f}, {1:12.7f}".format(nll, nllnew)

def swap_parameters13():
   """
   Same as swap_parameters04 except that it uses the 1,3 
   parameterization instead of the original on that we later called 0,4.
   """
   a1 = Double()
   c1 = Double()
   a2 = Double()
   c3 = Double()
   err = Double()
   fitter.GetParameter(0, a1, err)
   fitter.GetParameter(1, c1, err)
   fitter.GetParameter(2, a2, err)
   fitter.GetParameter(3, c3, err)
   param = numpy.array([a1, c1, a2, c3])
   grad = numpy.array([0., 0., 0., 0.])
   nll = Double()
   npars = Long(4)
   fitter.Eval(npars, grad, nll, param, 4)
   a_1 = c1 * c3 * a1**2
   c_1 = -2 * a1 * c3 * a1 / a_1
   a_2 = a2 * a1 / a_1
   c_3 = -1 / (2 * c3 * a1 * a_1)
   fitter.Command("set parameter 1 {0:20.12f}".format(a_1))
   fitter.Command("set parameter 2 {0:20.12f}".format(c_1))
   fitter.Command("set parameter 3 {0:20.12f}".format(a_2))
   fitter.Command("set parameter 4 {0:20.12f}".format(c_3))
   nllnew = Double()
   param2 = numpy.array([a_1, c_1, a_2, c_3])
   fitter.Eval(npars, grad, nllnew, param2, 4)
   print ""
   print "swapping param = {0:12.8f}, {1:12.8f}, {2:12.8f}, {3:12.8f}".format(
                              a1,        c1,        a2,        c3   )
   print "  for dual set = {0:12.8f}, {1:12.8f}, {2:12.8f}, {3:12.8f}".format(
                              a_1,       c_1,       a_2,       c_3  )
   print "nll old, new = {0:12.7f}, {1:12.7f}".format(nll, nllnew)

def contour(ph, pv, nbins=10, nbox=9, nsigma=1):
   """
   Plot the nll function in the region +/- nsigma about the current
   fit value in parameter space as a 2D histogram with dimensions nbins
   x nbins, plotting parameter ph on the horizontal axis and pv on the
   vertical. The current minimum of the nll function is subtracted
   before it is plotted. The nbox argument asks for each bin in the
   histogram to be replaced with a square nbox x nbox smaller bins and
   filled using a linear extrapolation from the center square using the
   gradient evaluated there. It works best when nbox is an odd number,
   so the center bin is centered in the box. Parameter numbers in the
   ph,pv arguments are based at 1.
   """
   a1 = Double()
   c1 = Double()
   a2 = Double()
   c3 = Double()
   a1err = Double()
   c1err = Double()
   a2err = Double()
   c3err = Double()
   fitter.GetParameter(0, a1, a1err)
   fitter.GetParameter(1, c1, c1err)
   fitter.GetParameter(2, a2, a2err)
   fitter.GetParameter(3, c3, c3err)
   param = numpy.array([a1, c1, a2, c3])
   parerr = numpy.array([a1err, c1err, a2err, c3err])
   grad = numpy.array([0., 0., 0., 0.])
   nll_min = Double()
   npars = Long(4)
   fitter.Eval(npars, grad, nll_min, param, 4)
   xlim = (param[ph-1] - nsigma * parerr[ph-1],
           param[ph-1] + nsigma * parerr[ph-1])
   ylim = (param[pv-1] - nsigma * parerr[pv-1],
           param[pv-1] + nsigma * parerr[pv-1])
   h = gROOT.FindObject("contour")
   if h:
      h.Delete()
   h = TH2D("contour", "nll contour in parameter " + str(pv) +
                       " vs " + str(ph), 
            nbins * nbox, xlim[0], xlim[1], nbins * nbox, ylim[0], ylim[1])
   nll = Double()
   for ix in range(1, nbox * nbins + 1, nbox):
      param[ph-1] = h.GetXaxis().GetBinCenter(ix + nbox/2)
      for iy in range(1, nbox * nbins + 1, nbox):
         param[pv-1] = h.GetYaxis().GetBinCenter(iy + nbox/2)
         fitter.Eval(npars, grad, nll, param, 2)
         for bx in range (0, nbox):
            dx = (bx - nbox/2) * h.GetXaxis().GetBinWidth(ix)
            for by in range (0, nbox):
               dy = (by - nbox/2) * h.GetYaxis().GetBinWidth(iy)
               val = nll - nll_min + grad[ph-1]*dx + grad[pv-1]*dy
               h.SetBinContent(ix + bx, iy + by, val)
   h.SetContour(100)
   init_canvas()
   h.Draw("colz")
   return h

def compare_costheta(phicut=(-math.pi,math.pi)):
   """
   Plot the costheta distribution for the real data (blue error bars),
   and superimpose on top the Monte Carlo phase space distribution 
   (green error bars) renormalized to match the total number of real
   events, and the PWA fit (red error bars). The interval defined in
   the optional argument phicut selects a region in phi within which
   the events are selected for inclusion in this plot.
   """
   hreal = tool.hist_real_costheta("real_costheta", phicut[0], phicut[1])
   hmcraw = tool.hist_mc_costheta("mcraw_costheta", 0, phicut[0], phicut[1])
   hmcwgt = tool.hist_mc_costheta("mcwgt_costheta", 1, phicut[0], phicut[1])
   hreal.SetLineColor(kBlue)
   init_canvas()
   hreal.Draw()
   hmcraw.Scale(hreal.Integral() / hmcraw.Integral())
   hmcraw.SetLineColor(kGreen)
   hmcwgt.SetLineColor(kRed)
   hmcraw.Draw("e,same")
   hmcwgt.Draw("e,same")
   c1.Update

def compare_phi(costhetacut=(-1,1)):
   """
   Plot the phi distribution for the real data (blue error bars), and
   superimpose on top the Monte Carlo phase space distribution (green
   error bars) renormalized to match the total number of real events,
   and the PWA fit (red error bars).  The interval defined in the
   optional argument costhetacut selects a region in costheta within
   which the events are selected for inclusion in this plot.
   """
   hreal = tool.hist_real_phi("real_phi", costhetacut[0], costhetacut[1])
   hmcraw = tool.hist_mc_phi("mcraw_phi",0, costhetacut[0], costhetacut[1])
   hmcwgt = tool.hist_mc_phi("mcwgt_phi", 1, costhetacut[0], costhetacut[1])
   hreal.SetLineColor(kBlue)
   init_canvas()
   hreal.Draw()
   hmcraw.Scale(hreal.Integral() / hmcraw.Integral())
   hmcraw.SetLineColor(kGreen)
   hmcwgt.SetLineColor(kRed)
   hmcwgt.Draw("e,same")
   hmcraw.Draw("e,same")
   c1.Update

def test_rho_svd():
   """
   This function contains various mathematical experiments that I did to
   prove to myself that the parameterization I adopted for the SDM makes
   sense. None of this needs to be repeated for new fits, it was just a
   one-time exercise. In particular, I was curious why certain values 
   for the SDM matrix elements that lie outside the limits imposed by 
   the Schwarz inequality for off-diagonal elements (but not too far)
   seemed to never produce negative W values anywhere in phase space.
   This turns out to be true in the unpolarized case, but probably not
   true for the more general case with a polarized beam. Nevermind,
   because the SDM still must respect the Schwarz inequality on its own,
   regardless of whether violating it does or does not cause W to dip
   below zero at some point or other in the space of decay kinematics
   for this particular channel.
   """
   rho = TMatrixD(3,3)
   svd = TDecompSVD(rho)
   init_randoms()
   while True:
      rho11 = rgen.Uniform(0.5)
      rho00 = rgen.Uniform(1)
      scale10 = math.sqrt(rho11 * rho00) * 3
      rho10 = rgen.Uniform(-scale10, scale10)
      scaleM1 = rho11 * 2
      rhoM1 = rgen.Uniform(-scaleM1, scaleM1)
      rho[0][0] = rho11
      rho[1][0] = -rho10
      rho[2][0] = rhoM1
      rho[0][1] = -rho10
      rho[1][1] = rho00
      rho[2][1] = rho10
      rho[0][2] = rhoM1
      rho[1][2] = rho10
      rho[2][2] = rho11
      print "rho:" 
      for i in range(0,3):
         print "{0:12.7f} {1:12.7f} {2:12.7f}".format(
                 rho[i][0], rho[i][1], rho[i][2])
      tool.set_rho(rho11, rho00, rhoM1, rho10)
      grad = numpy.array([0,0,0,0], dtype=numpy.float64)
      nll = Double()
      npars = Long(4)
      nll = tool(npars, grad, 4)
      if math.isnan(nll):
         if abs(rho10) <= math.sqrt(rho00 * rho11) and abs(rhoM1) <= rho11:
            print "simple boundary check failed!"
         else:
            continue
      else:
         if rho10 * rho10 > rho00 * rho11 or rhoM1 * rhoM1 > rho11 * rho11:
            print "simple boundary check foiled!"
         else:
            continue
      svd.SetMatrix(rho)
      svd.Decompose()
      sval = svd.GetSig()
      U = svd.GetU()
      V = svd.GetV()
      print "svd of rho:"
      for sig in range(0, len(sval)):
         if U[0][sig] * V[0][sig] < 0:
            sval[sig] = -sval[sig]
            V[0][sig] = -V[0][sig]
            V[1][sig] = -V[1][sig]
            V[2][sig] = -V[2][sig]
      for sig in range(0, len(sval)):
         print "{0:12.7f}   {1:10.7f} {2:10.7f} {3:10.7f}" + \
               " {4:10.7f} {5:10.7f} {6:10.7f}".format(
                sval[sig], U[sig][0], U[sig][1], U[sig][2],
                V[sig][0], V[sig][1], V[sig][2])

      resp = raw_input("Press return to continue, e to evaluate: ")
      if len(resp) == 0 or resp[0] != 'e':
         continue
      param = numpy.array([0,0,0,0], dtype=numpy.float64)
      i1 = 0
      for sig in range(0,3):
         if U[1][sig] == 0:
            i1 = sig
            break
      i2 = (i1 + 1) % 3
      i0 = (i1 + 2) % 3
      param[0] = sval[i0] * U[0][i0]
      param[1] = sval[i0] * U[1][i0]
      param[2] = sval[i1] * U[0][i1]
      param[3] = sval[i2] * U[1][i2]
      #tool.set_rho(param)
      tool.set_rho(rho11, rho00, rhoM1, rho10)
      grad = numpy.array([0,0,0,0], dtype=numpy.float64)
      nll = Double()
      npars = Long(4)
      nll = tool(npars, grad, 4)
      print "nll = ", nll

def fake(nevents=10000000):
   """
   A simple random event generator that was useful for generating fake
   data that could be used for some of the mathematical experiments
   implement in test_rho_svd above. Who knows, it might come in handy
   in another context.
   """
   global tool
   try:
      if tool:
         pass
   except NameError:
      tool = PWAtool()
   tool.gen_real_events(nevents)
   tool.gen_mc_events(nevents)


# The following methods rely on the existence of a set of log files from
# jobs that carried out a series of fits to individual event samples,
# and stored the output generated by MINUIT. The files are all contained
# within a directory tree pointed to by the symbol (or directory name)
# fitresults. The specific structure of the files under that directory
# is detailed in the code below.

def loadfit(sample, nfit=-1):
   """
   Browse through the fitresults directory for the stated sample and 
   load the best-fit parameters for W into the PWAtool object. If the
   fit did not converge, return 0, otherwise return 1. If called with
   argument nfit < 0 (default) then it scans through all of the fit
   results for the requested sample and loads the one with the n'th
   smallest NLL, with nfit=-1 getting the smallest one, nfit=-2 the 
   next-smallest, etc. Return value is the number of fits that returned
   this value.
   """
   global tool
   try:
      if tool:
         pass
   except NameError:
      tool = PWAtool()
   fitrange = range(nfit, nfit+1)
   if nfit < 0:
      fitrange = range(0,150)
   a0 = 0
   c0 = 0
   a1 = 0
   c2 = 0
   a0err = 0
   c0err = 0
   a1err = 0
   c2err = 0
   NLLmin = []
   for fit in fitrange:
      fitlogfilename = "fitresults/pwafits{0}/PWAfits{1}.out".format(sample, fit)
      try:
         fitlog = open(fitlogfilename)
      except IOError:
         print "file", fitlogfilename, "not found, skipping"
         continue
      for line in fitlog:
         match = re.match(r"^ FCN=([0-9.+-Ee]+) +FROM MIGRAD +" +
                          r"STATUS=([A-Z]+) +([0-9]+) CALLS +([0-9]+) TOTAL",
                          line)
         if match and match.group(2) == "CONVERGED":
            break
      if match:
         NLL = float(match.group(1))
         for line in fitlog:
            match = re.match(r"^ *([1-4]) +([ac][012]) +([.0-9Ee+-]+)" +
                             r" +([.0-9Ee+-]+) +([.0-9Ee+-]+)" +
                             r" +([.0-9Ee+-]+)",
                             line)
            if match and match.group(2) == "a0":
               a0 = float(match.group(3))
               a0err = float(match.group(4))
            elif match and match.group(2) == "c0":
               c0 = float(match.group(3))
               c0err = float(match.group(4))
            elif match and match.group(2) == "a1":
               a1 = float(match.group(3))
               a1err = float(match.group(4))
            elif match and match.group(2) == "c2":
               c2 = float(match.group(3))
               c2err = float(match.group(4))
               break
         if not match:
            sys.stderr.write("extracting parameters from file " + 
                             fitlogfilename + " failed, aborting loadfit\n")
            return 0
         for row, line in zip(range(0,2), fitlog):
            pass
         errmat = numpy.array([[0.]*4 for row in range(0,4)])
         for row, line in zip(range(0,4), fitlog):
            values = line.strip().split()
            for col in range(0,4):
               errmat[row][col] = float(values[col])
         minimum = 0
         for minimum in range(0, len(NLLmin)):
            if abs(NLL - NLLmin[minimum]["value"]) < 0.5:
               NLLmin[minimum]["count"] += 1
               minimum = -1
               break
            elif NLL < NLLmin[minimum]["value"]:
               NLLmin.insert(minimum, {})
               break
         if minimum == -1:
            continue
         elif  minimum == len(NLLmin):
            NLLmin.append({})
         NLLmin[minimum]["value"] = NLL
         NLLmin[minimum]["count"] = 1
         NLLmin[minimum]["a0"] = a0
         NLLmin[minimum]["a0err"] = a0err
         NLLmin[minimum]["c0"] = c0
         NLLmin[minimum]["c0err"] = c0err
         NLLmin[minimum]["a1"] = a1
         NLLmin[minimum]["a1err"] = a1err
         NLLmin[minimum]["c2"] = c2
         NLLmin[minimum]["c2err"] = c2err
         NLLmin[minimum]["errmat"] = errmat
         print "new minimum found with NLL", NLL
   if len(NLLmin) == 0:
      sys.stderr.write("sample " + str(sample) + " fit " + str(fit) +
                       " never converged, no solution available," +
                       " aborting loadfit\n")
      return 0
   if nfit < 0 and abs(nfit) > len(NLLmin):
      sys.stderr.write("loadfit error - requested minimum number " + 
                       str(-nfit) + " but only " + str(len(NLLmin)) +
                       " minima were found, aborting loadfit\n")
      return 0
   nfit = abs(nfit) - 1
   a0 = NLLmin[nfit]["a0"]
   c0 = NLLmin[nfit]["c0"]
   a1 = NLLmin[nfit]["a1"]
   c2 = NLLmin[nfit]["c2"]
   #param = numpy.array([a0, c0/a0, a1/a0, c2/a0])
   #tool.set_parameters13(param)
   param = numpy.array([a0, c0, a1, c2])
   tool.set_parameters04(param)
   tool.set_sfact(1)
   print "best-fit parameters loaded from fit log file", fitlogfilename, ":"
   print "  NLL =", NLLmin[nfit]["value"], "(seen",\
                    NLLmin[nfit]["count"], "times)"
   print "   a0 =", a0, "+/-", a0err;
   print "   c0 =", c0, "+/-", c0err;
   print "   a1 =", a1, "+/-", a1err;
   print "   c2 =", c2, "+/-", c2err;
   if fitter:
      fitter.Command("set param 1 {0:20.12e}".format(param[0]))
      fitter.Command("set param 2 {0:20.12e}".format(param[1]))
      fitter.Command("set param 3 {0:20.12e}".format(param[2]))
      fitter.Command("set param 4 {0:20.12e}".format(param[3]))
      fitter.Command("show fcn")
   return NLLmin[nfit]

def plot_rho(i,j):
   """
   Plot the i,j component of the rho matrix (real part)
   """
   init_canvas()
   iname = str(i)
   jname = str(j)
   if i == -1:
      iname = "M"
   if j == -1:
      jname = "M"
   shift = [0, 0.004, -0.004, 0.008, -0.008, 0.012, -0.012]
   hist = {}
   for order in range(0,3):
      name = "rho{0}{1}{2}".format(iname, jname, order)
      title = "rho[{0},{1}]".format(i, j)
      hist[order] = TH1D(name, title, 30, shift[order], 1.2 + shift[order])
      hist[order].GetXaxis().SetTitle("|t| (GeV^2)")
      hist[order].GetYaxis().SetTitle(title)
      hist[order].GetYaxis().SetTitleOffset(1.6)
      hist[order].SetLineColor(order+1)
      hist[order].SetStats(0)
      for tbin in range(0,30):
         result = loadfit(tbin, -order-1)
         if result:
            errmat = TMatrixD(4,4)
            for row in range(0,4):
               for col in range(0,4):
                  errmat[row][col] = result["errmat"][row][col]
            rhocovar = TMatrixD(tool.get_PartialD())
            rhocovar *= errmat
            rhocovar *= TMatrixD(TMatrixD.kTransposed, tool.get_PartialD())
            traceRho = 2 * tool.get_Re_rho(1,1) + tool.get_Re_rho(0,0)
            hist[order].SetBinContent(tbin + 1, tool.get_Re_rho(i,j) / traceRho)
            rho11 = tool.get_Re_rho(1,1) / traceRho
            rho00 = tool.get_Re_rho(0,0) / traceRho
            rhoM1 = tool.get_Re_rho(-1,1) / traceRho
            rho10 = tool.get_Re_rho(1,0) / traceRho
            N = TMatrixD(4,4)
            N[0][0] = 2
            N[0][1] = 1
            N[0][2] = 0
            N[0][3] = 0
            N[1][0] = -2 * rho00 / traceRho
            N[1][1] = 2 * rho11 / traceRho
            N[1][2] = 0
            N[1][3] = 0
            N[2][0] = -2 * rhoM1 / traceRho
            N[2][1] = -rhoM1 / traceRho
            N[2][2] = 1 / traceRho
            N[2][3] = 0
            N[3][0] = -2 * rho10 / traceRho
            N[3][1] = -rho10 / traceRho
            N[3][2] = 0
            N[3][3] = 1 / traceRho
            NT = TMatrixD(4,4)
            NT.Transpose(N)
            Nrhocovar = TMatrixD(4,4)
            Nrhocovar.Mult(N, rhocovar)
            rhocovar.Mult(Nrhocovar, NT)
            if (i == 1 and j == 1) or (i == -1 and j == -1):
               error = 0.5 * math.sqrt(rhocovar[1][1])
            elif i == 0 and j == 0:
               error = math.sqrt(rhocovar[1][1])
            elif (i == 1 and j == -1) or (i == -1 and j == 1):
               error = math.sqrt(rhocovar[2][2])
            else:
               error = math.sqrt(rhocovar[3][3])
            hist[order].SetBinError(tbin + 1, error)
      if order == 0:
         c1.SetLogy(0)
         hist[order].Draw()
      elif hist[order].GetMaximum() != 0 or hist[order].GetMinimum() != 0:
         hist[order].Draw("same")
      c1.Update()
   return hist

def plot_dsdt():
   """
   Plot dsigma/dt
   """
   init_canvas()
   hist = {}
   shift = [0, 0.004, -0.004, 0.008, -0.008, 0.012, -0.012]
   for order in range(0,3):
      name = "dsdt{0}".format(order)
      title = "differential cross section dsigma/dt"
      hist[order] = TH1D(name, title, 30, shift[order], 1.2 + shift[order])
      hist[order].GetXaxis().SetTitle("|t| (GeV^2)")
      hist[order].GetYaxis().SetTitle("d#sigma/dt (ub/GeV^2)")
      hist[order].GetYaxis().SetTitleOffset(1.6)
      hist[order].SetLineColor(order+1)
      hist[order].SetStats(0)
      for tbin in range(0,30):
         result = loadfit(tbin, -order-1)
         if result:
            errmat = TMatrixD(4,4)
            for row in range(0,4):
               for col in range(0,4):
                  errmat[row][col] = result["errmat"][row][col]
            rhocovar = TMatrixD(tool.get_PartialD())
            rhocovar *= errmat
            rhocovar *= TMatrixD(TMatrixD.kTransposed, tool.get_PartialD())
            traceRho = 2 * tool.get_Re_rho(1,1) + tool.get_Re_rho(0,0)
            traceVar = 4 * rhocovar[0][0] + 4 * rhocovar[1][0] + rhocovar[1][1]
            Wintegral = 8 * math.pi / 3
            dt = tbins[tbin][1] - tbins[tbin][0]
            dsdt = traceRho * Wintegral / dt
            dsdt *= Nmcgen[tbin] / (4 * math.pi * integrated_luminosity)
            dsdt /= decay_branching_ratio
            dsdt /= protons_per_atom
            hist[order].SetBinContent(tbin + 1, dsdt)
            error = math.sqrt(traceVar) * Wintegral / dt
            error *= Nmcgen[tbin] / (4 * math.pi * integrated_luminosity)
            error /= decay_branching_ratio
            error /= protons_per_atom
            hist[order].SetBinError(tbin + 1, error)
      if order == 0:
         c1.SetLogy()
         hist[order].Draw()
      elif hist[order].GetMaximum() > 0:
         hist[order].Draw("same")
      c1.Update()
   return hist


############################
## Initialization Section ##
############################

# read the MC statistics from an input file

Ebins = {}
tbins = {}
Nmcgen = {}
for line in open("PWAmcgen.txt"):
   match = re.match(r"^ +([0-9]+)", line)
   if match:
      sample = int(match.group(1))
      item = line.strip().split()
      tbins[sample] = (float(item[1]), float(item[2]))
      Ebins[sample] = (float(item[3]), float(item[4]))
      Nmcgen[sample] = float(item[5])

integrated_luminosity = 1.89e7 # inverse microbarns
decay_branching_ratio = 0.0828
protons_per_atom = 4
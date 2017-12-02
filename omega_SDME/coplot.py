#!/bin/env python
#
# coplot.py - utility functions for overlaying histograms on a canvas
#
# author: richard.t.jones at uconn.edu
# version: november 22, 2017

from ROOT import *
import math

def h1(hnames):
   """
   Take a variable-length list of histograms and overlay the contents
   on a single plot, adjusting the vertical scale to accommodate the
   largest and smallest values. Argument hnames must be an array of
   strings, which are the names of the histograms to be plotted.
   """
   global htemp
   htemp = []
   hcolo = []
   ymin = 1e99
   ymax = -1e99
   for n in range(0, len(hnames)):
      if len(hnames[n]) > 0:
         h = gROOT.FindObject(hnames[n])
         if h:
            h.GetYaxis().SetTitleOffset(1.4)
            tempname = "htemp" + str(n)
            oo = gROOT.FindObject(tempname)
            if oo:
               oo.Delete()
            htemp.append(h.Clone(tempname))
            htemp[-1].SetDirectory(0)
            htemp[-1].SetStats(0)
            hcolo.append(n)
            for b in range(1, htemp[-1].GetNbinsX()+1):
               y = htemp[-1].GetBinContent(b)
               ey = htemp[-1].GetBinError(b)
               if math.isnan(y) or math.isnan(ey):
                  htemp[-1].SetBinContent(b,0)
                  htemp[-1].SetBinError(b,0)
                  continue
               if y - ey < ymin:
                  ymin = y - ey
               elif y + ey > ymax:
                  ymax = y + ey
            htemp[-1].Draw()
   global c1
   c1 = gROOT.FindObject("c1")
   if not c1:
      c1 = TCanvas("c1", "c1")
   if len(htemp) > 0:
      dy = ymax - ymin
      htemp[0].SetMaximum(ymax + dy/10)
      if ymin < 0 or dy < ymax * 0.2:
         htemp[0].SetMinimum(ymin - dy/10)
      else:
         htemp[0].SetMinimum(0)
   htemp[0].Draw()
   for n in range(1, len(htemp)):
      htemp[n].SetLineColor(hcolo[n] + 1)
      htemp[n].Draw("same")


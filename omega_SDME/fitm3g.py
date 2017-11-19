#!/bin/env python

from ROOT import *
import m3gfit
import math

f = TFile("ntuple_real.root")
s = f.Get("m3g0").Clone("s")
for n in range(1,44):
  h = f.Get("m3g" + str(n))
  s.Add(h)
s.ProjectionX("sx", 1, 1)

f1 = TF1("f1", m3gfit.gaus_on_poly, 0.4, 1.2, 6)
f1.SetParameter(0, 0.8)
f1.SetParameter(1, 0.08)
f1.SetParameter(2, 1e5)
f1.SetParameter(3, 1e4)
f1.SetParameter(4, 1e4)
f1.SetParameter(5, 1e4)
sx.Fit(f1, "", "", 0.4, 1.2)
mean = f1.GetParameter(0)
sigma = f1.GetParameter(1)
peakheight = f1.GetParameter(2)
bgheight = f1.GetParameter(3)
wbin = s.GetXaxis().GetBinWidth(1)
print "integral is", (2 * math.pi)**0.5 * sigma * peakheight / wbin
print "peak centroid is", mean, "GeV"
print "peak sigma is", sigma, "GeV"
print "sig/bg is", peakheight / bgheight
raw_input("what? ")

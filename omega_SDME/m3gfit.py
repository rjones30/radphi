#!/bin/env python

from ROOT import *
import math

def gaus_on_poly(m, par):
   """
   Declares a root-style parametric fitting function describing
   a double gaussian peak on top of a polynomial background.
    par[0] = gaussian mean 1
    par[1] = gaussian sigma 1
    par[2] = max height of gaussian 1
    par[3] = gaussian mean 2
    par[4] = gaussian sigma 2
    par[5] = max height of gaussian 2
    par[6] = order 0 polynonmial coefficient
    ...
    par[n] = order n-3 polynonmial coefficient
   """
   y = par[2] * math.exp(-0.5 * ((m[0] - par[0]) / par[1])**2)
   y += par[5] * math.exp(-0.5 * ((m[0] - par[3]) / par[4])**2)
   for o in range(0,99):
      p = o + 6
      if p < len(par):
         y += par[p] * (m[0] - par[0])**o
   return y

f1 = TF1("f1", gaus_on_poly, 0.4, 1.2, 8)
f1.SetParameter(0, 0.80)
f1.SetParameter(1, 0.06)
f1.SetParameter(2, 1e5)
f1.SetParameter(3, 0.60)
f1.SetParameter(4, 0.05)
f1.SetParameter(5, 1e3)
f1.SetParameter(6, 1e4)
f1.SetParameter(7, 1e4)
#f1.SetParameter(8, 1e4)

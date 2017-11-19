#!/usr/bin/python
#
# PWAmcgen.py - quick script for counting the number of MC generated
#               events in a given s,t bin.
#
# authors: richard.t.jones at uconn.edu, fridah.mokaya at uconn.edu
# version: november 4, 2016

import os
import glob
import sys
import numpy

from ROOT import *

h1 = TChain("h1", "mctuple output chain")
rootdir = "/pnfs4/phys.uconn.edu/data/Gluex/radphi/omega_sim-9-2014"
for n in range(1,102):
   h1.Add(rootdir + "/omega2p_gen-{0:d}-0.root".format(n))


tbins = [[0.0, 0.04], \
         [0.04, 0.08], \
         [0.08, 0.12], \
         [0.12, 0.16], \
         [0.16, 0.20], \
         [0.20, 0.24], \
         [0.24, 0.28], \
         [0.28, 0.32], \
         [0.32, 0.36], \
         [0.36, 0.4], \
         [0.40, 0.44], \
         [0.44, 0.48], \
         [0.48, 0.52], \
         [0.52, 0.56], \
         [0.56, 0.60], \
         [0.60, 0.64], \
         [0.64, 0.68], \
         [0.68, 0.72], \
         [0.72, 0.76], \
         [0.76, 0.80], \
         [0.80, 0.84], \
         [0.84, 0.88], \
         [0.88, 0.92], \
         [0.92, 0.96], \
         [0.96, 1.00], \
         [1.00, 1.04], \
         [1.04, 1.08], \
         [1.08, 1.12], \
         [1.12, 1.16], \
         [1.16, 1.20]]

Ebins = [[4.39, 5.39],
         [4.39, 4.60],
         [4.60, 4.90],
         [4.90, 5.15],
         [5.15, 5.39]]

proton_mass = 0.938272
Nmc = h1.GetEntries() / 100
print Nmc
tlimits = {}
Elimits = {}
Nmcgen = {}

# sample numbers 0..29 are fine-grained in t, coarse-grained in E
# sample numbers 100..114 are fine-grained in t, first fine E bin
# and so on for samples 200..214, 300..314, 400..414.

for Ebin in range(0, len(Ebins)):
   for sample in range(0, len(tbins)):
      if Ebin == 0:
         tlimits[sample] = tbins[sample]
      elif sample < len(tbins) / 2:
         tlimits[Ebin*100+sample] = [tbins[2*sample][0], tbins[2*sample+1][1]]
      else:
         continue
      Elimits[Ebin*100+sample] = Ebins[Ebin]
      Nmcgen[Ebin*100+sample] = 0

# count the MC generated statistics for all bins in one pass

for n in range(0, Nmc):
   h1.GetEntry(n)
   E = h1.momi[0]
   #tabs = h1.momf[13]**2 + h1.momf[14]**2 + h1.momf[15]**2 - \
   #       (h1.momf[12] - proton_mass)**2
   tabs = (h1.momf[1] + h1.momf[5] + h1.momf[9] - h1.momi[1])**2 + \
          (h1.momf[2] + h1.momf[6] + h1.momf[10] - h1.momi[2])**2 + \
          (h1.momf[3] + h1.momf[7] + h1.momf[11] - h1.momi[3])**2 - \
          (h1.momf[0] + h1.momf[4] + h1.momf[8] - h1.momi[0])**2
   for sample in tlimits.keys():
      if E >= Elimits[sample][0] and E <= Elimits[sample][1] and \
         tabs >= tlimits[sample][0] and tabs <= tlimits[sample][1]:
            Nmcgen[sample] += 1

# print out the results as a table

print "Sample     |t|low    |t|high      Elow      Ehigh     Nmcgen"
for sample in tlimits.keys():
   print "{0:5d} {1:12f} {2:12f} {3:12f} {4:12f} {5:15d}".format(sample,
          tlimits[sample][0], tlimits[sample][1],
          Elimits[sample][0], Elimits[sample][1], Nmcgen[sample])

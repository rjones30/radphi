#!/usr/bin/python
#
# frames.py - compares the decay theta values between the s-channel
#             helicity (SH) and the Gottfried-Jackson (GJ) frames
#             for the reaction gamma p -> p omega -> p pi0 gamma.
#
# author: richard.t.jones at uconn.edu
# version: july 22, 2015

from ROOT import *
from array import array
import random
import math

# all energies, momenta, masses are in GeV, c=1
# all angles are in radians

mProton = 0.938
mOmega = 0.780
mPi0 = 0.135

costSH = array('d', [0])
costGJ = array('d', [0])
phiSH = array('d', [0])

rootfile = TFile("frames.root", "recreate")

def sMandelstam(labE):
   """
   Computes Mandelstam variable s for the reaction gamma,p -> omega,p
   where labE is the photon energy in the initial proton rest frame.
   """
   return (labE + mProton) ** 2 - labE ** 2

def tlimits(s):
   """
   Computes the kinematic limits on Mandelstam t for the reaction
   gamma,p -> omega,p for a given value of Mandelstam s.
   """
   qin = (s - mProton ** 2) / (2 * sqrt(s))
   qout = sqrt(((s + mOmega ** 2 - mProton ** 2) ** 2) / (4 * s) - mOmega ** 2)
   Ein = qin
   Eout = sqrt(qout ** 2 + mOmega ** 2)
   tmin = (Ein - Eout) ** 2 - (qin + qout) ** 2
   tmax = (Ein - Eout) ** 2 - (qin - qout) ** 2
   return tmin, tmax

def costhetaGJ(s, t, costhetaSH, phiSH):
   """
   Computes the cos(theta) of the decay pi0 in the Gottfried-Jackson
   frame, for the reaction gamma,p -> omega,p and omega -> pi0,gamma
   for the given s and t, and helicity frame decay angles of the pi0
   given by cos(thetaSH) and phiSH [radians].
   """
   qOmegaDecay = (mOmega ** 2 - mPi0 ** 2) / (2 * mOmega)
   sinthetaSH = sqrt(1 - costhetaSH ** 2)
   pPi0SH = [qOmegaDecay * sinthetaSH * cos(phiSH), \
             qOmegaDecay * sinthetaSH * sin(phiSH), \
             qOmegaDecay * costhetaSH, \
             sqrt(qOmegaDecay ** 2 + mPi0 ** 2)]
   qout = sqrt(((s + mOmega ** 2 - mProton ** 2) ** 2) / (4 * s) - mOmega ** 2)
   Eout = sqrt(qout ** 2 + mOmega ** 2)
   boostGamma = Eout / mOmega
   boostBeta = qout / Eout
   qin = (s - mProton ** 2) / (2 * sqrt(s))
   Ein = qin
   pBeamCM = [0, 0, qin, Ein]
   costhetaCM = (t - (Ein - Eout) ** 2 + qin ** 2 + qout ** 2) \
                / (2 * qin * qout)
   sinthetaCM = sqrt(1 - costhetaCM ** 2)
   pBeamCM = [pBeamCM[0] * costhetaCM - pBeamCM[2] * sinthetaCM, \
              pBeamCM[1], \
              pBeamCM[2] * costhetaCM + pBeamCM[0] * sinthetaCM, \
              pBeamCM[3]]
   pBeamSH = [pBeamCM[0], \
              pBeamCM[1], \
              boostGamma * (pBeamCM[2] - boostBeta * pBeamCM[3]), \
              boostGamma * (pBeamCM[3] - boostBeta * pBeamCM[2])]
   return (pBeamSH[0] * pPi0SH[0] + pBeamSH[1] * pPi0SH[1] + \
           pBeamSH[2] * pPi0SH[2]) / (qOmegaDecay * pBeamSH[3])

def costGJ_vs_costSH(s, t, nevents):
   """
   Creates a Monte Carlo sample of random omega decays omega -> pi0,gamma
   from the reaction gamma,p -> omega,p at fixed s and t. A count of nevents
   are generated uniformly in costhetaSH,phiSH and the results are returned
   as a tree containing costhetaSH,costhetaGJ,phiSH for each event.
   """
   tree = TTree("GJvsSH", "cos(thetaGJ) vs cos(thetaSH)")
   tree.Branch("costSH", costSH, "costSH/D")
   tree.Branch("costGJ", costGJ, "costGJ/D")
   tree.Branch("phiSH", phiSH, "phiSH/D")
   global mOmega
   for event in range(0, nevents):
      mOmega = random.uniform(0.6, 0.9)
      costSH[0] = random.uniform(-0.8, 0.8)
      phiSH[0] = random.uniform(3.1, 3.2)
      costGJ[0] = costhetaGJ(s, t, costSH[0], phiSH[0])
      tree.Fill()
   tree.Write()
   return tree

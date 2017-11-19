#!/bin/env python
#
# svdplayer.py - play with singular value decomposition of random
#                matrices in the form of a vector meson SDM.
#
# author: richard.t.jones at uconn.edu
# version: november 17, 2017

from ROOT import *

rnd = TRandom3(0)

while True:
   rho = TMatrixD(3,3)
   for n in range(0,3):
      aRe = 20*(rnd.Rndm() - 0.5)
      aIm = 20*(rnd.Rndm() - 0.5)
      bRe = 20*(rnd.Rndm() - 0.5)
      bIm = 20*(rnd.Rndm() - 0.5)
      cRe = 20*(rnd.Rndm() - 0.5)
      cIm = 20*(rnd.Rndm() - 0.5)
      rho[0][0] += aRe*aRe + aIm*aIm
      rho[0][1] += aRe*bRe + aIm*bIm
      rho[0][2] += aRe*cRe + aIm*cIm
      rho[1][1] += bRe*bRe + bIm*bIm
      rho[1][2] += bRe*cRe + bIm*cIm
      rho[2][2] += cRe*cRe + cIm*cIm
   rho[0][0] += rho[2][2]
   rho[0][1] -= rho[1][2]
   rho[0][2] += rho[0][2]
   rho[1][1] += rho[1][1]
   rho[1][2] = -rho[0][1]
   rho[2][2] = +rho[0][0]
   rho[1][0] = rho[0][1]
   rho[2][0] = rho[0][2]
   rho[2][1] = rho[1][2]
   print "rho: ", rho.Print()
   svd = TDecompSVD(3,3)
   svd.SetMatrix(rho)
   svd.Decompose()
   sig = svd.GetSig()
   print "sig: ", sig.Print()
   U = svd.GetU()
   print "U:", U.Print()
   V = svd.GetV()
   print "V:", V.Print()
   raw_input("ok? ")

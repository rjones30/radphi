#!/bin/env python
#
# logs2out.py - script to create Fridah-style PWAfitNNN.out logs
#               from pwa_omega.sub job stdout logs
#
# author: richard.t.jones at uconn.edu
# version: november 15, 2016

import re
import sys

for proc in range(0,90):
   if proc < 30:
      sample = proc
   elif proc < 45:
      sample = proc - 30 + 100
   elif proc < 60:
      sample = proc - 45 + 200
   elif proc < 75:
      sample = proc - 60 + 300
   elif proc < 90:
      sample = proc - 75 + 400
   out = 0
   for line in open("logs.d/stdout.{0}".format(proc)):
      if out:
         out.write(line)
         if re.match(r">>>>> Completed fit", line):
            out = 0
            continue
      match = re.match(r">>>>> Starting fit ([0-9]+) <<<<<", line)
      if match:
         fit = int(match.group(1))
         filename = "fitresults/pwafits{0}/PWAfits{1}.out".format(sample, fit)
         out = open(filename, "w")
         print "writing", filename

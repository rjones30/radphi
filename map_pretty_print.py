#!/bin/env python
#
# map_pretty_print.py - a simple script to convert the ascii representation
#                       of a mapmanager map (from map2ascii) into a human
#                       readable format.
#
# authors: richard.t.jones at uconn.edu and fridah.mokaya at uconn.edu
# version: march 1, 2017

import sys

if len(sys.argv) < 2:
   print "Usage: map_pretty_print.py <mapfile.ascii>"
   sys.exit(1)

try:
   mapfile = open(sys.argv[1])
except:
   print "Unable to open input map file", sys.argv[1]
   sys.exit(2)

for line in mapfile:
   fields = line.rstrip().split()
   if len(fields) == 0:
      continue
   elif line[0] == "S":
      system = ""
      for i in range(1, len(fields)):
         if int(fields[i], 16) == 0:
            break
         system += "{0:c}".format(int(fields[i], 16))
      print "System", system
   elif line[0] == "I":
      tablelen = int(fields[1])
      tabletype = int(fields[2])
      table = ""
      for i in range(3, len(fields)):
         if int(fields[i], 16) == 0:
            break
         table += "{0:c}".format(int(fields[i], 16))
      print "  table", table
   elif line[0] == "F":
      run = int(fields[1])
      print "    new run range starting at", run
   elif line[0] == "A":
      run = int(fields[1])
      tabletype = 91
      print "    new run range starting at", run
   elif line[0] == "C":
      run = int(fields[1])
      arraysize = int(fields[2])
      print "    new run range starting at", run
      print "      ",
   elif line[0] == "E":
      continue
   elif line[0] == " " and tabletype == 1:
      print "      ",
      for i in range(0, min(len(fields), 5)):
         print "{0:12.8f}".format(float(fields[i])),
      print ""
      print "      ",
      for i in range(5, len(fields)):
         print "{0:12.8f}".format(float(fields[i])),
      print ""
   elif line[0] == " " and tabletype == 2:
      for i in range(0, len(fields)):
         if int(fields[i], 16) == 0:
            print ""
            break
         elif int(fields[i], 16) == 10:
            print "\n      ",
         else:
            sys.stdout.write("{0:c}".format(int(fields[i], 16)))
   elif line[0] == " " and tabletype == 91:
      print "      ",
      for i in range(0, min(len(fields), 5)):
         print "{0:8d}".format(int(fields[i])),
      print ""
      print "      ",
      for i in range(5, len(fields)):
         print "{0:8d}".format(int(fields[i])),
      print ""
   else:
      print "Unknown descriptor letter", fields[0], ", quitting!"
      sys.exit(3)

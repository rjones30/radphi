//
// PWAmcgen.py - quick utility for counting the number of MC generated
//               events in a given s,t bin, ported from python to C++
//               for reasons of computational speed.
//
// authors: richard.t.jones at uconn.edu, fridah.mokaya at uconn.edu
// version: november 4, 2016

#include <iostream>
#include <sstream>
#include <string>
#include <map>

#include <TChain.h>

double tbins[30][2] =
        {{0.0, 0.04}, 
         {0.04, 0.08}, 
         {0.08, 0.12}, 
         {0.12, 0.16}, 
         {0.16, 0.20}, 
         {0.20, 0.24}, 
         {0.24, 0.28}, 
         {0.28, 0.32}, 
         {0.32, 0.36}, 
         {0.36, 0.4}, 
         {0.40, 0.44}, 
         {0.44, 0.48}, 
         {0.48, 0.52}, 
         {0.52, 0.56}, 
         {0.56, 0.60}, 
         {0.60, 0.64}, 
         {0.64, 0.68}, 
         {0.68, 0.72}, 
         {0.72, 0.76}, 
         {0.76, 0.80}, 
         {0.80, 0.84}, 
         {0.84, 0.88}, 
         {0.88, 0.92}, 
         {0.92, 0.96}, 
         {0.96, 1.00}, 
         {1.00, 1.04}, 
         {1.04, 1.08}, 
         {1.08, 1.12}, 
         {1.12, 1.16}, 
         {1.16, 1.20}};

double Ebins[5][2] =
        {{4.39, 5.39},
         {4.39, 4.60},
         {4.60, 4.90},
         {4.90, 5.15},
         {5.15, 5.39}};

int main()
{
   TChain h1("h1", "mctuple output chain");
   std::string rootdir = "/pnfs4/phys.uconn.edu/data/Gluex/radphi/omega_sim-9-2014";
   for (int n=1; n < 102; ++n) {
      std::stringstream rootfile;
      rootfile << rootdir << "/omega2p_gen-" << n << "-0.root";
      h1.Add(rootfile.str().c_str());
   }
   int Nmc = h1.GetEntries();
   std::cout << "Nmc=" << Nmc << std::endl;
   std::map<int, double[2]> tlimits;
   std::map<int, double[2]> Elimits;
   std::map<int, double> Nmcgen;

   // sample numbers 0..29 are fine-grained in t, coarse-grained in E
   // sample numbers 100..114 are fine-grained in t, first fine E bin
   // and so on for samples 200..214, 300..314, 400..414.

   for (int Ebin=0; Ebin < 5; ++Ebin) {
      for (int sample=0; sample < 30; ++sample) {
         if (Ebin == 0) {
            tlimits[sample][0] = tbins[sample][0];
            tlimits[sample][1] = tbins[sample][1];
         }
         else if (sample < 15) {
            tlimits[Ebin*100+sample][0] = tbins[2*sample][0];
            tlimits[Ebin*100+sample][1] = tbins[2*sample+1][1];
         }
         else {
            continue;
         }
         Elimits[Ebin*100+sample][0] = Ebins[Ebin][0];
         Elimits[Ebin*100+sample][1] = Ebins[Ebin][1];
         Nmcgen[Ebin*100+sample] = 0;
      }
   }

   // count the MC generated statistics for all bins in one pass

   float momi[1][4];
   float momf[4][4];
   h1.SetBranchAddress("momi", momi);
   h1.SetBranchAddress("momf", momf);
   for (int n=0; n < Nmc; ++n) {
      h1.GetEntry(n);
      double E = momi[0][0];
      double tabs = pow(momf[0][1] + momf[1][1] + momf[2][1] - momi[0][1], 2) + 
                    pow(momf[0][2] + momf[1][2] + momf[2][2] - momi[0][2], 2) + 
                    pow(momf[0][3] + momf[1][3] + momf[2][3] - momi[0][3], 2) - 
                    pow(momf[0][0] + momf[1][0] + momf[2][0] - momi[0][0], 2);
      std::map<int, double[2]>::iterator iter;
      for (iter = tlimits.begin(); iter != tlimits.end(); ++iter) {
         if (tabs >= iter->second[0] && tabs <= iter->second[1] &&
             E >= Elimits[iter->first][0] && E <= Elimits[iter->first][1])
         {
            ++Nmcgen[iter->first];
         }
      }
   }

   // print out the results as a table

   std::cout << "Sample     |t|low    |t|high      Elow      Ehigh     Nmcgen" << std::endl;
   std::map<int, double[2]>::iterator iter;
   for (iter = tlimits.begin(); iter != tlimits.end(); ++iter) {
      std::cout << iter->first << " " << iter->second[0] << " " << iter->second[1] << " "
                << Elimits[iter->first][0] << " " << Elimits[iter->first][1] << " "
                << Nmcgen[iter->first] << std::endl;
   }

   return 0;
}

#include <TH1D.h>
#include <math.h>

double r1 = 7.;
double r2 = 35.;
double dr = 0.5;
double phi_dr = -30 * 3.1416/180;

TH1D *dphi()
{
   TH1D *h = new TH1D("dphi","delta phi / phi",360,-180,180);
   double lastphi1=0, lastphi2=0;
   for (int i=1; i<=360; ++i) {
      double phi0 = h->GetBinCenter(i) * (M_PI/180);
      double x1 = r1*cos(phi0);
      double y1 = r1*sin(phi0);
      double x2 = -r2*cos(phi0);
      double y2 = -r2*sin(phi0);
      x1 -= dr*cos(phi_dr);
      x2 -= dr*cos(phi_dr);
      y1 -= dr*sin(phi_dr);
      y2 -= dr*sin(phi_dr);
      double phi1 = atan2(y1,x1);
      double phi2 = atan2(y2,x2);
      double deltaphi = (phi1 - lastphi1) * (180/M_PI);
      lastphi1 = phi1;
      lastphi2 = phi2;
      deltaphi += (deltaphi < 0)? 360 : 0;
      if (i > 1)
         h->SetBinContent(i,4750*deltaphi);
   }
   return h;
}

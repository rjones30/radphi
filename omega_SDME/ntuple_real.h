//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov  9 16:21:02 2017 by ROOT version 6.06/08
// from TTree h1/Radphi
// found on file: /pnfs4/phys.uconn.edu/data/Gluex/radphi/pass-6-2014/r7415-0.root
//////////////////////////////////////////////////////////

#ifndef ntuple_real_h
#define ntuple_real_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>

#include <vector>
#include <string>

// Headers needed by this particular selector


class ntuple_real : public TSelector {
public :
   //TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Output trees and histograms
   TH2D *h2mass3g[100];
   TH2D *h2dphi[100];
   TH2D *h2E3g[100];
   TH2D *h2heli[100];

   TFile *dumpfile;
   TTree *tkine;
   struct kinetree {
      double Ebeam;
      double plabfwd[4]; // stored in TLorentzVector order px,py,pz,E 
      double plabpi0[4];
      double plabbac[4];
      double plabrec[4];
      double plabini[4];
      double plabtot[4];
      double preafwd[4];
      double preapi0[4];
      double preabac[4];
      double prearec[4];
      double preaini[4];
      double preatot[4];
      double preabea[4];
      double phelfwd[4];
      double phelpi0[4];
      double phelbac[4];
      double phelrec[4];
      double phelini[4];
      double pheltot[4];
      double phelbea[4];
      double theta0;
      double theta1;
      double costhetahel;
      double phihel;
   } tkinerow;

   TFile *skimfile;
   TTree *tskim;
   struct skimtree {
      double runevent;
      double Ebeam;
      double Emiss;
      double m3g;
      double dphi;
      double tabs;
      double costhetahel;
      double phihel;
      double weight;
   } tskimrow;

   // Utility member functions
   TLorentzVector GetPpi0Helicity(double Ebeam);

   // Declaration of leaf types
   Int_t           runno;
   Int_t           eventno;
   Int_t           ismc;
   Int_t           evtype;
   Int_t           memam;
   Int_t           cmam;
   Int_t           emam;
   Int_t           m2mam;
   Int_t           nbsd;
   Int_t           ibsd[400];   //[nbsd]
   Int_t           absd[400];   //[nbsd]
   Int_t           tbsd[400];   //[nbsd]
   Int_t           nbgv;
   Int_t           ibgv[400];   //[nbgv]
   Int_t           abgv[400];   //[nbgv]
   Int_t           tbgv[400];   //[nbgv]
   Int_t           ncpv;
   Int_t           icpv[250];   //[ncpv]
   Int_t           acpv[250];   //[ncpv]
   Int_t           tcpv[250];   //[ncpv]
   Int_t           nupv;
   Int_t           iupv[100];   //[nupv]
   Int_t           aupv[100];   //[nupv]
   Int_t           tupv[100];   //[nupv]
   Int_t           nlgd;
   Int_t           ilgd[640];   //[nlgd]
   Int_t           jlgd[640];   //[nlgd]
   Int_t           algd[640];   //[nlgd]
   Int_t           ntag;
   Int_t           itag[320];   //[ntag]
   Int_t           ttag[320];   //[ntag]
   Int_t           nhbsd;
   Int_t           chbsd[400];   //[nhbsd]
   Float_t         Ebsd[400];   //[nhbsd]
   Int_t           ntbsd[400];   //[nhbsd]
   Int_t           t1bsd[400];   //[nhbsd]
   Int_t           npix;
   Int_t           ipix[200];   //[npix]
   Int_t           rpix[200];   //[npix]
   Int_t           lpix[200];   //[npix]
   Int_t           spix[200];   //[npix]
   Float_t         zpix[200];   //[npix]
   Float_t         phipix[200];   //[npix]
   Float_t         rpixt[200];   //[npix]
   Float_t         lpixt[200];   //[npix]
   Float_t         spixt[200];   //[npix]
   Float_t         rpixe[200];   //[npix]
   Float_t         lpixe[200];   //[npix]
   Float_t         spixe[200];   //[npix]
   Int_t           nhbgv;
   Int_t           chbgv[400];   //[nhbgv]
   Float_t         Ebgvdwn[400];   //[nhbgv]
   Float_t         Ebgvup[400];   //[nhbgv]
   Int_t           ntbgvdwn[400];   //[nhbgv]
   Int_t           t1bgvdwn[400];   //[nhbgv]
   Int_t           ntbgvup[400];   //[nhbgv]
   Int_t           t1bgvup[400];   //[nhbgv]
   Int_t           nhcpv;
   Int_t           chcpv[250];   //[nhcpv]
   Float_t         Ecpv[250];   //[nhcpv]
   Int_t           ntcpv[250];   //[nhcpv]
   Int_t           t1cpv[250];   //[nhcpv]
   Int_t           nhupv;
   Int_t           chupv[50];   //[nhupv]
   Float_t         Eupv[50];   //[nhupv]
   Int_t           ntupv[50];   //[nhupv]
   Int_t           t1upv[50];   //[nhupv]
   Int_t           nhtag;
   Int_t           chtag[320];   //[nhtag]
   Float_t         Etag0[320];   //[nhtag]
   Float_t         Etag1[320];   //[nhtag]
   Int_t           nttagl[320];   //[nhtag]
   Int_t           nttagr[320];   //[nhtag]
   Int_t           t1tagl[320];   //[nhtag]
   Int_t           t1tagr[320];   //[nhtag]
   Int_t           ntimes;
   Float_t         le[500];   //[ntimes]
   Int_t           nrec;
   Float_t         trec0;
   Float_t         therec[100];   //[nrec]
   Float_t         phirec[100];   //[nrec]
   Float_t         derec[100];   //[nrec]
   Float_t         Erec[100];   //[nrec]
   Float_t         trec[100];   //[nrec]
   Int_t           mrec[100];   //[nrec]
   Int_t           ncoin;
   Int_t           cochan[30];   //[ncoin]
   Float_t         cotime[30];   //[ncoin]
   Float_t         coenergy[30];   //[ncoin]
   Float_t         tagweight[30];   //[ncoin]
   Int_t           nhlgd;
   Int_t           chlgd[640];   //[nhlgd]
   Float_t         Elgd[640];   //[nhlgd]
   Int_t           clust[640];   //[nhlgd]
   Int_t           nphot;
   Int_t           nfrwd;
   Float_t         Efrwd;
   Float_t         pvect[400][4];   //[nphot]
   Int_t           nbcl;
   Float_t         bce[24];   //[nbcl]
   Float_t         bcphi[24];   //[nbcl]
   Float_t         bcz[24];   //[nbcl]
   Float_t         bct[24];   //[nbcl]
   Float_t         bcse[24];   //[nbcl]
   Int_t           nmes;
   Int_t           mtype[500];   //[nmes]
   Float_t         ptot[500][4];   //[nmes]
   Float_t         amass[500];   //[nmes]
   Int_t           idtype[500][2];   //[nmes]
   Int_t           ichild[500][2];   //[nmes]

   // List of branches
   TBranch        *b_runno;   //!
   TBranch        *b_eventno;   //!
   TBranch        *b_ismc;   //!
   TBranch        *b_evtype;   //!
   TBranch        *b_memam;   //!
   TBranch        *b_cmam;   //!
   TBranch        *b_emam;   //!
   TBranch        *b_m2mam;   //!
   TBranch        *b_nbsd;   //!
   TBranch        *b_ibsd;   //!
   TBranch        *b_absd;   //!
   TBranch        *b_tbsd;   //!
   TBranch        *b_nbgv;   //!
   TBranch        *b_ibgv;   //!
   TBranch        *b_abgv;   //!
   TBranch        *b_tbgv;   //!
   TBranch        *b_ncpv;   //!
   TBranch        *b_icpv;   //!
   TBranch        *b_acpv;   //!
   TBranch        *b_tcpv;   //!
   TBranch        *b_nupv;   //!
   TBranch        *b_iupv;   //!
   TBranch        *b_aupv;   //!
   TBranch        *b_tupv;   //!
   TBranch        *b_nlgd;   //!
   TBranch        *b_ilgd;   //!
   TBranch        *b_jlgd;   //!
   TBranch        *b_algd;   //!
   TBranch        *b_ntag;   //!
   TBranch        *b_itag;   //!
   TBranch        *b_ttag;   //!
   TBranch        *b_nhbsd;   //!
   TBranch        *b_chbsd;   //!
   TBranch        *b_Ebsd;   //!
   TBranch        *b_ntbsd;   //!
   TBranch        *b_t1bsd;   //!
   TBranch        *b_npix;   //!
   TBranch        *b_ipix;   //!
   TBranch        *b_rpix;   //!
   TBranch        *b_lpix;   //!
   TBranch        *b_spix;   //!
   TBranch        *b_zpix;   //!
   TBranch        *b_phipix;   //!
   TBranch        *b_rpixt;   //!
   TBranch        *b_lpixt;   //!
   TBranch        *b_spixt;   //!
   TBranch        *b_rpixe;   //!
   TBranch        *b_lpixe;   //!
   TBranch        *b_spixe;   //!
   TBranch        *b_nhbgv;   //!
   TBranch        *b_chbgv;   //!
   TBranch        *b_Ebgvdwn;   //!
   TBranch        *b_Ebgvup;   //!
   TBranch        *b_ntbgvdwn;   //!
   TBranch        *b_t1bgvdwn;   //!
   TBranch        *b_ntbgvup;   //!
   TBranch        *b_t1bgvup;   //!
   TBranch        *b_nhcpv;   //!
   TBranch        *b_chcpv;   //!
   TBranch        *b_Ecpv;   //!
   TBranch        *b_ntcpv;   //!
   TBranch        *b_t1cpv;   //!
   TBranch        *b_nhupv;   //!
   TBranch        *b_chupv;   //!
   TBranch        *b_Eupv;   //!
   TBranch        *b_ntupv;   //!
   TBranch        *b_t1upv;   //!
   TBranch        *b_nhtag;   //!
   TBranch        *b_chtag;   //!
   TBranch        *b_Etag0;   //!
   TBranch        *b_Etag1;   //!
   TBranch        *b_nttagl;   //!
   TBranch        *b_nttagr;   //!
   TBranch        *b_t1tagl;   //!
   TBranch        *b_t1tagr;   //!
   TBranch        *b_ntimes;   //!
   TBranch        *b_le;   //!
   TBranch        *b_nrec;   //!
   TBranch        *b_trec0;   //!
   TBranch        *b_therec;   //!
   TBranch        *b_phirec;   //!
   TBranch        *b_derec;   //!
   TBranch        *b_Erec;   //!
   TBranch        *b_trec;   //!
   TBranch        *b_mrec;   //!
   TBranch        *b_ncoin;   //!
   TBranch        *b_cochan;   //!
   TBranch        *b_cotime;   //!
   TBranch        *b_coenergy;   //!
   TBranch        *b_tagweight;   //!
   TBranch        *b_nhlgd;   //!
   TBranch        *b_chlgd;   //!
   TBranch        *b_Elgd;   //!
   TBranch        *b_clust;   //!
   TBranch        *b_nphot;   //!
   TBranch        *b_nfrwd;   //!
   TBranch        *b_Efrwd;   //!
   TBranch        *b_pvect;   //!
   TBranch        *b_nbcl;   //!
   TBranch        *b_bce;   //!
   TBranch        *b_bcphi;   //!
   TBranch        *b_bcz;   //!
   TBranch        *b_bct;   //!
   TBranch        *b_bcse;   //!
   TBranch        *b_nmes;   //!
   TBranch        *b_mtype;   //!
   TBranch        *b_ptot;   //!
   TBranch        *b_amass;   //!
   TBranch        *b_idtype;   //!
   TBranch        *b_ichild;   //!


   ntuple_real(TTree *tree =0);
   virtual ~ntuple_real();
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(ntuple_real,0);

};

#endif

#ifdef ntuple_real_cxx
ntuple_real::ntuple_real(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

   TString rootfile_default("/pnfs4/phys.uconn.edu/data/Gluex/radphi/pass-6-2014/r7415-0.root");
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(rootfile_default);
      if (!f || !f->IsOpen()) {
         f = new TFile(rootfile_default);
      }
      f->GetObject("h1",tree);

   }
   Init(tree);
}

ntuple_real::~ntuple_real()
{
   //if (!fChain) return;
   //delete fChain->GetCurrentFile();
}

void ntuple_real::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   //fReader.SetTree(tree);

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runno", &runno, &b_runno);
   fChain->SetBranchAddress("eventno", &eventno, &b_eventno);
   //fChain->SetBranchAddress("ismc", &ismc, &b_ismc);
   //fChain->SetBranchAddress("evtype", &evtype, &b_evtype);
   //fChain->SetBranchAddress("memam", &memam, &b_memam);
   //fChain->SetBranchAddress("cmam", &cmam, &b_cmam);
   //fChain->SetBranchAddress("emam", &emam, &b_emam);
   //fChain->SetBranchAddress("m2mam", &m2mam, &b_m2mam);
   //fChain->SetBranchAddress("nbsd", &nbsd, &b_nbsd);
   //fChain->SetBranchAddress("ibsd", ibsd, &b_ibsd);
   //fChain->SetBranchAddress("absd", absd, &b_absd);
   //fChain->SetBranchAddress("tbsd", tbsd, &b_tbsd);
   //fChain->SetBranchAddress("nbgv", &nbgv, &b_nbgv);
   //fChain->SetBranchAddress("ibgv", ibgv, &b_ibgv);
   //fChain->SetBranchAddress("abgv", abgv, &b_abgv);
   //fChain->SetBranchAddress("tbgv", tbgv, &b_tbgv);
   //fChain->SetBranchAddress("ncpv", &ncpv, &b_ncpv);
   //fChain->SetBranchAddress("icpv", icpv, &b_icpv);
   //fChain->SetBranchAddress("acpv", acpv, &b_acpv);
   //fChain->SetBranchAddress("tcpv", tcpv, &b_tcpv);
   //fChain->SetBranchAddress("nupv", &nupv, &b_nupv);
   //fChain->SetBranchAddress("iupv", iupv, &b_iupv);
   //fChain->SetBranchAddress("aupv", aupv, &b_aupv);
   //fChain->SetBranchAddress("tupv", tupv, &b_tupv);
   //fChain->SetBranchAddress("nlgd", &nlgd, &b_nlgd);
   //fChain->SetBranchAddress("ilgd", ilgd, &b_ilgd);
   //fChain->SetBranchAddress("jlgd", jlgd, &b_jlgd);
   //fChain->SetBranchAddress("algd", algd, &b_algd);
   //fChain->SetBranchAddress("ntag", &ntag, &b_ntag);
   //fChain->SetBranchAddress("itag", itag, &b_itag);
   //fChain->SetBranchAddress("ttag", ttag, &b_ttag);
   //fChain->SetBranchAddress("nhbsd", &nhbsd, &b_nhbsd);
   //fChain->SetBranchAddress("chbsd", chbsd, &b_chbsd);
   //fChain->SetBranchAddress("Ebsd", Ebsd, &b_Ebsd);
   //fChain->SetBranchAddress("ntbsd", ntbsd, &b_ntbsd);
   //fChain->SetBranchAddress("t1bsd", t1bsd, &b_t1bsd);
   //fChain->SetBranchAddress("npix", &npix, &b_npix);
   //fChain->SetBranchAddress("ipix", ipix, &b_ipix);
   //fChain->SetBranchAddress("rpix", rpix, &b_rpix);
   //fChain->SetBranchAddress("lpix", lpix, &b_lpix);
   //fChain->SetBranchAddress("spix", spix, &b_spix);
   //fChain->SetBranchAddress("zpix", zpix, &b_zpix);
   //fChain->SetBranchAddress("phipix", phipix, &b_phipix);
   //fChain->SetBranchAddress("rpixt", rpixt, &b_rpixt);
   //fChain->SetBranchAddress("lpixt", lpixt, &b_lpixt);
   //fChain->SetBranchAddress("spixt", spixt, &b_spixt);
   //fChain->SetBranchAddress("rpixe", rpixe, &b_rpixe);
   //fChain->SetBranchAddress("lpixe", lpixe, &b_lpixe);
   //fChain->SetBranchAddress("spixe", spixe, &b_spixe);
   //fChain->SetBranchAddress("nhbgv", &nhbgv, &b_nhbgv);
   //fChain->SetBranchAddress("chbgv", chbgv, &b_chbgv);
   //fChain->SetBranchAddress("Ebgvdwn", Ebgvdwn, &b_Ebgvdwn);
   //fChain->SetBranchAddress("Ebgvup", Ebgvup, &b_Ebgvup);
   //fChain->SetBranchAddress("ntbgvdwn", ntbgvdwn, &b_ntbgvdwn);
   //fChain->SetBranchAddress("t1bgvdwn", t1bgvdwn, &b_t1bgvdwn);
   //fChain->SetBranchAddress("ntbgvup", ntbgvup, &b_ntbgvup);
   //fChain->SetBranchAddress("t1bgvup", t1bgvup, &b_t1bgvup);
   fChain->SetBranchAddress("nhcpv", &nhcpv, &b_nhcpv);
   //fChain->SetBranchAddress("chcpv", chcpv, &b_chcpv);
   //fChain->SetBranchAddress("Ecpv", Ecpv, &b_Ecpv);
   fChain->SetBranchAddress("ntcpv", ntcpv, &b_ntcpv);
   fChain->SetBranchAddress("t1cpv", t1cpv, &b_t1cpv);
   //fChain->SetBranchAddress("nhupv", &nhupv, &b_nhupv);
   //fChain->SetBranchAddress("chupv", chupv, &b_chupv);
   //fChain->SetBranchAddress("Eupv", Eupv, &b_Eupv);
   //fChain->SetBranchAddress("ntupv", ntupv, &b_ntupv);
   //fChain->SetBranchAddress("t1upv", t1upv, &b_t1upv);
   //fChain->SetBranchAddress("nhtag", &nhtag, &b_nhtag);
   //fChain->SetBranchAddress("chtag", chtag, &b_chtag);
   //fChain->SetBranchAddress("Etag0", Etag0, &b_Etag0);
   //fChain->SetBranchAddress("Etag1", Etag1, &b_Etag1);
   //fChain->SetBranchAddress("nttagl", nttagl, &b_nttagl);
   //fChain->SetBranchAddress("nttagr", nttagr, &b_nttagr);
   //fChain->SetBranchAddress("t1tagl", t1tagl, &b_t1tagl);
   //fChain->SetBranchAddress("t1tagr", t1tagr, &b_t1tagr);
   //fChain->SetBranchAddress("ntimes", &ntimes, &b_ntimes);
   fChain->SetBranchAddress("le", le, &b_le);
   fChain->SetBranchAddress("nrec", &nrec, &b_nrec);
   fChain->SetBranchAddress("trec0", &trec0, &b_trec0);
   fChain->SetBranchAddress("therec", therec, &b_therec);
   fChain->SetBranchAddress("phirec", phirec, &b_phirec);
   //fChain->SetBranchAddress("derec", derec, &b_derec);
   //fChain->SetBranchAddress("Erec", Erec, &b_Erec);
   //fChain->SetBranchAddress("trec", trec, &b_trec);
   //fChain->SetBranchAddress("mrec", mrec, &b_mrec);
   fChain->SetBranchAddress("ncoin", &ncoin, &b_ncoin);
   //fChain->SetBranchAddress("cochan", cochan, &b_cochan);
   fChain->SetBranchAddress("cotime", cotime, &b_cotime);
   fChain->SetBranchAddress("coenergy", coenergy, &b_coenergy);
   fChain->SetBranchAddress("tagweight", tagweight, &b_tagweight);
   //fChain->SetBranchAddress("nhlgd", &nhlgd, &b_nhlgd);
   //fChain->SetBranchAddress("chlgd", chlgd, &b_chlgd);
   //fChain->SetBranchAddress("Elgd", Elgd, &b_Elgd);
   //fChain->SetBranchAddress("clust", clust, &b_clust);
   fChain->SetBranchAddress("nphot", &nphot, &b_nphot);
   fChain->SetBranchAddress("nfrwd", &nfrwd, &b_nfrwd);
   fChain->SetBranchAddress("Efrwd", &Efrwd, &b_Efrwd);
   fChain->SetBranchAddress("pvect", pvect, &b_pvect);
   //fChain->SetBranchAddress("nbcl", &nbcl, &b_nbcl);
   //fChain->SetBranchAddress("bce", bce, &b_bce);
   //fChain->SetBranchAddress("bcphi", bcphi, &b_bcphi);
   //fChain->SetBranchAddress("bcz", bcz, &b_bcz);
   //fChain->SetBranchAddress("bct", bct, &b_bct);
   //fChain->SetBranchAddress("bcse", bcse, &b_bcse);
   //fChain->SetBranchAddress("nmes", &nmes, &b_nmes);
   //fChain->SetBranchAddress("mtype", mtype, &b_mtype);
   //fChain->SetBranchAddress("ptot", ptot, &b_ptot);
   //fChain->SetBranchAddress("amass", amass, &b_amass);
   //fChain->SetBranchAddress("idtype", idtype, &b_idtype);
   //fChain->SetBranchAddress("ichild", ichild, &b_ichild);
   Notify();
}

Bool_t ntuple_real::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef ntuple_real_cxx

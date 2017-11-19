//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov  9 16:16:09 2017 by ROOT version 6.06/08
// from TTree h1/Radphi
// found on file: /pnfs4/phys.uconn.edu/data/Gluex/radphi/omega_sim-9-2014/omega2p_gen-99-0.root
//////////////////////////////////////////////////////////

#ifndef mctuple_gen_h
#define mctuple_gen_h

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

// Headers needed by this particular selector


class mctuple_gen : public TSelector {
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
   TTree *tkinemc;
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
   TTree *tskimmc;
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
   TLorentzVector GetPpi0HelicityMC(double Ebeam);

   // Declaration of leaf types
   Int_t           runno;
   Int_t           eventno;
   Int_t           ismc;
   Int_t           ni;
   Int_t           type;
   Float_t         wgt;
   Int_t           kindi[30];   //[ni]
   Float_t         momi[30][4];   //[ni]
   Int_t           nv;
   Int_t           kindv[30];   //[nv]
   Float_t         origv[30][4];   //[nv]
   Int_t           nf;
   Int_t           kindf[30];   //[nf]
   Int_t           vertf[30];   //[nf]
   Float_t         momf[30][4];   //[nf]

   // List of branches
   TBranch        *b_runno;   //!
   TBranch        *b_eventno;   //!
   TBranch        *b_ismc;   //!
   TBranch        *b_ni;   //!
   TBranch        *b_type;   //!
   TBranch        *b_wgt;   //!
   TBranch        *b_kindi;   //!
   TBranch        *b_momi;   //!
   TBranch        *b_nv;   //!
   TBranch        *b_kindv;   //!
   TBranch        *b_origv;   //!
   TBranch        *b_nf;   //!
   TBranch        *b_kindf;   //!
   TBranch        *b_vertf;   //!
   TBranch        *b_momf;   //!

   mctuple_gen(TTree *tree =0);
   virtual ~mctuple_gen();
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

   ClassDef(mctuple_gen,0);

};

#endif

#ifdef mctuple_gen_cxx
mctuple_gen::mctuple_gen(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   TString rootfile_default("/pnfs4/phys.uconn.edu/data/Gluex/radphi/omega_sim-9-2014/omega2p_gen-99-0.root");
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(rootfile_default);
      if (!f || !f->IsOpen()) {
         f = new TFile(rootfile_default);
      }
      f->GetObject("h1",tree);
   }
   Init(tree);
}

mctuple_gen::~mctuple_gen()
{
   //if (!fChain) return;
   //delete fChain->GetCurrentFile();
}

void mctuple_gen::Init(TTree *tree)
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
   fChain->SetBranchAddress("ismc", &ismc, &b_ismc);
   fChain->SetBranchAddress("ni", &ni, &b_ni);
   fChain->SetBranchAddress("type", &type, &b_type);
   fChain->SetBranchAddress("wgt", &wgt, &b_wgt);
   fChain->SetBranchAddress("kindi", kindi, &b_kindi);
   fChain->SetBranchAddress("momi", momi, &b_momi);
   fChain->SetBranchAddress("nv", &nv, &b_nv);
   fChain->SetBranchAddress("kindv", kindv, &b_kindv);
   fChain->SetBranchAddress("origv", origv, &b_origv);
   fChain->SetBranchAddress("nf", &nf, &b_nf);
   fChain->SetBranchAddress("kindf", kindf, &b_kindf);
   fChain->SetBranchAddress("vertf", vertf, &b_vertf);
   fChain->SetBranchAddress("momf", momf, &b_momf);
   Notify();
}

Bool_t mctuple_gen::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef mctuple_gen_cxx

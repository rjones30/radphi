#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TString.h>

TString rootbase("/home/radphi/work1/");

TChain *chain_real()
{
   TChain *ch = new TChain("skim");
   for (int i=1; i < 21; ++i) {
      TString fname;
      fname.Form(rootbase + "ntuple_real_skim-%2.2d.root", i);
      ch->Add(fname);
   }
   return ch;
}

TChain *chain_mcgen()
{
   TChain *ch = new TChain("skimmc");
   for (int i=1; i < 13; ++i) {
      TString fname;
      fname.Form(rootbase + "mctuple_gen_skim-%d.root", i);
      ch->Add(fname);
   }
   return ch;
}

TChain *chain_mcsim()
{
   TChain *ch = new TChain("skim");
   TChain *chmc = new TChain("skimmc");
   for (int i=1; i < 13; ++i) {
      TString fname;
      fname.Form(rootbase + "mctuple_sim_skim-%d.root", i);
      ch->Add(fname);
      chmc->Add(fname);
   }
   ch->AddFriend(chmc);
   return ch;
}

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec 13 12:01:56 2021 by ROOT version 6.24/06
// from TTree tcl/Cluster data tree
// found on file: list_alpha_clusters.root
//////////////////////////////////////////////////////////

#ifndef TPX3Cluster_h
#define TPX3Cluster_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <algorithm>

// Header file for the classes stored in the TTree if any.

class TPX3Cluster {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           n;
   Short_t         x[5000];   //[n]
   Short_t         y[5000];   //[n]
   Long64_t        t[5000];   //[n]
   Short_t         e[5000];   //[n]
   Double_t        tof[5000];   //[n]

   // List of branches
   TBranch        *b_n;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_t;   //!
   TBranch        *b_e;   //!
   TBranch        *b_tof;   //!

   TPX3Cluster(TTree *tree=0);
   virtual ~TPX3Cluster();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TPX3Cluster_cxx
TPX3Cluster::TPX3Cluster(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/Project/2023/Napoli/fIaD_000000_clusters.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data/Project/2023/Napoli/fIaD_000000_clusters.root");
      }
      f->GetObject("tcl",tree);

   }
   Init(tree);
}

TPX3Cluster::~TPX3Cluster()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TPX3Cluster::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TPX3Cluster::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TPX3Cluster::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("n", &n, &b_n);
   fChain->SetBranchAddress("x", x, &b_x);
   fChain->SetBranchAddress("y", y, &b_y);
   fChain->SetBranchAddress("t", t, &b_t);
   fChain->SetBranchAddress("e", e, &b_e);
   fChain->SetBranchAddress("tof", tof, &b_tof);
   Notify();
}

Bool_t TPX3Cluster::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TPX3Cluster::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TPX3Cluster::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TPX3Cluster_cxx

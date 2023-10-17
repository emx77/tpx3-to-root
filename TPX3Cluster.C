#define TPX3Cluster_cxx
#include "TPX3Cluster.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>



void TPX3Cluster::Loop()
{

std::ofstream ofs ("cluster_data.txt", std::ofstream::out);

std::cout << "# ID" << ' ' << "size" << ' ' << "first_GToA " << ' ' << "SumToT " << ' ' << " mean_x " << ' ' << " mean_y " << std::endl;
std::cout << "#   " << ' ' << "    " << ' ' << "(1.5625 ns)" << ' ' << "(25 ns)" << ' ' << "(pixels)" << ' ' << "(pixels)" << std::endl;

ofs << "# ID" << ' ' << "size" << ' ' << "first_GToA " << ' ' << "SumToT " << ' ' << " mean_x " << ' ' << " mean_y " << std::endl;
ofs << "#   " << ' ' << "    " << ' ' << "(1.5625 ns)" << ' ' << "(25 ns)" << ' ' << "(pixels)" << ' ' << "(pixels)" << std::endl;

TH1F *hn = new TH1F("hn","",250,0,250);
//TH1F *ht = new TH1F("ht","",200,0,200);
//TH1F *ht = new TH1F("ht","",200,0,200);
TH2F *hxy = new TH2F("hxy","",256,0,256,256,0,256);


//   In a ROOT session, you can do:
//      root> .L TPX3Cluster.C
//      root> TPX3Cluster t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   cout << "nentries" << ' ' << nentries << endl;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   
      //nbytes += nb;
      if (jentry>=1000000) break;    
      // if (Cut(ientry) < 0) continue;

      Long64_t tmin = *std::min_element(t,t+n);
      Long64_t tmax = *std::max_element(t,t+n);

      Short_t xmin = *std::min_element(x,x+n);
      Short_t xmax = *std::max_element(x,x+n);

      Short_t ymin = *std::min_element(y,y+n);
      Short_t ymax = *std::max_element(y,y+n);

      Short_t esum = std::accumulate(e,e+n,0);

      hn->Fill(n);
      hxy->Fill(0.5*(xmax+xmin)-258,0.5*(ymax+ymin));	
      if (jentry%100000==0) {
          std::cout << jentry << ' ' << n << ' ' << tmin << ' ' << esum << ' ' << 0.5*(xmax+xmin)-258 << ' ' << 0.5*(ymax+ymin) << std::endl;
      }
      ofs << jentry << ' ' << n << ' ' << tmin << ' ' << esum << ' ' << 0.5*(xmax+xmin)-258 << ' ' << 0.5*(ymax+ymin) << std::endl;
      //if (n>7 && tmax-tmin < 200) { // ~0.3 uS
          //for (int i=0; i<n; i++) {
             //ht->Fill(t[i]-tmin);
          //}
      //}

      
       
        
   }

   // cluster area, compare with disk area

   // cluster min time
   // cluster mean time
   // cluster weighted mean time
   // cluster centre area time

   // pixel occurence

   // mean time deviation per pixel

   // time spread/accuracy within a cluster


   

   TCanvas *c1 = new TCanvas("c1","");
   c1->Divide(2,1);
   c1->cd(1)  ;
   hn->Draw();
   c1->cd(2);
   hxy->Draw("colz");	
}


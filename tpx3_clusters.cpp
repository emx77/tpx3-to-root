#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "TH2F.h"

#include "TROOT.h"

int clfind(int ihit, int clusid, int nsubset,
	   double *x, double *y, double *t, int *clusnr);

double t0find(Long64_t ntdc, double *tdc_time, double pixelhit_time);


int tpx3_clusters(string filename, long nhits=-1) {

   gROOT->Reset();

   TFile *f = new TFile(filename.c_str());
   if (!f || f->IsZombie()) {
     cout << "input file unavailable" << endl;
     return 1;
   }
   string ofile = filename.substr(0,filename.size()-5)+"_clusters.root";

   cout << "input file = " << filename << endl;
   cout << "output file = " << ofile << endl;

   TTree *t2 = (TTree*)f->Get("t2");
   long nentries=t2->GetEntries();
   cout << "Number of pixelhit entries: " << nentries << endl;
   
   TTree *ttdc = (TTree*)f->Get("ttrig");
   Long64_t ntrig = ttdc->GetEntries();
   cout << "Number of tdc entries: " << ntrig << endl; 
   
   
   // select rising edge TDCs
   ttdc->Draw("ts","type==0||type==2","goff");
   double *tdc = ttdc->GetV1();
   
   const Int_t kMaxPixel=300; // maximum allowed cluster size
   Int_t npix;
   Short_t xpix[kMaxPixel];
   Short_t ypix[kMaxPixel];
   Long64_t tpix[kMaxPixel];
   Short_t epix[kMaxPixel];
   //Float_t mx;
   //Float_t my;
   //Int_t etot;
   //Int_t dtpix[kMaxPixel];
   //Int_t tmin;
   Double_t tof[kMaxPixel];
   
   TFile *clFile = new TFile(ofile.c_str(),"recreate");
   TTree *tcl = new TTree("tcl","Cluster data tree");
   tcl->Branch("n",&npix,"n/I");
   tcl->Branch("x",xpix,"x[n]/S");
   tcl->Branch("y",ypix,"y[n]/S");
   tcl->Branch("t",tpix,"t[n]/L");
   tcl->Branch("e",epix,"e[n]/S");
   tcl->Branch("tof",tof,"tof[n]/D");
   //tcl->Branch("mx",&mx,"mx/F");
   //tcl->Branch("my",&my,"my/F");
   //tcl->Branch("etot",&etot,"etot/I");
   //tcl->Branch("tmin",&tmin,"tmin/I");
   //tcl->Branch("dtpix",dtpix,"dtpix/I");
   //tcl->Branch("dt",&dt,"dt/I");
      
   int stepsize=16384;
   //int nsteps=(nentries/stepsize)+1;
   //nsteps=100;
   
   int *clusnr = new int[stepsize];

   long nprocessed = 0;

   if (nhits<0 || nhits>nentries) nhits = nentries;
   cout << "number of hits to process: " << nhits << endl;
   
   int istep = 0;
   while (nprocessed<nhits) {
     if ( (nhits-nprocessed)<stepsize ) {
       stepsize = nhits-nprocessed;
     }
     // cout << "selecting tree data ... " << endl;
     int nsubset = t2->Draw("ypix:xpix:GToA:ToT", "","goff", stepsize, nprocessed);
     nprocessed+=nsubset;
     // cout << "number of entries in subset: " << nsubset << ' ' << istep << endl;
     
     //long nentries=t2->GetEntries();
 
     // cout << "cluster finding ... " << endl;
   
     double *x = t2->GetV1();
     double *y = t2->GetV2();
     double *t = t2->GetV3();
     double *e = t2->GetV4();
     
     // int *clusnr = new int[nsubset];

     
     
     
     // clusters IDs need to be resetted
     for (int i=0; i<nsubset; i++) {
       clusnr[i]=-1;
       // cout << x[i] << ' ' <<y[i] << ' ' << t[i] << endl;
     }
   
     int clusid=0;
     
     for (int i=0; i<nsubset; i++) {
       if (clusnr[i]==-1) {
         clfind(i, clusid, nsubset, x, y, t, clusnr);
         clusid++;
       }
     }
     cout << "subsection: " << istep << " number of clusters found: " << clusid << endl;
     istep++;
   //for (int i=0; i<nentries; i++) {
    // cout << x[i] << ' ' <<y[i] << ' ' << t[i] << ' ' << clusnr[i] << endl;
      
   //}

     // loop over found clusters
  
     for (int i=0; i<clusid; i++) {
       //cout << x[i] << ' ' <<y[i] << ' ' << t[i] << ' ' << clusnr[i] << endl;
       npix=0;
       //mx=0;
       //my=0;
       //etot=0;
       //tmin=0x3ffff;
       for (int j=0; j<nsubset; j++) {
         //cout << j << endl; 
         if (clusnr[j]==i && npix < kMaxPixel) {
            xpix[npix]=(Short_t)(x[j]+0.5);
            ypix[npix]=(Short_t)(y[j]+0.5);
            tpix[npix]=(Long64_t)(t[j]+0.5);
	   //if (tpix[npix]<tmin) {
	   //tmin=tpix[npix];
	   //}
            epix[npix]=(Short_t)(e[j]+0.5);
           //cout << npix << ' ' << e[j] << ' ' << epix[npix] << endl;
	   //mx+=xpix[npix];
	   //my+=ypix[npix];
	   //etot+=epix[npix];
            tof[npix] = (Double_t) (t[j]*1.5625E-9 - t0find(ntrig, tdc, t[j]*1.5625E-9) ) ;
            npix++;
        }
       }
       //if (npix!=0) {
       //mx = 0.5+mx/npix;
       //my = 0.5+my/npix;
	 //cout << etot << endl << endl;
       //}
       tcl->Fill();
     } // next cluster
   
   } // loop over all data subsets.

   delete[] clusnr;

   cout << "processed " << nprocessed << " entries (" << 100.*nprocessed/nentries << "%)" << endl;
   
   tcl->Write();

   clFile->Close();

   f->Close();


   cout << "end of program" << endl;
  return 0;
}

// cluster finding function
int clfind(int ihit, int clusid, int nsubset,
	   double *x, double *y, double *t, int *clusnr) {
  
  clusnr[ihit]=clusid;

  for (int i=0; i<nsubset; i++) {
    if (clusnr[i]<0) {
      if (  TMath::Abs(x[i]-x[ihit])<=1 &&  TMath::Abs(y[i]-y[ihit])<=1
	    && TMath::Abs(t[i]-t[ihit])<=640 ) {
	clfind(i,clusid,nsubset,x,y,t,clusnr);
      }
    }
  }
 
  return 0;
}

double t0find(Long64_t ntdc, double *tdc_time, double pixelhit_time) {
    Long64_t i=0;
    while (i<ntdc) {
      if (pixelhit_time<tdc_time[i]) {
          break;
      }
      i++;
    }
    return tdc_time[i-1];
};


#include <iostream>
#include  <iomanip>
#include <fstream>

#include "TROOT.h"

#include "TH1F.h"
#include "TH2F.h"

#include "TFile.h"
#include "TTree.h"

using namespace std;

int tpx3_to_root(string filename) {
  
  gROOT->Reset();

  int debug=0;

  TH1F *h1 = new TH1F("h1","erik",256,0,256);

  //TH2F *h2c0 = new TH2F("h2c0","",256,0,256,256,0,256);
  //TH2F *h2c1 = new TH2F("h2c1","",256,0,256,256,0,256);
  //TH2F *h2c2 = new TH2F("h2c2","",256,0,256,256,0,256);
  //TH2F *h2c3 = new TH2F("h2c3","",256,0,256,256,0,256);

  TH2F *h2quad = new TH2F("h2quad","",516,0,516,516,0,516);

    // frame
	  // chip
	  // pix x
	  // pix y
	  // ToT
	  // ToA
	  // FToA
	  // CToA
	  // SpiderTime

   cout << " input file = " << filename << endl;
   string ofile = filename.substr(0,filename.size()-4)+"root";
   cout << " output file = " << ofile << endl;
     
   TFile *f = new TFile(ofile.c_str(),"recreate");

   TTree *t2 = new TTree("t2","");
   
   UShort_t xpix, ypix, ToT, ToA, spidrTime;
   UChar_t chipnr,FToA; 
   UInt_t framenr;
   
   UInt_t CToA;

   // todo: use correct types
   
   t2->Branch("framenr",&framenr,"framenr/i");
   t2->Branch("chipnr",&chipnr,"chipnr/b");
   t2->Branch("xpix",&xpix,"xpix/s");
   t2->Branch("ypix",&ypix,"ypix/s");
   t2->Branch("ToT",&ToT,"ToT/s");
   t2->Branch("ToA",&ToA,"ToA/s");
   t2->Branch("FToA",&FToA,"FtoA/b");
   t2->Branch("CToA",&CToA,"CToA/i");  
   t2->Branch("SpidrTime",&spidrTime,"SpidrTime/s");  

   // cout << " opening file: " << filename << endl;

  ifstream infi(filename.c_str(), ifstream::binary);
  if (!infi) {
    cout << " file: " << filename << " not found." << endl;
  }
  else {

    const int hl=8; // length of header packet
    UChar_t *buffer = new UChar_t[hl];   
    const int dl=9000; // max length of data packet 
    unsigned long *databuffer = new unsigned long[dl];

    int nheaders=10000000;
    int maxhits =2000000000;
    
    // loop over the headers in the file

    // number of data packets
    long count=0;

    //number of pixel hits 
    long hitcount=0;
    
    // counts per chip
    long chipcount[4];

    // frame counter per chip
    int frame[4];
    
    for (int i=0;i<4;i++) {
      chipcount[i]=0;
      frame[i]=0;
    }
    
    while ((infi.good()) && count<nheaders) {
      if (count%10000==0) cout << "header: " << count << endl;
      count++;

      // reset buffer
      for (int i=0; i<hl; i++) {
        buffer[i]=0;
      }

      infi.read((char*)buffer, hl);

      chipnr = buffer[4];

      int mode = buffer[5];

      // 1b 90 00 03 33 58 50 54
      //     
      
      if (count==1) cout << "mode=" << mode << endl;
      
      int size = ((0xff & buffer[7]) <<8) | (0xff & buffer[6]);

      int pixdatasize = size/8; 
      
      ////for (int i=0; i<hl; i++) {
      ////  cout << ' ' << count << ' '  << buffer[i] << ' ' << (int)(unsigned char)buffer[i];
      ////}
      //cout << endl;

      //cout << " next packet size = " << size << endl;

      if (size>dl) {
        cout << "data size > " << dl << endl;
        break;
      }
         
      infi.read((char*)databuffer,size);
      
      for (int i=0; i<pixdatasize; i++) {
	//cout << hex << ' ' << databuffer[i] << dec << ' ';
	ULong64_t temp = databuffer[i];

	spidrTime = (UShort_t) (temp & 0xffff);

	//cout << hex << (temp>>56) << dec << ' ';
	
	int hdr = (int)(temp>>56);
	h1->Fill(hdr);

	//cout << hdr << endl;
	//if (hdr<150) {
	//  cout << count << ' ' << hex << temp << dec << endl;
	//}

	int h4 = temp>>48;

	if (h4==0x71bf) {
	  int chipID = (int) (temp >> 16) & 0xffff;
	  cout << count << ' ' << Form("EndOfCommand on %04x at %5d",chipID, spidrTime) << endl;
	}

	if (h4==0x71b0) {
	  int chipID = (int) (temp >> 16) & 0xffff;
	  cout << count << ' ' << Form("EndOfReadOut on %04x at %5d",chipID, spidrTime) << endl;
	  cout << (int)chipnr << ' ' << chipcount[chipnr] << endl;
	  chipcount[chipnr] = 0;
	  frame[chipnr]++;
	}

	if (h4==0x71ef) {
	  int chipID = (int) (temp >> 16) & 0xffff;
          cout << count << ' ' << Form("EndOfResetSequentialCommand on %04x at %5d",chipID, spidrTime) << endl;
	}
	
	int h2 = temp>>60;  
	if (h2==0x6) { 
	  cout << hex << temp << dec << endl;
	
	  long coarsetime = temp>>12 & 0xFFFFFFFF;
	
	  //cout << coarsetime*25e-9 << endl;
	  int tmpfine = (temp >> 5 ) & 0xF;   // 12 phases of 320 MHz clock in bits 5 to 8
	  tmpfine = ((tmpfine-1) << 9) / 12;     // subtract 1 (as fractions 0 to 11 are numbered 1 to 12), then shift by 9 positions to reduce the error coming from the integer division by 12 
	  int trigtime_fine = (temp & 0x0000000000000E00) | (tmpfine & 0x00000000000001FF);   // combine the 3 bits with a size of 3.125 ns with the rest of the fine time from the 12 clock phases
	  double time_unit=25./4096;
	  cout << setprecision(15) << (coarsetime*25E-9 + trigtime_fine*time_unit*1E-9) << endl;
        }
	if (temp>>60==0xb && hitcount < maxhits) {
	  if (hitcount%1000000==0) cout << "hit: " << hitcount << endl;
	  hitcount++;
	   
          chipcount[chipnr]++;
	  
	  // data driven pixel data
          // doublecolumn * 2
	  long dcol = (temp & 0x0FE0000000000000L) >> 52; //(16+28+9-1)
	  // superpixel * 4
	  long spix = (temp & 0x001F800000000000L) >> 45; //(16+28+3-2)
	  // pixel
	  long pix = (temp & 0x0000700000000000L) >> 44; //(16+28)
	  int x = (int) (dcol + pix / 4);
	  int y = (int) (spix + (pix & 0x3));
	  //int yx = y * 256 + x;
	  if (debug) cout << count << ' ' << hitcount << ' ' << chipnr << ' ' << x << ' ' << y << endl;

	  framenr=frame[chipnr];

          ToA  = (UShort_t) ((temp >> (16 + 14)) & 0x3fff);
	  ToT  = (UShort_t) ((temp >> (16 + 4)) & 0x3ff);
          FToA =  (UChar_t) ((temp >> 16) & 0xf);
	  CToA = (ToA << 4) | (~FToA & 0xf);
	  
	  if (chipnr==0) {
            // h2c0->Fill(x,y);
	    xpix=x+260;
	    ypix=y; 
	  }
	  if (chipnr==1) {
	    // h2c1->Fill(255-x,255-y);
	    xpix=255-x+260;
	    ypix=255-y+260;
	  }
	  if (chipnr==2) {
	     // h2c2->Fill(255-x,255-y);
	     xpix=255-x;
	     ypix=255-y+260;
	  }
	  if (chipnr==3) {
            // h2c3->Fill(x,y);
	    xpix=x;
	    ypix=y;
	  }

	  h2quad->Fill(xpix,ypix);
                      
          t2->Fill();
	 
	} // ==> if pixel data
			
      } // end of data packet loop
    
    } // end of headers loop

    cout << count << " data packets found " << endl;
    cout << hitcount << " pixel hits found " << endl;

    delete buffer;
    delete databuffer;

  } // file exists
  infi.close();
  
  h1->Write();
  h2quad->Write();  

  t2->Write();
  f->Close();

  delete h1;
  delete h2quad;
  
  return 0;
}

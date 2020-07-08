#include <iostream>
#include <iomanip>
#include <fstream>

#include "TROOT.h"

#include "TH1F.h"
#include "TH2F.h"

#include "TFile.h"
#include "TTree.h"

#include "TMath.h"
#include "TTimeStamp.h"

using namespace std;

int tpx3_to_root(string filename, unsigned long nrawpixelhits=0, unsigned first_frame=0) {
    
    gROOT->Reset();
    
    TTimeStamp ts_start;

    TTimeStamp ts_sort_start;
   
    TTimeStamp ts_sort_end;
    
    int debug=0;
    
    int sort=0;

    TH1F *h1 = new TH1F("h1","erik",256,0,256);

    //TH2F *h2c0 = new TH2F("h2c0","",256,0,256,256,0,256);
    //TH2F *h2c1 = new TH2F("h2c1","",256,0,256,256,0,256);
    //TH2F *h2c2 = new TH2F("h2c2","",256,0,256,256,0,256);
    //TH2F *h2c3 = new TH2F("h2c3","",256,0,256,256,0,256);

    TH1F *h11 = new TH1F("h11","",1125,0,1125);
    TH1F *h12 = new TH1F("h12","dt",2000,-1.5625e-3,1.5625e-3);
    TH1F *h13 = new TH1F("h13","dt",3000,-2.0e-3,1.0e-3);
    TH1F *h14 = new TH1F("h14","dt",3000,-2.0e-3,1.0e-3);
    

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
    // cout << " number of raw pixel hits = " << nrawpixelhits << endl;
    cout << " first frame number = " << first_frame << endl;

    string ofile = filename.substr(0,filename.size()-4)+"root";
    cout << " output file = " << ofile << endl;
    
    TFile *f = new TFile(ofile.c_str(),"recreate");

    TTree *t2 = new TTree("t2","");

    UShort_t xpix, ypix, ToT, ToA, spidrTime;
    UChar_t chipnr,FToA; 
    UInt_t framenr;

    UInt_t CToA;

    // todo: use correct types
    // todo: CToA is double information
    // todo: store tdc stamps
    // todo: check for roll over in data conversion
    // todo: reduce hard coded values

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
    
    // todo: open multiple .tpx3 files

    ifstream infi(filename.c_str(), ifstream::binary);
    if (!infi) {
        cout << " file: " << filename << " not found." << endl;
    }
    else {
        const int hl=8; // length of header packet
        UChar_t *buffer = new UChar_t[hl];   
        const int dl=9000; // max length of data packet 
        unsigned long *databuffer = new unsigned long[dl];

        int nheaders= 10000000;
        unsigned long maxhits =2000000000;
        if (nrawpixelhits!=0) maxhits=nrawpixelhits;

        unsigned long Timer_LSB32 = 0;
        unsigned long Timer_MSB16 = 0;
        unsigned long long timemaster = 0;
    
        // loop over the headers in the file

        // number of data packets
        long count=0;

        //number of pixel hits 
        unsigned long hitcount=0;
    
        // counts per chip
        long chipcount[4];

        // frame counter per chip
        int frame[4];
    
        for (int i=0;i<4;i++) {
            chipcount[i]=0;
            frame[i]=0;
        }

        double tdc_time = -1;
        int trigcnt = -1;
        
        Double_t t_min=TMath::Power(2,34);
        Double_t t_max=0;
        
        Double_t t_prevhit = 0;
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

            h11->Fill(pixdatasize);
    
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
            
            //if (chipnr!=1) continue;
            
            //UShort_t prev_time = 0;
            
            for (int i=0; i<pixdatasize; i++) {
                //cout << hex << ' ' << databuffer[i] << dec << ' ';
                ULong64_t temp = databuffer[i];

                spidrTime = (UShort_t) (temp & 0xffff);
                 //if (i>0) {
                //   h12->Fill(spidrTime-prev_time);    
                //}
                // prev_time=spidrTime;

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
                    int h3 = temp>>56;
                    int tdc_nr=-1; // can be 1 or 2
                    int edge_type=-1; // 0 for rising edge, 1 for falling edge	
            
                    if (h3 == 0x6f) {
                        tdc_nr=1;
                        edge_type=0;		
                    }
                    else if (h3 == 0x6a) {
                        tdc_nr=1;
                        edge_type=1;		
                    }
                    else if (h3 == 0x6e) {
                        tdc_nr=2;
                        edge_type=0;		
                    }
                    else if (h3 == 0x6b) {
                        tdc_nr=2;
                        edge_type=1;		
                    } 	
                    trigcnt = (int) (temp>>44) & 0xfff;   

                    long coarsetime = temp>>12 & 0xFFFFFFFF;	
        
                    //cout << coarsetime*25e-9 << endl;
                    int tmpfine = (temp >> 5 ) & 0xF;   // 12 phases of 320 MHz clock in bits 5 to 8
                    tmpfine = ((tmpfine-1) << 9) / 12;     // subtract 1 (as fractions 0 to 11 are numbered 1 to 12), then shift by 9 positions to reduce the error coming from the integer division by 12 
                    int trigtime_fine = (temp & 0x0000000000000E00) | (tmpfine & 0x00000000000001FF);   // combine the 3 bits with a size of 3.125 ns with the rest of the fine time from the 12 clock phases
                    double time_unit=25./4096;
                    tdc_time = (coarsetime*25E-9 + trigtime_fine*time_unit*1E-9);
                    if (count<20) { 
                        cout << count << ' ' << trigcnt << ' ' << hex << temp << dec << endl; 
                        cout << "tdc_nr: " << tdc_nr << " edge_type: " << edge_type << " tdc_time: " << setprecision(15) <<  tdc_time << endl;
                    }
            
                }
                if (h2==0x4) { 

                    //cout << hex << temp << dec << endl;
            
                    //spidr time is already decoded.

                    if (((temp >> 56) & 0xF) == 0x4) {
                        Timer_LSB32 = (temp >> 16) & 0xFFFFFFFF;
                    } 
                    else if (((temp >> 56) & 0xF) == 0x5) {
                        Timer_MSB16 = (temp >> 16) & 0xFFFF;
                        timemaster = Timer_MSB16;						
                        timemaster = (timemaster << 32) & 0xFFFF00000000;
                        timemaster = timemaster | Timer_LSB32;
                        //todo: check for jumps for specific settings
                        //int diff = (spidrTime >> 14) - ((Timer_LSB32 >> 28) & 0x3);
                        //if ((spidrTime >> 14) == ((Timer_LSB32 >> 28) & 0x3)) { 
                        //} 
                        //else { 
                        //  Timer_MSB16 = Timer_MSB16 - diff;
                        //}
                        cout << "Global time:  " << timemaster *25e-9 << endl; //converted Global timestamp
                    }
        
                }
        
                if (h2==0xb && hitcount < maxhits) {
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
                    
                    Double_t t_hit = (Double_t)(409.6e-6)*spidrTime+(Double_t)(1.5625e-9)*((Int_t)CToA-15); 
                    if (t_hit>t_max) t_max = t_hit;
                    if (t_hit<t_min) t_min = t_hit;
                    if (t_hit>7.0) {cout << "hc: " << hitcount << endl;}
                
                    if (hitcount%10000==0 || hitcount%10000==1) {
                        
                        
                        cout << "hitcount: " << hitcount 
                             << " pixel hit time " 
                             << setprecision(12) 
                             << t_hit << " " 
                             << (spidrTime*262144+ToA*16-FToA) << endl;
                        cout << " global pixel hit time " << setprecision(12) 
                             << (timemaster >> 30)*2*13.4217728 + spidrTime*409.6e-6 + CToA*1.5625e-9 << endl;
                        cout << " tdc trig cnt: " << trigcnt << " last tdc event time: " << tdc_time << endl; 
                    }
                    
                    if (hitcount>1) {
                        h12->Fill(t_hit-t_prevhit);    
                    }
                    t_prevhit = t_hit;

                    //todo: check for jumps for specific settings

                    //todo: check alignment of pixel hits and tdc time stamps with global time stamp 
                    
                    //int diff = ( spidrTime >> 14 ) - ( (Timer_LSB32 >> 28) & 0x3 );
                    //if (diff!=0) {
                    //  cout << diff << ' ' << hex << ((Timer_LSB32 >> 28) & 0x3) << " " << (spidrTime>>14) << dec << endl;
                    //}
                    //else {                                                                                                                                                    
                    //	//cout << (timemaster >> 29)*13.4217728 +(spidrTime&0x7FFF)*409e-6 + CToA*1.5625e-9 << endl;
                    //}
                    //if ((spidrTime >> 14) == ((Timer_LSB32 >> 28) & 0x3)) { 
                    //} 
                    //else { 
                    //  Timer_MSB16 = Timer_MSB16 - diff;
                    //}

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

                    if (framenr >= first_frame) {
                        h2quad->Fill(xpix,ypix);           
                        t2->Fill();
                    }
    
                } // ==> if pixel data
            
            } // end of data packet loop
    
        } // end of headers loop

        cout << count << " data packets found " << endl;
        cout << hitcount << " pixel hits found " << endl;
        
        cout << "ToA of first hit (s): " << t_min << endl;
        cout << "ToA of last hit (s): "  << t_max << endl;
        if (t_max>t_min) {
           cout << "pixelhits/s (MHz) : " << hitcount/(t_max-t_min)/1e6 << endl;
        }
        
        for (ULong_t i=0; i<hitcount; i++) {
        Double_t t_expected = t_min+ (t_max-t_min)*i/hitcount;
		t2->GetEntry(i);
        Double_t t_hit = (Double_t)(1.5625e-9)*(spidrTime*262144+ToA*16-FToA);
        //cout << t_expected << ' ' << t_hit << endl;
        h13->Fill(t_hit - t_expected);
        }
    
        TTree *tsorted = (TTree*)t2->CloneTree(0);
        
        ts_sort_start.Set();
    
        if (sort) {
            cout << hitcount << " sorting " << endl;
            t2->SetEstimate(hitcount);
            t2->Draw("SpidrTime*262144+ToA*16-FToA","","goff");
            ULong_t *index = new ULong_t[hitcount];
            TMath::Sort(hitcount, t2->GetV1(),index,kFALSE);

            for (ULong_t i=0;i<hitcount;i++) {
                t2->GetEntry(index[i]);
                tsorted->Fill();
                Double_t t_expected = t_min+ (t_max-t_min)*i/hitcount;
                Double_t t_hit = (Double_t)(1.5625e-9)*(spidrTime*262144+ToA*16-FToA);
                h14->Fill(t_hit - t_expected);
            }
            
            ts_sort_end.Set();
        
            delete [] index;
        
            tsorted->Write();
        }
        else {
            t2->Write();
        }
	

        delete[] buffer;
        delete[] databuffer;

    } // file exists
    infi.close();

    h1->Write();
    h11->Write(); 
    h12->Write();
    h2quad->Write();  
    h13->Write();
    h14->Write();
     
    //ULong_t nentries = t2->GetEntries();

    //t2->SetEstimate(nentries);
    
    //t2->Draw("SpidrTime*262144+ToA*16-FToA","","goff");
    
    
    //ULong_t *index = new ULong_t[nentries];

    //cout << nentries << " sorting " << endl;
    
	//TMath::Sort(nentries, t2->GetV1(),index,kFALSE);
	
	
	//Create an empty clone of the original tree
    
	//TTree *tsorted = (TTree*)t2->CloneTree(0);
	//for (ULong_t i=0;i<nentries;i++) {
	//	t2->GetEntry(index[i]);
	//	tsorted->Fill();
	//}
	
	
	//delete [] index;
    
   
    //t2->Write();
    f->Close();

    delete h1;
    delete h2quad;
    delete h11;
    delete h12;
    delete h13;
    delete h14;
    
    TTimeStamp ts_end;
    
    cout << "total processing time: " << ts_end-ts_start << " (s)" <<endl;
    if (sort)  cout << " sorting time: " << ts_sort_end - ts_sort_start << " (s)" <<endl;

    return 0;
}

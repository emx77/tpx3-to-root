#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "TROOT.h"

#include "TH1F.h"
#include "TH2F.h"

#include "TFile.h"
#include "TTree.h"

#include "TMath.h"

using namespace std;

int tpx3_to_root(string filename, unsigned long nrawpixelhits=0) {

    int debug=0;
    
    gROOT->Reset();

    TH1F *h1 = new TH1F("h1","erik",256,0,256);

    TH2F *h2quad = new TH2F("h2quad","",516,0,516,516,0,516);

    cout << " input file = " << filename << endl;
    
    int pos = filename.rfind('.');
    
    string extension = filename.substr(pos+1,filename.size());
    
    vector<string> tpx3_file;
    
    if (extension=="txt") { // list of files
        ifstream f(filename);
        string tmp_string;
        while (!f.eof() ) {
          f >> tmp_string;
          // cout << tmp_string << endl;
          if (f.eof()) break;
          tpx3_file.push_back(tmp_string); 
        }
    }
    else {
        tpx3_file.push_back(filename);
    }
    
    
    int nfiles = tpx3_file.size();
    
    cout << "nr of files to process: " << tpx3_file.size() << endl;
    
    string ofile = filename.substr(0,pos+1)+"root";
    cout << " output file = " << ofile << endl;
    TFile *f = new TFile(ofile.c_str(),"recreate");

    TTree *t2 = new TTree("t2","");
   
    UShort_t xpix, ypix, ToT, ToA, spidrTime;
    UChar_t chipnr,FToA; 
    UInt_t framenr;
   
    Int_t CToA;
    Long_t GToA;

    t2->Branch("framenr",&framenr,"framenr/i");
    t2->Branch("chipnr",&chipnr,"chipnr/b");
    t2->Branch("xpix",&xpix,"xpix/s");
    t2->Branch("ypix",&ypix,"ypix/s");
    t2->Branch("ToT",&ToT,"ToT/s");
    t2->Branch("ToA",&ToA,"ToA/s");
    t2->Branch("FToA",&FToA,"FtoA/b");
    // CToA now obsolete because of glocal ToA which inlcudes the spidrTime rollovers
    // roll over detection requires at least one pixelhit every ~ 3s. 
    // t2->Branch("CToA",&CToA,"CToA/I");  
    t2->Branch("GToA",&GToA,"GToA/L");
    t2->Branch("SpidrTime",&spidrTime,"SpidrTime/s");  

    // A tree fr TDC data

    Double_t tdc_ts; // timestamp is s
    UChar_t tdc_type;
    UInt_t tdc_nr;
    
    TTree *ttdc = new TTree("ttrig","");
   	ttdc->Branch("ts",&tdc_ts,"ts/D");
    ttdc->Branch("type",&tdc_type,"type/b");
    ttdc->Branch("nr",&tdc_nr,"nr/i"); 
    
    
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
    double prev_tdc_time = -1;
    
    unsigned long tdc_1r = 0;
    unsigned long tdc_1f = 0;
    unsigned long tdc_2r = 0;
    unsigned long tdc_2f = 0;  

    double maxTDC = 3.125E-9*TMath::Power(2,35);
    int ro_tdc_count = 0;   
    //int ro_tdc_state = 0;   
    
    
	// variables for roll over detection and correction
    int ro_count=0;
    int ro_state=0;
    long maxGToA = TMath::Power(2,34);
    int late_hit=0;
    
    for (int ifile=0; ifile<nfiles; ifile++) { 
        cout << " opening file: " << tpx3_file[ifile] << endl;
        ifstream infi(tpx3_file[ifile].c_str(), ifstream::binary);
        if (!infi) {
            cout << " file: " << filename << " not found." << endl;
        }
        else {

           
        while ((infi.good()) && count<nheaders &&hitcount<maxhits) {
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
                cout << size << " data size > " << dl << endl;
                break;
            }
            
            infi.read((char*)databuffer,size);
            
            // int tmp_max=0;
            // int tmp_min=100000;
            int npixhits=0;
            
	    	

            for (int i=0; i<pixdatasize; i++) {
                //cout << hex << ' ' << databuffer[i] << dec << ' ';
                ULong64_t temp = databuffer[i];
                
                spidrTime = (UShort_t) (temp & 0xffff);
                
                //if (spidrTime>tmp_max) tmp_max=spidrTime;
                //if (spidrTime<tmp_min) tmp_min=spidrTime;
               
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
                    if (debug) cout << count << ' ' << Form("EndOfCommand on %04x at %5d",chipID, spidrTime) << endl;
                }
                
                if (h4==0x71b0) {
                    int chipID = (int) (temp >> 16) & 0xffff;
                    if (debug) cout << count << ' ' << Form("EndOfReadOut on %04x at %5d",chipID, spidrTime) << endl;
                    if (debug) cout << (int)chipnr << ' ' << chipcount[chipnr] << endl;
                    chipcount[chipnr] = 0;
                    frame[chipnr]++;
                }
                
                if (h4==0x71ef) {
                    int chipID = (int) (temp >> 16) & 0xffff;
                    if (debug)  cout << count << ' ' << Form("EndOfResetSequentialCommand on %04x at %5d",chipID, spidrTime) << endl;
                }
                
                int h2 = temp>>60;  
                if (h2==0x6) {
                    int h3 = temp>>56;
                    int tdc_chan=-1; // can be 1 or 2
                    int edge_type=-1; // 0 for rising edge, 1 for falling edge	
                    
                    if (h3 == 0x6f) {
                        tdc_chan=1;
                        edge_type=0;
                        tdc_1r++; 
                        tdc_type = 0;
                    }
                    else if (h3 == 0x6a) {
                        tdc_chan=1;
                        edge_type=1;
                        tdc_1f++;
                        tdc_type = 1; 
                    }
                    else if (h3 == 0x6e) {
                        tdc_chan=2;
                        edge_type=0;
                        tdc_2r++;
                        tdc_type = 2; 
                    }
                    else if (h3 == 0x6b) {
                        tdc_chan=2;
                        edge_type=1;
                        tdc_2f++;
                        tdc_type = 3;
                    } 	
                    trigcnt = (int) (temp>>44) & 0xfff;   
                    tdc_nr = trigcnt;
                    
                    
                    prev_tdc_time = tdc_time;
                    
		    // 32 bits 
                    long coarsetime = temp>>12 & 0xFFFFFFFF;	
                    
                    //cout << coarsetime*25e-9 << endl;
                    int tmpfine = (temp >> 5 ) & 0xF;   // 12 phases of 320 MHz clock in bits 5 to 8
                    tmpfine = ((tmpfine-1) << 9) / 12;     // subtract 1 (as fractions 0 to 11 are numbered 1 to 12), then shift by 9 positions to reduce the error coming from the integer division by 12 
                    int trigtime_fine = (temp & 0x0000000000000E00) | (tmpfine & 0x00000000000001FF);   // combine the 3 bits with a size of 3.125 ns with the rest of the fine time from the 12 clock phases
                    double time_unit=25./4096;
                    tdc_time = ((double)coarsetime*25E-9 + trigtime_fine*time_unit*1E-9);
                    //if (count<20) { 
                    //    cout << count << ' ' << trigcnt << ' ' << hex << temp << dec << endl; 
                    //    cout << "tdc_chan: " << tdc_chan << " edge_type: " << edge_type << " tdc_time: " << setprecision(15) <<  tdc_time << endl;
                    //}
                    
                    if (tdc_time<prev_tdc_time) {
                        ro_tdc_count+=1;
                    }
                    
                    tdc_ts = tdc_time+(ro_tdc_count)*maxTDC;
                    ttdc->Fill();
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
                        cout << "Global time:  " << setprecision(15) << timemaster *25e-9 << endl; //converted Global timestamp
                    }
                    
                }
                
                if (h2==0xb && hitcount < maxhits) {
                    if (hitcount%1000000==0) cout << "hit: " << hitcount << endl;
                    hitcount++;
                    npixhits++;
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
                    // CToA calculation to keep CToA>=0
                    //CToA = (ToA << 4) | (~FToA & 0xf);
                    // uncorrected CToA calculation 
                    CToA = (Int_t) (ToA << 4) - FToA;
                                                  
                    GToA = ((Long_t(spidrTime)) << 18 ) + Long_t(CToA);
                    late_hit = 0;
                    if (1.0*GToA>0.95*maxGToA && ro_state==0) { 
                         ro_state=1;
                         //cout << "ro1 ";
                    }
                    if (ro_state==1 && 1.0*GToA<0.05*maxGToA) {
                        ro_state=2;
                        //cout << "ro2 "; 
                        ro_count++;
                    }
                    if (ro_state==2) {
                            if (1.0*GToA>0.95*maxGToA) {
                                late_hit=1;
                                //cout << "lh " << ro_count << ' ';
                            }
                            else if (1.0*GToA>0.05*maxGToA) {
                                ro_state=0;
                                //cout << "ro0 ";
                            }                                                        
                    }
                    GToA = GToA + (ro_count-late_hit)*maxGToA;
                                 
                    
                    //if (hitcount%10000==0 || hitcount%10000==1) {
                    //  cout << "hitcount: " << hitcount << " pixel hit time " << spidrTime*409.6e-6+CToA*1.5625e-9 <<  endl;
                    //   cout << " global pixel hit time " << (timemaster >> 30)*2*13.4217728 + spidrTime*409.6e-6 + CToA*1.5625e-9 << endl;
                    //  cout << " tdc trig cnt: " << trigcnt << " last tdc event time: " << tdc_time << endl; 
                    //}
                    
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
                        xpix=x+260;
                        ypix=y; 
                    }
                    if (chipnr==1) {
                        xpix=255-x+260;
                        ypix=255-y+260;
                    }
                    if (chipnr==2) {
                        xpix=255-x;
                        ypix=255-y+260;
                    }
                    if (chipnr==3) {
                        xpix=x;
                        ypix=y;
                    }
                    
                    h2quad->Fill(xpix,ypix);
                    
                    t2->Fill();
                      
                } // ==> if packet is a pixelhit
                
            } // next datapacket entry 
        
            // check for roll overs inside the datapacket
            // if (npixhits>0) cout << "packet " << count << " " << npixhits << ' ' << tmp_min << ' ' << tmp_max << endl;
        
        
        } // next datapacket
        
        
        
        
        
        } // file exists and good
    
        infi.close();
  
    } // next file in list
    
    cout << count << " data packets found " << endl;
    cout << hitcount << " pixel hits found " << endl;
        
    cout << "tdc: " << tdc_1r << ' ' << tdc_1f << ' ' << tdc_2r << ' ' << tdc_2f << ' ' <<endl;  

    delete[] buffer;
    delete[] databuffer;
    
    h1->Write();
    h2quad->Write();  

    t2->Write();
    ttdc->Write();
    f->Close();

    delete h1;
    delete h2quad;
      
    return 0;
}

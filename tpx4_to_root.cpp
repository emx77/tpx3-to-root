
#include <iostream>
#include <fstream>
#include <vector>

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

int tpx4_to_root(string filename, unsigned long nrawpixelhits=0) {

    gROOT->Reset();

    ifstream infi(filename.c_str(), ifstream::binary);

    if (!infi) {
        cout << " file: " << filename << " not found." << endl;
        return 1; 
    }

    int pos = filename.rfind('.');

    string extension = filename.substr(pos+1,filename.size());
    
    vector<string> tpx4_file;
    if (extension=="txt") { // list of files
        ifstream f(filename);
        string tmp_string;
        while (!f.eof() ) {
          f >> tmp_string;
          // cout << tmp_string << endl;
          if (f.eof()) break;
          tpx4_file.push_back(tmp_string); 
        }
    }
    else {
        tpx4_file.push_back(filename);
    } 
    int nfiles = tpx4_file.size();
    
    cout << "nr of files to process: " << tpx4_file.size() << endl;


    string ofile = filename.substr(0,pos+1)+"root";
    cout << " output file = " << ofile << endl;
    TFile *f = new TFile(ofile.c_str(),"recreate");

    TH1F *h1hb = new TH1F("h1hb","h1hb",11000,-0.5,10.5);  


    TH1F *h1 = new TH1F("h1","eoc",256,0,256);  
    TH1F *h1eoc2 = new TH1F("h1eoc2","eoc2",256,0,256);  

    TH1F *h1tot = new TH1F("h1tot","tot",256,0,2048);  
    TH1F *h1toa = new TH1F("h1toa","toa",65536,0,65536);  
    TH1F *h1pileup = new TH1F("h1pileup","pileup",2,0,2);  

    TH1F *h1ufToA_start = new TH1F("h1ufToA_start","h1ufToA_start",16,0,16);  
    TH1F *h1ufToA_stop = new TH1F("h1ufToA_stop","h1ufToA_stop",16,0,16);  
    
    TH1F *h1fToA_rise = new TH1F("h1fToA_rise","h1fToA_rise",32,0,32);  
    TH1F *h1fToA_fall = new TH1F("h1fToA_fall","h1fToA_fall",32,0,32);  
   
    TH2F *h2 = new TH2F("h2","xy",448,0,448,512,0,512);
    TH2F *h2tot = new TH2F("h2tot","xytot",448,0,448,512,0,512);
    TH2F *h2toa = new TH2F("h2toa","xytoa",448,0,448,512,0,512);


    UChar_t top; 

    UShort_t Col, Row;
    UInt_t ToT, ToA;
    UChar_t Pileup;
    UChar_t fToA_rise, fToA_fall; 
    UChar_t ufToA_start, ufToA_stop;

    ULong64_t heartbeat;
    ULong64_t last_hb;  

    ULong64_t last_heartbeat_top = 0;   
    ULong64_t last_heartbeat_bottom = 0;
      

    TTree *t2 = new TTree("t2",""); 
    //t2->Branch("framenr",&framenr,"framenr/i");
    //t2->Branch("chipnr",&chipnr,"chipnr/b");
    t2->Branch("xpix",&Col,"xpix/s");
    t2->Branch("ypix",&Row,"ypix/s");
    t2->Branch("ToT",&ToT,"ToT/s");
    t2->Branch("ToA",&ToA,"ToA/s");
    t2->Branch("Pilepup",&Pileup,"Pileup/b");
    t2->Branch("fToA_rise",&fToA_rise,"fToA_rise/b");
    t2->Branch("fToA_fall",&fToA_fall,"fToA_fall/b");

    t2->Branch("ufToA_start",&ufToA_start,"ufToA_start/b");
    t2->Branch("ufToA_stop",&ufToA_stop,"ufToA_stop/b");

    t2->Branch("last_hb",&last_hb,"last_hb/L"); 

    TTree *hb = new TTree("hb",""); 
    hb->Branch("hb",&heartbeat,"hb/L"); 
    hb->Branch("top",&top,"top/b"); 
     
    
    //t2->Branch("FToA",&FToA,"FtoA/b");
    // CToA now obsolete because of glocal ToA which inlcudes the spidrTime rollovers
    // roll over detection requires at least one pixelhit every ~ 3s. 
    // t2->Branch("CToA",&CToA,"CToA/I");  
    //t2->Branch("GToA",&GToA,"GToA/L");
    //t2->Branch("SpidrTime",&spidrTime,"SpidrTime/s");      


   
    const Int_t dpl=8; // length of data packet
    
    ULong64_t *data_packet = new ULong64_t; 
    ULong64_t hb_packet = 0; 
    // temperary char buffer 
    Char_t *buffer = new Char_t[dpl];

      
    ULong64_t npixelhits = 0; 
    ULong64_t ntimestamps = 0;
    ULong64_t ntpx4markers = 0;
    

    ULong64_t maxhits =2000000000;
    if (nrawpixelhits!=0) maxhits=nrawpixelhits;


    for (int ifile=0; ifile<nfiles; ifile++) { 
        cout << " opening file: " << tpx4_file[ifile] << endl;
        ifstream infi(tpx4_file[ifile].c_str(), ifstream::binary);
        if (!infi) {
            cout << " file: " << filename << " not found." << endl;
        }
        else {
            ULong64_t i=0; 
            infi.seekg (0, infi.end);
            ULong64_t infi_length = infi.tellg();
            infi.seekg (0, infi.beg);
            ULong64_t g = infi.tellg();

            while (i<maxhits && g<infi_length) {
              
                infi.read((Char_t*)data_packet, dpl); 
                i++; 
                g = infi.tellg();
                UInt_t eoc = ((*data_packet)>>55) & 0xFF;
                top=(*data_packet>>63)&0x1;
                //if (top) { continue; }

                h1->Fill(eoc);

                if (eoc==0xE0) {
                    ntimestamps++;
                    hb_packet = *data_packet;
                    heartbeat = (*data_packet & 0xFFFFFFFFFFFF);
                    cout << ntimestamps << ' ' << (int)top << ' ' << heartbeat << ' ' << heartbeat*25e-9 <<  endl;
                    if (top) last_heartbeat_top = heartbeat;
                    else last_heartbeat_bottom = heartbeat; 
                    h1hb->Fill(heartbeat*25e-9);
                    hb->Fill(); 
                    continue;  
                }
        
                buffer = (Char_t*) data_packet;
                TString s((Char_t*)buffer); 
                //cout << dec << g << ' ' << s << '-' << s.SubString(0,8) << endl;
         
                if (s.Contains("TPX4")) {
                    ntpx4markers++;
                    continue;
                } 

                // there are additional packets with high eoc values;
                if (eoc>223) {
                    cout << (int)top << " eoc:" << eoc << ' ' << *data_packet << ' ' << (((*data_packet)>>52) & 0x7) << ' ' << ((*data_packet) &0x7ffffffffffff) << endl; 
                    continue;
                } 
                // 

                //if (s.Contains("TPX4")) {
                //    cout << dec << g << ' ' << hex << *data_packet << " (TPX4 packet) " << endl;
                //}

                ToA = (*data_packet>>30) & 0xFFFF;
                // grey to binary count conversion
                uint16_t temp = ToA;
                uint16_t inv = 0;
                while (temp != 0) {
	                inv ^= temp;
	                temp >>= 1;
                }
                ToA = inv;

                //ULong64_t extToa = (hb_packet & ~0xFFFFl) | toa;
       
                //uint16_t val = ToA;




                //val ^= val >> 8;
                //val ^= val >> 4;
                //val ^= val >> 2;
                //val ^= val >> 1;
                //ToA=val;
          
                // cout << ToA << endl;
                //ToA = inv;
                ToT = (*data_packet>>1) & 0x7FF;
                Pileup = (*data_packet) & 0x1;    
   
                ULong64_t addr=(*data_packet>>46) & 0x3ffff;

                UInt_t Pix=addr&0x1F;
                UInt_t Sp=(addr>>5)&0xF;
                UInt_t EoC=(addr>>9)&0xFF;
                UInt_t Top=(addr>>17)&0x1;

                h1eoc2->Fill(EoC); 

                if (Top) {
                    EoC = 223-EoC;
                    Sp=15-Sp;
                    Pix=31-Pix;
                }
                Col = (UShort_t) (2*EoC + (Int_t)(Pix%8/4));
                Row = (UShort_t) (Pix%4 + (Int_t)(Pix/8)*4 + Sp*16 + Top*256);


                if (eoc>223) { 
                    cout << dec << Top << ' ' << g << ' ' << hex << *data_packet << dec << " eoc " << eoc << ' ' <<  endl;
                    if (Row>447) {
                    cout << dec << g << ' ' << hex << *data_packet << dec << ' ' << eoc << ' '<< Col << ' ' << Row << ' ' << ToA << ' ' << ToT << ' ' << Pileup << endl;
                    cout << Top << ' ' << EoC << ' ' << Sp << ' ' << Pix << endl;
                    } 
                } 
            
                //cout << "    " << Top << ' ' << ToA << ' ' << ToA*25e-9 << endl;    
        
                ufToA_start = (*data_packet>>26) & 0xF;
                ufToA_stop  = (*data_packet>>22) & 0xF;
                fToA_rise   = (*data_packet>>17) & 0x1F;
                fToA_fall   = (*data_packet>>12) & 0x1F;

                UInt_t uftoa_list[8] = {15, 14, 12, 8, 0, 1, 3, 7};
                UInt_t j=0; 
                while (j<8) {
                    if (uftoa_list[j]==ufToA_start) {
                        ufToA_start = j;
                        break;
                    } 
                    j++;
                    if (j==8) cout << " wrong ufToA_start value: " << ufToA_start << endl; 
                } 
                j=0;
                while (j<8) {
                    if (uftoa_list[j]==ufToA_stop) {
                        ufToA_stop = j;
                        break;
                    } 
                    j++;
                    if (j==8) cout << " wrong ufToA_stop value: " << ufToA_stop << endl; 
                } 

                if (top) last_hb = last_heartbeat_top; 
                else last_hb = last_heartbeat_bottom;  
                

                cout << "pixelhit: top " << (int) top << " last hb: "  << last_hb << " ToA "  << ToA << ' ' << endl;  

                // UInt_t ratio_VCO_CKDLL=16;



                // Double_t toa =  

                // return (uftoa_list.index(uftoa))

                

                h2->Fill(Col,Row,1);
                h2tot->Fill(Col,Row,ToT);
                h2toa->Fill(Col,Row,ToA);

                h1tot->Fill(ToT);
                h1toa->Fill(ToA);
                h1pileup->Fill(Pileup);

                h1ufToA_start->Fill(ufToA_start);   
                h1ufToA_stop->Fill(ufToA_stop);
                h1fToA_rise->Fill(fToA_rise);
                h1fToA_fall->Fill(fToA_fall);

                t2->Fill();
     
                npixelhits++;
            }  
            cout << g << '/' << infi_length << endl;    
            cout << npixelhits << ' ' << ntimestamps << ' ' << ntpx4markers << endl;
        }
    }
    cout << "npixelhits: " << npixelhits << " ntimestamps: " << ntimestamps << " ntpx4markers: " << ntpx4markers << endl;

    h1->Write();  
    h1hb->Write();
    h1eoc2->Write();  

    h1tot->Write();  
    h1toa->Write();  
    h1pileup->Write();  

    h1ufToA_start->Write();  
    h1ufToA_stop->Write();  
    
    h1fToA_rise->Write();  
    h1fToA_fall->Write();  
   
    h2->Write();
    h2tot->Write();
    h2toa->Write();

    t2->Write();
    hb->Write();

    // close .root file
    f->Close();

    return 0;
}


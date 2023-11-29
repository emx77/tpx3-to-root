
#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

int tpx4_to_root(string filename, unsigned long nrawpixelhits=0) {

    gROOT->Reset();

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

          
    ifstream infi(filename.c_str(), ifstream::binary);

    if (!infi) {
        cout << " file: " << filename << " not found." << endl;
        return 1; 
    }

    infi.seekg (0, infi.end);
    ULong64_t infi_length = infi.tellg();
    infi.seekg (0, infi.beg);

    const Int_t dpl=8; // length of data packet
    
    ULong64_t *data_packet = new ULong64_t; 
    ULong64_t hb_packet = 0; 
    // temperary char buffer 
    Char_t *buffer = new Char_t[dpl];

    ULong64_t i=0;   
    ULong64_t npixelhits = 0; 
    ULong64_t ntimestamps = 0;
    ULong64_t ntpx4markers = 0;
    ULong64_t g = infi.tellg();

    ULong64_t maxhits =2000000000;
    if (nrawpixelhits!=0) maxhits=nrawpixelhits;
    while (i<maxhits && g<infi_length) {
              
        infi.read((Char_t*)data_packet, dpl); 
        i++; 
        g = infi.tellg();
        UInt_t eoc = ((*data_packet)>>55) & 0xFF;
        UInt_t top=(*data_packet>>63)&0x1;
        //if (top) { continue; }

        h1->Fill(eoc);

        if (eoc==0xE0) {
            ntimestamps++;
            hb_packet = *data_packet;
            ULong64_t heartbeat = (*data_packet & 0xFFFFFFFFFFFF);
            //cout << top << ' ' << heartbeat << ' ' << heartbeat*25e-9 <<  endl;
            continue;  
        }
        
        buffer = (Char_t*) data_packet;
        TString s((Char_t*)buffer); 
        //cout << dec << g << ' ' << s << '-' << s.SubString(0,8) << endl;
         
        if (s.Contains("TPX4")) {
            ntpx4markers++;
            continue;
        } 

        //if (s.Contains("TPX4")) {
        //    cout << dec << g << ' ' << hex << *data_packet << " (TPX4 packet) " << endl;
        //}

        uint16_t ToA = (*data_packet>>30) & 0xFFFF;
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
          
        //cout << ToA << ' ' << inv << endl;
        //ToA = inv;
        UInt_t ToT = (*data_packet>>1) & 0x7FF;
        UInt_t Pileup = (*data_packet) & 0x1;    
   
        ULong64_t addr=(*data_packet>>46) & 0x3ffff;
        UInt_t Pix=addr&0x1F;
        UInt_t Sp=(addr>>5)&0xF;
        UInt_t EoC=(addr>>9)&0xFF;
        UInt_t Top=(addr>>17)&0x1;
        if (Top) {
            EoC = 223-EoC;
            Sp=15-Sp;
            Pix=31-Pix;
        }
        UInt_t Col=2*EoC + (int)(Pix%8/4);
        UInt_t Row=Pix%4 + (int)(Pix/8)*4 + Sp*16 + Top*256;


        if (eoc>223) { 
            cout << dec << Top << ' ' << g << ' ' << hex << *data_packet << dec << " eoc " << eoc << ' ' <<  endl;
            //cout << dec << g << ' ' << hex << *data_packet << dec << ' ' << eoc << ' '<< Col << ' ' << Row << ' ' << ToA << ' ' << ToT << ' ' << Pileup << endl;
        }
        //cout << "    " << Top << ' ' << ToA << endl;    
        
        UInt_t ufToA_start = (*data_packet>>26) & 0xF;
        UInt_t ufToA_stop  = (*data_packet>>22) & 0xF;
        UInt_t fToA_rise   = (*data_packet>>17) & 0x1F;
        UInt_t fToA_fall   = (*data_packet>>12) & 0x1F;

        // uftoa_list = [15, 14, 12, 8, 0, 1, 3, 7]
        // return (uftoa_list.index(uftoa))

        h1eoc2->Fill(eoc); 

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
     
        npixelhits++;
    }  
    cout << g << '/' << infi_length << endl;    
    cout << npixelhits << ' ' << ntimestamps << ' ' << ntpx4markers << endl;

    return 0;
}


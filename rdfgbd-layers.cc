#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

int main(int argc, char* argv[]) {

   string sdcd, srecpqr, infile, sbegframe, sendframe;
   if (argc == 6) { sdcd = argv[1]; srecpqr = argv[2]; infile = argv[3]; sbegframe = argv[4]; sendframe = argv[5];  }
      else { cout << ">> Please provide dcd trajectory, receptor.pqr, rdf input file, and first & last frame to be analyzed\n>> Usage: ./exec traj.dcd receptor.pqr rdf.in first_frame last_frame\n"; return 0; }

/// Read an input file with atom numbers of ref pts and prx atoms, etc. ///   
   ifstream inf;
   inf.open(infile.c_str(), ios::in);
   string infline;
   int i=0; 

  while (inf.is_open()) {
     getline (inf, infline);
     if (inf.eof()) { break; }
     i++;
   }
   int nline = i;
   inf.clear(); inf.seekg(0, ios::beg);

   int k = 0;
   string lines[nline]; 

   while (inf.is_open() && k < nline) {
    getline (inf, infline);
    lines[k] = infline;  
    k++;
   }

   int nrep = atoi(lines[1].c_str());
   int freq = atoi(lines[3].c_str());
   int nligatoms = atoi(lines[6].c_str());
   int ligrefatom = atoi(lines[8].c_str());
   int dcdfreq = atoi(lines[10].c_str());
   int ts = atoi(lines[12].c_str());   
   string indivrdfs = lines[14];
   double start = atof(lines[16].c_str());
   double width = atof(lines[18].c_str());
   cout << ">> Reading " << nrep << " ligand replicates in trajectory\n>> Reading every " << freq << " frames\n";

   inf.close();

   int begframe = atoi(sbegframe.c_str());
   int endframe = atoi(sendframe.c_str());
   
/// Get number of atoms in receptor ///
   ifstream recpqr; recpqr.open(srecpqr.c_str(), ios::in); 
   string recline;
   int nrecatoms = 0;
   while (recpqr.is_open()) {
     getline (recpqr, recline);
     if (recpqr.eof()) { break; }
     if (recline.find("ATOM") != std::string::npos) { nrecatoms++; }
   }
   recpqr.close();
   cout << ">> " << nrecatoms << " atoms in receptor\n";

/// Read and store stationary receptor coordinates from pqr ///
   double rec[3][nrecatoms];
   ifstream recpqr2; recpqr2.open(srecpqr.c_str(), ios::in);

   for (int i=0; i<nrecatoms; i++) {
    getline (recpqr2, recline);
    rec[0][i] = atof(recline.substr(33, 7).c_str());
    rec[1][i] = atof(recline.substr(43, 7).c_str());
    rec[2][i] = atof(recline.substr(53, 7).c_str());
    }
    recpqr2.close();

/// Reading DCD trajectory ///
   ifstream dcd;
   dcd.open(sdcd.c_str(), ios::in|ios::binary|ios::ate);
   long length; 
   length = dcd.tellg();

   cout << ">> Size of file is " << length/(1024*1024) << " MB.\n";
   char * buffer = new char[length];
   dcd.seekg(0, ios::beg);
   cout << "** Reading DCD file\n";
   long headlen = 0;
   int tmp, natoms;
   dcd.read((char*)&tmp, sizeof(int));
   headlen = headlen + 8 + tmp;
   dcd.seekg(headlen);
   dcd.read((char*)&tmp, sizeof(int));
   headlen = headlen + 8 + tmp;
   dcd.seekg(headlen);
   dcd.read((char*)&tmp, sizeof(int));
   headlen = headlen + 8 + tmp;
   dcd.read((char*)&natoms, sizeof(int));
   long framelen = natoms*12 + 80;
   long nframe = ((length-headlen)/framelen);

   cout << ">> " << natoms << " atoms in dcd file (is this right? It may not matter.\n>> " << nframe << " frames in dcd file\n";
  
   double cx, cy, cz; // Unit cell dimensions
   float coords[natoms][3][nframe];
   for (int i=0; i<endframe; /*i<nframe;*/ i++) {
     dcd.seekg(long(headlen) + framelen*i + 4); 
     dcd.read((char*)&cx, sizeof(double));
     dcd.seekg(long(headlen) + framelen*i + 20);
     dcd.read((char*)&cy, sizeof(double));
     dcd.seekg(long(headlen) + framelen*i + 44);
     dcd.read((char*)&cz, sizeof(double));

     dcd.seekg(long(headlen) + framelen*i + 56);
     dcd.read((char*)coords, sizeof(float)); // What is coords here, could it just as well be some other float var? Try

     for (int j=0;j<3;j++) {
       for (int k=0;k<natoms;k++) {
         dcd.read((char*)&coords[k][j][i], sizeof(float));
         }
      dcd.read((char*)&tmp, sizeof(float));
      dcd.read((char*)&tmp, sizeof(float));
      }
     }

   dcd.close();
   delete[] buffer;
   int bins[nrep][32]; // Create our bins
   for (int z=0; z<nrep; z++) { 
      for (int y=0; y<32; y++) { bins[z][y] = 0; } // Set all bin values to 0 first
    }
   int nmeasured = 0;
   double analyzed = 0;
   int assocframes[nrep];
   for (int i=0; i<nrep; i++) { assocframes[i] = 0; }

   if (endframe > nframe) { endframe = nframe; }   

   for (int j=begframe; j<endframe; j+=freq) {
     double totdist = 0;
     int nassoc = 0;
     //cout << "---------------------------\n** Reading frame " << j+1 << endl << "---------------------------\n";   
      //for (int k = 0; k<nrep; k++) {
      for (int k = 0; k<nrep; k++) {
       double mindist = 10000;
       nmeasured++;
        for (int m=0; m<nligatoms; m++) {  // Find min z dist between surface and any ligand atom
          
         double dz = coords[(k*nligatoms)+m][2][j] - 15.0;
         if (dz < mindist) { mindist = dz; }        
        }
      
     if      (mindist <= start)                                          { bins[k][0]++;  }
     else if (mindist > start && mindist <= start+width)                 { bins[k][1] ++; }
     else if (mindist > start+width && mindist <= start+(2*width))       { bins[k][2] ++; }
     else if (mindist > start+(2*width) && mindist <= start+(3*width))   { bins[k][3] ++; }
     else if (mindist > start+(3*width) && mindist <= start+(4*width))   { bins[k][4] ++; }
     else if (mindist > start+(4*width) && mindist <= start+(5*width))   { bins[k][5] ++; }
     else if (mindist > start+(5*width) && mindist <= start+(6*width))   { bins[k][6] ++; }
     else if (mindist > start+(6*width) && mindist <= start+(7*width))   { bins[k][7] ++; }
     else if (mindist > start+(7*width) && mindist <= start+(8*width))   { bins[k][8] ++; }
     else if (mindist > start+(8*width) && mindist <= start+(9*width))   { bins[k][9] ++; }
     else if (mindist > start+(9*width) && mindist <= start+(10*width))  { bins[k][10] ++;}
     else if (mindist > start+(10*width) && mindist <= start+(11*width)) { bins[k][11] ++; }
     else if (mindist > start+(11*width) && mindist <= start+(12*width)) { bins[k][12] ++; }
     else if (mindist > start+(12*width) && mindist <= start+(13*width)) { bins[k][13] ++; }
     else if (mindist > start+(13*width) && mindist <= start+(14*width)) { bins[k][14] ++; }
     else if (mindist > start+(14*width) && mindist <= start+(15*width)) { bins[k][15] ++; }
     else if (mindist > start+(15*width) && mindist <= start+(16*width)) { bins[k][16] ++; }
     else if (mindist > start+(16*width) && mindist <= start+(17*width)) { bins[k][17] ++; }
     else if (mindist > start+(17*width) && mindist <= start+(18*width)) { bins[k][18] ++; }
     else if (mindist > start+(18*width) && mindist <= start+(19*width)) { bins[k][19] ++; }
     else if (mindist > start+(19*width) && mindist <= start+(20*width)) { bins[k][20] ++; }
     else if (mindist > start+(20*width) && mindist <= start+(21*width)) { bins[k][21] ++; }
     else if (mindist > start+(21*width) && mindist <= start+(22*width)) { bins[k][22] ++; }
     else if (mindist > start+(22*width) && mindist <= start+(23*width)) { bins[k][23] ++; }
     else if (mindist > start+(23*width) && mindist <= start+(24*width)) { bins[k][24] ++; }
     else if (mindist > start+(24*width) && mindist <= start+(25*width)) { bins[k][25] ++; }
     else if (mindist > start+(25*width) && mindist <= start+(26*width)) { bins[k][26] ++; }
     else if (mindist > start+(26*width) && mindist <= start+(27*width)) { bins[k][27] ++; }
     else if (mindist > start+(27*width) && mindist <= start+(28*width)) { bins[k][28] ++; }
     else if (mindist > start+(28*width) && mindist <= start+(29*width)) { bins[k][29] ++; }
     else if (mindist > start+(29*width) && mindist <= start+(30*width)) { bins[k][30] ++; }
     else if (mindist > start+(30*width) && mindist <= start+(31*width)) { bins[k][31] ++; }
     //cout << "REP " << k+1 << " is " << mindist << "A away from DNA in frame " << j+1 << endl;
     }
   analyzed++;
   }

   if (indivrdfs == "YES") {
   for (int i=0; i<nrep; i++) { cout << "MOL " << i+1 << " rdf: " << endl;
       for (int j=0; j<32; j++) { cout << j*5 << " to " << (j+1)*5 << ": " << bins[i][j] << endl; }
     }
   }

   // Calculate 22 bins that hold sum of all MOLs in each bin
   int totbin[32];
   for (int t=0; t<32; t++) { totbin[t] = 0; } // Set all bins to 0 first
   for (int i=0; i<32; i++) { 
     for (int j=0; j<nrep; j++) { totbin[i] += bins[j][i]; }
    cout << ">> Bin " << i+1 << " raw total: " << totbin[i] << endl;
    }
   cout << "Averages for each bin:\n";
   for (int i=0; i<32; i++) {
    cout << totbin[i]/analyzed << endl;
   }
   cout << "Total MOL distances measured: " << nmeasured << endl;
 
 return 0;
}

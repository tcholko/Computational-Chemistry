#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

int main(int argc, char* argv[]) {

   string sdcd, infile, sbegframe, sendframe;
   if (argc == 5) { sdcd = argv[1]; infile = argv[2]; sbegframe = argv[3]; sendframe = argv[4];  }
      else { cout << ">> Please provide dcd trajectory, rdf input file, and first and last frame to be analyzed\n>> Usage: ./exec traj.dcd rdf.in first_frame last_frame\n"; return 0; }

///// Read an input file with atom numbers of ref pts and prx atoms, etc. /////   
   ifstream inf;
   inf.open(infile.c_str(), ios::in);
   string infline;
   int i=0; 
   double ref[10];// Number of ref pts is 10 for now

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

   int refbeg = atoi(lines[6].c_str());
   int refend = atoi(lines[8].c_str());
   int nligatoms = atoi(lines[10].c_str());
   int refnatom = (refend-refbeg)+1;
   cout << "refbeg= " << refbeg<<" refend= "<< refend<<"refnatom= "<<refnatom<<endl;
   int nmol = atoi(lines[1].c_str());
   int freq = atoi(lines[3].c_str());   
   cout << ">> " << nmol << " molecules in rdf\n>> Reading every " << freq << " frames\n";

   double mol[nmol];
   for (int j=0; j<nmol; j++) { 
     mol[j] = atof(lines[j+13].c_str());  // Reads in RDF molecule  atom numbers
     cout << "mol " <<j+1<< " atom is " << mol[j] << endl; 
    }

   inf.close();

///// Reading DCD trajectory ////// 
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

   cout << ">> " << natoms << " atoms in dcd file\n>> " << nframe << " frames in dcd file\n";
  
   double cx, cy, cz; // Unit cell dimensions
   float coords[natoms][3][nframe];
   for (int i=0;i<nframe;i++) {
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
   float bins[nmol][22]; // Create our bins
   for (int z=0; z<nmol; z++) { 
      for (int y=0; y<22; y++) { bins[z][y] = 0; } // Set all bin values to 0 first
    }
   int nmeasured = 0, analyzed = 0;
   int begframe = atoi(sbegframe.c_str());
   int endframe = atoi(sendframe.c_str());
   if (endframe > nframe) { endframe = nframe; }   

   for (int j=begframe; j<endframe; j+=freq) {
     double totdist = 0;
     cout << "---------------------------\n** Reading frame " << j+1 << endl << "---------------------------\n";
     for (int k = 0; k<nmol; k++) { 
       double mindist = 1000; 
       int molfirstatom = mol[k]-1;   
       nmeasured++;
      for (int m=0; m<nligatoms; m++) {
       for (int l=refbeg-1 ; l<refend; l++) {
         //int refpt = ref[l]-1;
         double dx2 = pow((coords[l][0][j]-coords[molfirstatom+m][0][j]), 2);
         double dy2 = pow((coords[l][1][j]-coords[molfirstatom+m][1][j]), 2);
         double dz2 = pow((coords[l][2][j]-coords[molfirstatom+m][2][j]), 2);
         double dist = sqrt(dx2 + dy2 + dz2);
         if (dist < mindist) { mindist = dist; }        
       }
      }
     if      (mindist <= 5)                   { bins[k][0] ++; }
     else if (mindist > 5  && mindist <= 10)  { bins[k][1] ++; }
     else if (mindist > 10 && mindist <= 15)  { bins[k][2] ++; }
     else if (mindist > 15 && mindist <= 20)  { bins[k][3] ++; }
     else if (mindist > 20 && mindist <= 25)  { bins[k][4] ++; }
     else if (mindist > 25 && mindist <= 30)  { bins[k][5] ++; }
     else if (mindist > 30 && mindist <= 35)  { bins[k][6] ++; }
     else if (mindist > 35 && mindist <= 40)  { bins[k][7] ++; }
     else if (mindist > 40 && mindist <= 45)  { bins[k][8] ++; }
     else if (mindist > 45 && mindist <= 50)  { bins[k][9] ++; }
     else if (mindist > 50 && mindist <= 55)  { bins[k][10] ++;}
     else if (mindist > 55 && mindist <= 60)  { bins[k][11] ++; }
     else if (mindist > 60 && mindist <= 65)  { bins[k][12] ++; }
     else if (mindist > 65 && mindist <= 70)  { bins[k][13] ++; }
     else if (mindist > 70 && mindist <= 75)  { bins[k][14] ++; }
     else if (mindist > 75 && mindist <= 80)  { bins[k][15] ++; }
     else if (mindist > 80 && mindist <= 85)  { bins[k][16] ++; }
     else if (mindist > 85 && mindist <= 90)  { bins[k][17] ++; }
     else if (mindist > 90 && mindist <= 95)  { bins[k][18] ++; }
     else if (mindist > 95 && mindist <= 100) { bins[k][19] ++; }
     else if (mindist > 100 && mindist <=105) { bins[k][20] ++; }
     else if (mindist > 105 && mindist <=110) { bins[k][21] ++; }
     totdist += mindist;
     //cout << "MOL " << k+1 << " is " << mindist << "A away from DNA in frame " << j+1 << endl;
     }
     analyzed++;
     double avgdist = totdist/nmol;
     //cout << ">> Average MOL distance from DNA in frame " <<j+1<< " is " << avgdist << " A\n";
    }

   for (int i=0; i<nmol; i++) { cout << "MOL " << i+1 << " rdf: " << endl;
       for (int j=0; j<22; j++) { cout << j*5 << " to " << (j+1)*5 << ": " << bins[i][j] << endl; }
   }
   // Calculate 22 bins that hold sum of all MOLs in each bin
   float totbin[22];
   float avgbin;
   for (int t=0; t<22; t++) { totbin[t] = 0; } // Set all bins to 0 first
   for (int i=0; i<22; i++) { 
     for (int j=0; j<nmol; j++) { totbin[i] += bins[j][i]; avgbin = bins[j][i]/analyzed;}
    cout << ">> Bin " << i+1 << " Total: " << totbin[i] << " | Average: " << avgbin << endl;
   }
   cout << "Total MOL distances measured: " << nmeasured << endl << "Frames analyzed: " << analyzed << endl;
   
return 0;
}

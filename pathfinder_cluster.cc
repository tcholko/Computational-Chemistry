#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <bitset>
#include <sstream>
using namespace std;

int main(int argc, char* argv[]) {

   string sdcd, scomatom, snrep, sfreq, infile;
   if (argc == 6) { sdcd = argv[1], scomatom = argv[2], snrep = argv[3], sfreq = argv[4], infile = argv[5]; }
   else { cout << "**Command line error\n**Usage: ./exec traj.dcd COM_atom_# #_replicates frame_read_frequency rec_input_coords\n"; return 0;}

  float rec_coords[5][2];

/// Read input file ///
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

   rec_coords[0][0] = atof(lines[1].c_str());
   rec_coords[0][1] = atof(lines[2].c_str());
   rec_coords[1][0] = atof(lines[3].c_str());
   rec_coords[1][1] = atof(lines[4].c_str()); 
   rec_coords[2][0] = atof(lines[5].c_str());
   rec_coords[2][1] = atof(lines[6].c_str());
   rec_coords[3][0] = atof(lines[7].c_str());
   rec_coords[3][1] = atof(lines[8].c_str());
   rec_coords[4][0] = atof(lines[9].c_str());
   rec_coords[4][1] = atof(lines[10].c_str());

   inf.close();
  
/// Reading DCD trajectory ///
  int freq = atoi(sfreq.c_str());
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
   float coords[natoms][3][nframe/freq];
   for (int i=0; i<nframe; i++) {
    dcd.seekg(long(headlen) + framelen*i*freq + 4);
    dcd.read((char*)&cx, sizeof(double));
    dcd.seekg(long(headlen) + framelen*i*freq + 20);
    dcd.read((char*)&cy, sizeof(double));
    dcd.seekg(long(headlen) + framelen*i*freq + 44);
    dcd.read((char*)&cz, sizeof(double));
    dcd.seekg(long(headlen) + framelen*i*freq + 56);
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

   int bind_rec_cnt[5]; for (int i=0; i<5; i++) { bind_rec_cnt[i] = 0; }
   int COMatom = atoi(scomatom.c_str()), nrep = atoi(snrep.c_str()), currentrep = 1, twoD = 0, threeD = 0;
   int ligatoms = natoms/nrep;
   for (int n=COMatom; n<=COMatom+((nrep-1)*ligatoms); n+=ligatoms) { // Iterates over all replicates
     cout << "** Finding binders for replicate " << currentrep <<"..."<< endl;
     for (int m=0; m<(nframe/freq)-1; m++) {
       float dz = coords[n-1][2][m+1]-coords[n-1][2][m]; 

      if (dz > 200 /*&& coords[n-1][m] < 150*/) {
      int adsorbed = 0;
       for (int f=0; f<20; f++) {  // Checking backwards f frames for adsorption
         for (int l=0; l<ligatoms; l++) {
          float disttoSAM = coords[((currentrep-1)*ligatoms)+l][2][m-f] - 15.46;
          if (disttoSAM <= 10) { adsorbed++; break; } 
         }                                            
        }

        if (adsorbed >= 1) { twoD++; /* Here is where code to count lifetime of traj should go. Or maybe just build in a timer for all trajs to begin with? */ } 
        else { threeD++; } // formerly if (coords[n-1][m-20] < 70) 

        cout <<"dz = "<<dz<<" at frame "<<(m*freq)+1<< " comatom= " <<n<< "\n2D binders: "<<twoD<<"  3D binders: "<<threeD<<endl<<"Total: "<<twoD+threeD<<endl;   

        float mindist = 100000, avg1 = 0, avg2 = 0, avg3 = 0, avg4 = 0, avg5 = 0;
        int bind_rec;

        for (int r=0; r<ligatoms; r+=360) {
     // May have to avg over many atoms to get avg dist to each probe and the smallest avg is the one it bound too
        //float mindist = 100000, avg1 = 0, avg2 = 0, avg3 = 0, avg4 = 0, avg5 = 0;
        //int bind_rec;
        float drec1_x2 = pow((coords[r][0][m] - rec_coords[0][0]), 2);
        float drec1_y2 = pow((coords[r][1][m] - rec_coords[0][1]), 2);
        float drec1    = sqrt(drec1_x2 + drec1_y2);
        float drec2_x2 = pow((coords[r][0][m] - rec_coords[1][0]), 2);
        float drec2_y2 = pow((coords[r][1][m] - rec_coords[1][1]), 2);
        float drec2    = sqrt(drec2_x2 + drec2_y2);
        float drec3_x2 = pow((coords[r][0][m] - rec_coords[2][0]), 2);
        float drec3_y2 = pow((coords[r][1][m] - rec_coords[2][1]), 2);
        float drec3    = sqrt(drec3_x2 + drec3_y2);
        float drec4_x2 = pow((coords[r][0][m] - rec_coords[3][0]), 2);
        float drec4_y2 = pow((coords[r][1][m] - rec_coords[3][1]), 2);
        float drec4    = sqrt(drec4_x2 + drec4_y2);
        float drec5_x2 = pow((coords[r][0][m] - rec_coords[4][0]), 2);
        float drec5_y2 = pow((coords[r][1][m] - rec_coords[4][1]), 2);
        float drec5    = sqrt(drec5_x2 + drec5_y2);
        avg1 += drec1; avg2 += drec2; avg3 += drec3; avg4 += drec4; avg5 += drec5;
       }
        avg1 /= 5; avg2 /= 5; avg3 /= 5; avg4 /= 5; avg5 /= 5;        
        if (avg1 < mindist) { mindist = avg1; bind_rec = 0; }
        if (avg2 < mindist) { mindist = avg2; bind_rec = 1; }
        if (avg3 < mindist) { mindist = avg3; bind_rec = 2; }
        if (avg4 < mindist) { mindist = avg4; bind_rec = 3; }
        if (avg5 < mindist) { mindist = avg5; bind_rec = 4; }


/*        if (drec1 < mindist) { mindist = drec1; bind_rec = 0; }
        if (drec2 < mindist) { mindist = drec2; bind_rec = 1; }
        if (drec3 < mindist) { mindist = drec3; bind_rec = 2; }
        if (drec4 < mindist) { mindist = drec4; bind_rec = 3; }
        if (drec5 < mindist) { mindist = drec5; bind_rec = 4; }
*/       
        cout << "Binding occured at receptor" << bind_rec << endl;
        bind_rec_cnt[bind_rec]++;
      }
    }
     
     currentrep++;
    }
   cout << "Binding total at each receptor: \n 1: " << bind_rec_cnt[0] << endl<< "2: " << bind_rec_cnt[1] << endl << "3: " << bind_rec_cnt[2] << endl << "4: " << bind_rec_cnt[3] << endl << "5: " << bind_rec_cnt[4] << endl;
return 0;
}

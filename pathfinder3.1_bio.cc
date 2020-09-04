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

   string sdcd, scomatom, snrep, sfreq;
   if (argc == 5) { sdcd = argv[1], scomatom = argv[2], snrep = argv[3], sfreq = argv[4]; }
   else { cout << "**Command line error\n**Usage: ./exec traj.dcd COM_atom_# #_replicates frame_read_frequency\n"; return 0;}
   
   ifstream dcd;
   dcd.open(sdcd.c_str(), ios::in|ios::binary|ios::ate);
   long length; 
   length = dcd.tellg();

   cout << "**Size of file is " << length/(1024*1024) << " MB.\n";
   char * buffer = new char[length];
   dcd.seekg(0, ios::beg);
   long headlen = 0;
   int tmp, natoms, freq = atoi(sfreq.c_str());
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
   cout << natoms << " atoms in dcd file\n" << nframe << " frames in dcd file\n";
  
   double cx, cy, cz; // Unit cell dimensions
   float coords[natoms][nframe/freq];  //**********Change back to just nframe if this doesnt work*********!!!!!!!!!!
   // Coordinates stored as: frame 1 -> all x, all y, all z (ordered by atom index #), frame 2 -> all x, all y, all z, frame 3... etc
   cout << "** Reading ligand coordinates...\n";
   for (int i=0;i<nframe/freq;i++) {
     dcd.seekg(long(headlen) + framelen*i*freq + 4); 
     dcd.read((char*)&cx, sizeof(double));
     dcd.seekg(long(headlen) + framelen*i*freq + 20);
     dcd.read((char*)&cy, sizeof(double));
     dcd.seekg(long(headlen) + framelen*i*freq + 44);
     dcd.read((char*)&cz, sizeof(double));
     dcd.seekg(long(headlen) + framelen*i*freq + 56);
     dcd.read((char*)coords, sizeof(float)); // What is coords here, could it just as well be some other float var? Try
     // Instead of reading for all 3 coords, after reading up to here, seek up to beginning of z coords by reading 2*natoms floats
     for (int j=0;j<2*natoms;j++) { 
       dcd.read((char*)&tmp, sizeof(float)); 
       }
     dcd.read((char*)&tmp, sizeof(float)); // Then read past headers that come after each dimension of coords, 2 floats per dimension
     dcd.read((char*)&tmp, sizeof(float));
     dcd.read((char*)&tmp, sizeof(float));
     dcd.read((char*)&tmp, sizeof(float));
     for (int k=0;k<natoms;k++) {
       dcd.read((char*)&coords[k][i], sizeof(float));
       }
     }
   dcd.close();
   delete[] buffer;
   
   int COMatom = atoi(scomatom.c_str()), nrep = atoi(snrep.c_str()), currentrep = 1, twoD = 0, threeD = 0; // Make COMatom and nrep input vars
   int ligatoms = natoms/nrep;
   for (int n=COMatom; n<=COMatom+((nrep-1)*ligatoms); n+=ligatoms) { // Iterates over all replicates
     cout << "** Finding binders for replicate " << currentrep <<"..."<< endl;
     for (int m=0; m<(nframe/freq)-1; m++) {
       float dz = coords[n-1][m+1]-coords[n-1][m]; 

      if (dz > 200 /*&& coords[n-1][m] < 150*/) {
      int adsorbed = 0;
       for (int f=0; f<20; f++) {  // Checking backwards f frames for adsorption
         for (int l=0; l<ligatoms; l++) {
          float disttoSAM = coords[((currentrep-1)*ligatoms)+l][m-f] - 15.46;
          if (disttoSAM <= 10) { adsorbed++; break; } 
         }                                            
        }

        if (adsorbed >= 1) { twoD++; /* Here is where code to count lifetime of traj should go. Or maybe just build in a timer for all trajs to begin with? */ } 
        else { threeD++; } // formerly if (coords[n-1][m-20] < 70) 

        cout <<"dz = "<<dz<<" at frame "<<(m*freq)+1<< " comatom= " <<n<< "\n2D binders: "<<twoD<<"  3D binders: "<<threeD<<endl<<"Total: "<<twoD+threeD<<endl;   
        }
      }
     
     currentrep++;
    }

return 0;
}

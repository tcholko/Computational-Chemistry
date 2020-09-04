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

   string sdcd;
   if (argc == 2) { sdcd = argv[1]; }
   
   ifstream dcd;
   dcd.open(sdcd.c_str(), ios::in|ios::binary|ios::ate); //ate mean open w/ cursor at end of file
   long length; 
   length = dcd.tellg();

   cout << "Size of file is " << length/(1024*1024) << " MB.\n";
   char * buffer = new char[length];
   dcd.seekg(0, ios::beg);
   cout << "Reading DCD file\n";
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

   cout << natoms << " atoms in dcd file\n" << nframe << " frames in dcd file\n One frame is " << framelen << " bits\nHeadlen is " <<headlen<<" bits\n";
  
   double cx, cy, cz; // Unit cell dimensions
   float coords[natoms][3][nframe];
   // Coordinates stored as: frame 1 -> all x, all y, all z (ordered by atom index #), frame 2 -> all x, all y, all z, frame 3... etc
   for (int i=0;i<nframe;i++) {
     dcd.seekg(long(headlen) + framelen*i + 4);
     //cout << dcd.tellg() << endl; 
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
         cout << coords[k][j][i] << endl;
         }
       dcd.read((char*)&tmp, sizeof(float));
       dcd.read((char*)&tmp, sizeof(float));
       }
     }

   dcd.close();
   delete[] buffer;

return 0;
}

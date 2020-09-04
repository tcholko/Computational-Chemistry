#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
using namespace std;

int main(int argc, char* argv[]) {
   string sdcd;
   string savef;
   string comatom;
   string sndim;
   if (argc == 5) { sdcd = argv[1]; savef = argv[2]; comatom = argv[3]; sndim = argv[4]; } 
    else { cout << "Please specify your trajectory file, frame save frequency, C.O.M atom #, and diffusion dimensionality \nUsage: ./diffusion trajectory.dcd save_freq COM_atom_# ndim\n"; 
   return 0; }
   cout << "Warning: receptor center coordinates assumed to be 0, 0, 0. Change coordinates in code if this is not correct.\n";
   cout << "User Input: Trajectory: " << sinput <<", Trajectory save freq: "<<savef<<", Ligand COM atom #: "<<comatom<<", Diffusion process dimensionality: "<< sndim << endl; 
   double save = atof(savef.c_str());
   double datom = atof(comatom.c_str());
   double ndim = atof(sndim.c_str());

////////// Read DCD ///////////////////
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
   // Read past headers & other info
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
   // Reading unit cell info
     for (int i=0;i<nframe;i++) {
       dcd.seekg(long(headlen) + framelen*i + 4); 
       dcd.read((char*)&cx, sizeof(double));
       dcd.seekg(long(headlen) + framelen*i + 20);
       dcd.read((char*)&cy, sizeof(double));
       dcd.seekg(long(headlen) + framelen*i + 44);
       dcd.read((char*)&cz, sizeof(double));
       dcd.seekg(long(headlen) + framelen*i + 56);
       dcd.read((char*)coords, sizeof(float)); // What is coords here, could it just as well be some other float var? Try
   // Reading all coordinates into 3D array. k = natom; j = dimension; i = nframe
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
/////// Now calculate timsept and displacement between each frame for COM atom /////////////

START HERE NEXT TIM




/////// Determine traj length ////////////////////////
   ifstream nextinput;
   nextinput.open(sinput.c_str(), ios::in);
   string strcoords;
   int l = 0;
   double atom, xt, yt, zt;
   double t = 0;   

   while (nextinput.is_open() && l < nframe) {
     getline (nextinput, strcoords);
     if (strcoords.find("END") != std::string::npos) { continue; } //Do nothing
        else { atom = atof(strcoords.substr(5, 6).c_str()); }
  
      if (atom == datom) {
       xt = atof(strcoords.substr(30, 8).c_str());
       yt = atof(strcoords.substr(38, 8).c_str());
       zt = atof(strcoords.substr(46, 8).c_str());
       l++;
       
       double comx = 0.; double comy = 0.; double comz = 0.;
       double distfromCOM2 = pow((xt-comx), 2) + pow((yt-comy), 2) + pow((zt-comz), 2);
       double distfromCOM = sqrt(distfromCOM2);
       double dt = (distfromCOM-75)*0.00212; //ts range divided by scale distance range 
        if (dt > 1.0) { dt = 1.0; }
        else if (dt < 0.1) { dt = 0.1; }
       cout << "xt, yt, zt, dt: " << xt <<" "<<yt<<" "<<zt<<" "<<dt<<endl; 
       t += dt; 
      }
    }
  t*= save;
  nextinput.close();
/////////////////////////////////////////////////////////////////////
/*
   input.open(sinput.c_str(), ios::in);
   string strfr;
   int count = 0;
   cout << datom;
   while (getline (input, strfr)) {
      if (strfr.find("END") != std::string::npos) { count++;
      }
     }
   int nframe = count;
 double trajlength = nframe*ts;   
*/
   ifstream traj;
   traj.open(sinput.c_str(), ios::in);
   string strline[nframe];
   double x[nframe]; double y[nframe]; double z[nframe];
   //double atom;
   int i = 0;

   while (traj.is_open() && i < nframe) {
     getline (traj, strline[i]);
     if (strline[i].find("END") != std::string::npos) { continue; } //Do nothing
      //if (strline[i].find(" 125 ") != std::string::npos) {
      else { atom = atof(strline[i].substr(5, 6).c_str()); }
       //cout << atom; 
      
      //if (strline[i].substr(8, 3).c_str() == comatom) {
      if (atom == datom) {
       double x1 = atof(strline[i].substr(30, 8).c_str());
       double y1 = atof(strline[i].substr(38, 8).c_str());
       double z1 = atof(strline[i].substr(46, 8).c_str());
      
       x[i] = x1; y[i] = y1; z[i] = z1;  
       i++; 
      }
     }
   traj.close();
   
   double tdisp = 0;
   for (int j = 0; j < nframe-1; j++) {
      double dx = x[j+1] - x[j];
        if (dx > 200) {
          dx -= 316; }  // Adjusts by 2*exbound
        else if (dx < -200) {
          dx += 316; }
      double dx2 = dx*dx;
      double dy = y[j+1] - y[j];
        if (dy > 200) {
          dy -= 316; }  
        else if (dy < -200) {
          dy += 316; }
      double dy2 = dy*dy;
      double dz = z[j+1] - z[j];
        if (dz > 200) {
          dz -= 350; }  // Adjust for binding reset to b-plane
        else if (dz < -200) {
          dz += 316; }
      double dz2 = dz*dz;
      double disp = sqrt(dx2 + dy2 + dz2);
      
      tdisp += disp;
      cout << tdisp << "\n"; 
      }
  
   cout <<"Trajectory length: "<< t <<" ps"<<"\n"<<"D = "<< (tdisp/(t*2*ndim))*1e12*1e-16 <<" cm^2/s"<<"\n";
  
   return 0;
}
  

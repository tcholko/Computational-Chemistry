#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
using namespace std;

int main(int argc, char* argv[]) {
   ifstream input;
   string sinput;
   string savef;
   string comatom;
   string sndim;
   if (argc == 5) { sinput = argv[1]; savef = argv[2]; comatom = argv[3]; sndim = argv[4]; } 
    else { cout << "Please specify your trajectory file, frame save frequency, C.O.M atom #, and diffusion dimensionality \nUsage: ./diffusion trajectory.pdb save_freq COM_atom_# ndim\n"; 
   return 0; }
   cout << "Warning: receptor center coordinates assumed to be 0, 0, 0. Change coordinates in code if this is not correct.\n";
   cout << "User Input: Trajectory: " << sinput <<", Trajectory save freq: "<<savef<<", Ligand COM atom #: "<<comatom<<", Diffusion process dimensionality: "<< sndim << endl; 
   double save = atof(savef.c_str());
   double datom = atof(comatom.c_str());
   double ndim = atof(sndim.c_str());

 input.open(sinput.c_str(), ios::in);
   string strfr;
   int count = 0;
   cout << "COM atom = " << datom;
   while (getline (input, strfr)) {
      if (strfr.find("END") != std::string::npos) { count++;
      }
     }
   int nframe = count;
   cout << "# frames = " << nframe;
 input.close();

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
       
       double comx = -60.; double comy = -60.; double comz = 0.;
       double distfromCOM2 = pow((xt-comx), 2) + pow((yt-comy), 2) + pow((zt-comz), 2);
       double distfromCOM = sqrt(distfromCOM2);
       double dt = 0.1 + ((distfromCOM-75)*0.00212); //0.00212 = ts range divided by scale distance range 
        if (dt > 1.0) { dt = 1.0; }
        else if (dt < 0.1) { dt = 0.1; }
       cout << "xt, yt, zt, dt: " << xt <<" "<<yt<<" "<<zt<<" "<<dt<<endl; 
       t += dt; 
      }
    }
  t*= save;
  nextinput.close();

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
  

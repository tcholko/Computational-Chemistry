#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
using namespace std;

int main(int argc, char* argv[]) {
   ifstream input;
   string sinput;
   string timestep;
   string comatom;
   string sndim;
   if (argc == 5) { sinput = argv[1]; timestep = argv[2]; comatom = argv[3]; sndim = argv[4]; } 
    else { cout << "Please specify your trajectory file, traj timestep, C.O.M atom #, and diffusion dimensionality \nUsage: ./diffusion trajectory.pdb timestep COM_atom_# ndim\n"; 
   return 0; }

   double ts = atof(timestep.c_str());
   double datom = atof(comatom.c_str());
   double ndim = atof(sndim.c_str());

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
 input.close();   

   ifstream traj;
   traj.open(sinput.c_str(), ios::in);
   string strline[nframe];
   double x[nframe]; double y[nframe]; double z[nframe];
   double atom;
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
          dx -= 400; }  // Adjusts by 2*exbound
        else if (dx < -200) {
          dx += 400; }
      double dx2 = dx*dx;
      double dy = y[j+1] - y[j];
        if (dy > 200) {
          dy -= 400; }  
        else if (dy < -200) {
          dy += 400; }
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
  
   cout <<"Trajectory length: "<< trajlength <<" ps"<<"\n"<<"D = "<< (tdisp/(trajlength*2*ndim))*1e12*1e-16 <<" cm^2/s"<<"\n";
  
   return 0;
}
  

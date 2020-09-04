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

   string sdcd, scomatom, snrep, sfreq, ssave;
   if (argc == 6) { sdcd = argv[1], ssave = argv[2], scomatom = argv[3], snrep = argv[4], sfreq = argv[5]; }
   else { cout << "**Command line error\n**Usage: ./exec traj.dcd traj_save_freq COM_atom_# #_replicates frame_read_frequency\n"; return 0;}
   
   ifstream dcd;
   dcd.open(sdcd.c_str(), ios::in|ios::binary|ios::ate);
   long length; 
   length = dcd.tellg();

   cout << "**Size of file is " << length/(1024*1024) << " MB.\n";
   char * buffer = new char[length];
   dcd.seekg(0, ios::beg);
   long headlen = 0;
   int tmp, natoms, ifreq = atoi(sfreq.c_str());
   double freq = atof(sfreq.c_str());
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
   float coords[natoms][3][nframe/ifreq];  //**********Change back to just nframe if this doesnt work*********!!!!!!!!!!
   // Coordinates stored as: frame 1 -> all x, all y, all z (ordered by atom index #), frame 2 -> all x, all y, all z, frame 3... etc
   cout << "** Reading ligand coordinates...\n";
   for (int i=0;i<nframe/ifreq;i++) {
     dcd.seekg(long(headlen) + framelen*i*freq + 4); 
     dcd.read((char*)&cx, sizeof(double));
     dcd.seekg(long(headlen) + framelen*i*freq + 20);
     dcd.read((char*)&cy, sizeof(double));
     dcd.seekg(long(headlen) + framelen*i*freq + 44);
     dcd.read((char*)&cz, sizeof(double));
     dcd.seekg(long(headlen) + framelen*i*freq + 56);
     dcd.read((char*)coords, sizeof(float)); // What is coords here, could it just as well be some other float var? Try
     // Instead of reading for all 3 coords, after reading up to here, seek up to beginning of z coords by reading 2*natoms floats
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

   double save = atof(ssave.c_str()), numhitsurface = 0.;
   float timer, surfhit_timer, totaltime = 0.0, twoDtotaltime = 0.0, threeDtotaltime = 0.0, totaltime2surf = 0.0;
   int COMatom = atoi(scomatom.c_str()), nrep = atoi(snrep.c_str()), currentrep = 1, twoD = 0, threeD = 0;
   int ligatoms = natoms/nrep, hitsurface, left, leftsurface = 0, trajcount = 1;
   for (int n=COMatom; n<=COMatom+((nrep-1)*ligatoms); n+=ligatoms) { // Iterates over all replicates
     cout << "** Finding binders for replicate " << currentrep <<"..."<< endl;
     timer = 0.0; hitsurface = 0; left = 0;
     for (int m=0; m<(nframe/ifreq)-1; m++) {
       float dz = coords[n-1][2][m+1]-coords[n-1][2][m]; 
       float dcenterx = coords[n-1][0][m] - 150.2;  // !! Change rec COM x, y, and z if not 0.0
       float dcentery = coords[n-1][1][m] - 67.5;
       float dcenterz = coords[n-1][2][m] - 8.3;
       float dcenter2 = pow(dcenterx, 2) + pow(dcentery, 2) + pow(dcenterz, 2);
       float ts = (0.1 + ((dcenter2-(75*75)) *((1.0 - 0.1)/((500*500)-(75*75))) ));  // ts range / scale dist range
       if (ts > 1.0) {ts = 1.0;}
       if (ts < 0.1) {ts = 0.1;}
      //cout << "Timestep: " << ts << endl;
      ts*=(save*freq);
      timer += ts;

     if (hitsurface ==1 && left ==0) {        
         float check_leftsurface = coords[n-1][2][m] - 15.46;
          if (check_leftsurface > 75) { leftsurface++; left = 1; } // If ligand COM is 75 or more above surface, count as leaving surface, then mark as "left" so we don't overcount
        }

     if (hitsurface == 0) {
      surfhit_timer += ts;
      for (int p=0; p<ligatoms; p+=5) {  // Check every 5th atom to make it quick
         float dz_to_SAM = coords[((currentrep-1)*ligatoms)+p][2][m] - 15.46; 
         if (dz_to_SAM < 10) { cout << "Trajectory " << trajcount << " time to first surface contact: " << timer/1E6 << " us\n"; 
            hitsurface = 1; totaltime2surf += surfhit_timer; surfhit_timer = 0.0; numhitsurface++; break; } // Break loop once we find any contact 
         }
      }
      if (dz > 200) {
      int adsorbed = 0;
       for (int f=0; f<30; f++) {  // Checking backwards f frames for steady adsorption
         for (int l=0; l<ligatoms; l++) {
          float disttoSAM = coords[((currentrep-1)*ligatoms)+l][2][m-f] - 15.46;
          if (disttoSAM <= 10) { adsorbed++; break; }  // Increment adsorbed each we find an adsorbed frame
         }                                            
        }

        if (adsorbed >= 10) { twoD++; twoDtotaltime += timer; } 
        else { threeD++; threeDtotaltime += timer; }  
       cout <<"Binding event detected. dz = "<<dz<<" at frame "<<(m*ifreq)+1<< " comatom= " <<n<< "\nTraj lifetime = " << timer/1E6 << " us\n";
       cout <<"2D binders: "<<twoD<<"  3D binders: "<<threeD<<endl<<"Total: "<<twoD+threeD<<endl;
       totaltime += timer;
       timer = 0.0;
       surfhit_timer = 0.0;
       hitsurface = 0;
       left = 0;
       trajcount++;   
        }
      }
     
     currentrep++;
    }
 
  cout << "Average trajectory binding time: " << totaltime/(twoD+threeD) << endl << "Average 2D binder time: " << twoDtotaltime/twoD << endl;
  cout << "Average 3D binder time: " << threeDtotaltime/threeD << endl << "Average first surface conctact " << totaltime2surf/(twoD+threeD) << endl;
  cout << ( numhitsurface /(twoD+threeD))*100. << "%" << " of the trajectories touched surface before binding\n" << leftsurface << " trajectories left the surface after touching\n";
return 0;
}

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
using namespace std;

int main(int argc, char* argv[]) {

   string sdcd, scomatom, snrep, sfreq;
   if (argc == 5) { sdcd = argv[1], scomatom = argv[2], snrep = argv[3], sfreq = argv[4]; }
   else { cout << "**Command line error\n**Usage: ./exec traj.dcd COM_atom_# n_replicates frame_read_frequency\n"; return 0;}
   cout << " * WARNING: Assuming writetraj = 5000 by default; change code if necessary\n";
   ifstream dcd;
   dcd.open(sdcd.c_str(), ios::in|ios::binary|ios::ate); //ate mean open w/ cursor at end of file
   long length; 
   length = dcd.tellg();

   cout << "> Size of file is " << length/(1024*1024) << " MB.\n";
   char * buffer = new char[length];
   dcd.seekg(0, ios::beg);
   long headlen = 0;
   int tmp, natoms; 
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
   cout << "> " << natoms << " atoms in dcd file\n" << "> " << nframe << " frames in dcd file\n";
  
   double cx, cy, cz; // Unit cell dimensions
   float coords[natoms][3][nframe];
   cout << "** Reading trajectory...\n";
   for (int i=0;i<nframe/freq;i++) {
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

/// Begin calculation of surface residence times ///
   cout << ">> Beginning calcultion of surface residence times\n";
   double eventrestime[5000];
   for (int i=0; i<5000; i++) { eventrestime[i] = 0; }
   double ttot = 0, tdwellevent = 0.0, dt = 0, reccenterx = 150., reccentery = 60., des_timer = 0.; // !! Assuming 0,0,0 receptor coords
   int nrep = atoi(snrep.c_str()), COMatom = atoi(scomatom.c_str()), ligatoms = natoms/nrep, nassoc = 0, ifreq = freq, adsorbed, e = 0;
   for (int n=0; n<nrep; n++) { 
     for (int f=0; f<nframe/freq; f++) {
      adsorbed = 0;
       for (int a=0; a<ligatoms; a++) {
          double surfdist = coords[a+n*ligatoms][2][f*ifreq] - 15.4; // !! Use correct surface z coordinate
          if (surfdist <= 5.0) {
            if (des_timer >= 4.0) {  // 4.0 requies >= 20ns desorb time to count new assoc. event if freq = 10
               //nassoc++; 
               if (nassoc >= 1) { 
                 eventrestime[nassoc-1] = tdwellevent;
                 tdwellevent = 0;
                 cout << "last event res time " << eventrestime[nassoc-1] << endl;
                }             
              nassoc++;
              } // END des_timer loop
            adsorbed = 1; des_timer = 0;
            double distx2 = pow((coords[(COMatom-1)+(n*ligatoms)][0][f] - reccenterx), 2);
            double disty2 = pow((coords[(COMatom-1)+(n*ligatoms)][1][f] - reccentery), 2);
            double centerdist2 = (distx2 + disty2);
            tdwellevent += dt;
            dt = 0.1 + ( (centerdist2 - (75*75) ) *( (1.0-0.1) /( (500*500) - (75*75) ) ) );
              if (dt < 0.1) { dt = 0.1; } if (dt > 1.0) { dt = 1.0; }
            dt *= (freq*5000); // !! 5000 is default 'writetraj' keyword
            //cout << "t_residence = " << tdwell << ", dt = " << dt <<" and adsorbed = "<<adsorbed<<" , nassoc = "<<nassoc<<" desb time= "<<des_timer<< endl;
            a = ligatoms;
           } // END surfdist loop   
       } // END atom loop
      if (adsorbed == 0) { 
         //cout << "t_residence replicate " << n+1 << " (desorbed) = " << tdwellevent << " ps\n"; 
         des_timer+=(freq*0.1); 
        }
      } // END frame loop
   } // END replicate loop

  for (int i=0; i<nassoc; i++) {
      ttot+=eventrestime[i];
    }

  cout << ">> Total surface res time: " << ttot/1e6 << " us, Num unique adsorptions: " << nassoc << ", Overall avg res time: " << ttot/(nassoc*1e6) << " us\n";
  double avg = ttot/nassoc; 
  double sumd2, sd;
  for (int i=0; i<nassoc; i++) { 
     sumd2 += pow(avg - eventrestime[i], 2);
     //cout << "event " << i+1 << " res time: " << eventrestime[i]/1e3 <<" ns\n";
    }
   sumd2 /= nassoc;
   sd = sqrt(sumd2);   
  //cout << "SD residence time: " << sd/1e3 << " ns\n";

/// Record and average only over long adosorptions, not just touches ///
/*  int nlong = 0;
  float longevents[nassoc];

  for (int i=0; i<nassoc; i++) { 
    if (eventrestime[i] > 1e3) { longevents[nlong] = eventrestime[i]; nlong++; } // if adsorbed longer than 1 ns, record it
   }

  double avg_2 = ttot/nassoc;
  double sumd2_2, sd_2;
  for (int i=0; i<nlong; i++) {
     sumd2_2 += pow(avg_2 - longevents[i], 2);
     cout << "Long event " << i+1 << " res time: " << longevents[i]/1e6 <<" us\n";
    }
  sumd2_2 /= nlong;
  sd_2 = sqrt(sumd2_2);
  cout << nlong << " long adsorption events found\n";
  cout << "SD residence time: " << sd_2/1e6 << " us\n";
*/

/// Throw out outliers
  int normal = 0, outlier = 0;
  float normevents[nassoc], outliers[nassoc];
  for (int i=0; i<nassoc; i++) {
    if ( (eventrestime[i]) < avg+(4*sd) ) { normevents[normal] = eventrestime[i]; normal++; } // if w/i 4 stdev of avg, record event
    else { eventrestime[i] = outliers[outlier]; outlier++; }
   }

  double norm_avg = 0;
   for (int i=0; i<normal; i++) { norm_avg += normevents[i]; }
  norm_avg /= normal;
  double sumd2_2, sd_2;
  for (int i=0; i<normal; i++) {
     sumd2_2 += pow(norm_avg - normevents[i], 2);
     cout << "Event " << i+1 << " res time: " << normevents[i]/1e3 <<" ns\n";
    }
  cout << "Average residence time minus outliers: " << norm_avg/1e3 << " ns\n";
  sumd2_2 /= normal;
  sd_2 = sqrt(sumd2_2);
  cout << nassoc-normal << " outlier events found\n";
  cout << "SD normal residence time: " << sd_2/1e3 << " ns\n";

 
return 0;
}

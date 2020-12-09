#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
using namespace std;

// Read input file
int main(int argc, char* argv[]) {
   string sinput, strfr, ligin, ligline, recin, recline, graline, selec_cut, siconc, slj_cut, sscaleLJ, sbeg, sfreq, send;
   if (argc == 11) { sinput = argv[1]; ligin = argv[2]; recin = argv[3]; selec_cut=argv[4]; siconc = argv[5]; slj_cut=argv[6]; sscaleLJ=argv[7];
                     sbeg=argv[8]; sfreq=argv[9]; send=argv[10]; }  
   else { cout << "Usage: ./E2.1 traj.pdb ligand.pqr receptor.pqr elec_cutoff ionic_conc (M) vdw_cutoff LJ_scale_coef (0 to 1) begin freq end" << endl;
          cout << "Incorrect number of input arguments. Proceeding with default settings (cutoff=12A, iconc=0.0M, scaleLJ=1.0)\n";
                   selec_cut="12"; slj_cut="12"; siconc="0.0"; sscaleLJ="1.0";
           cout << " ! Exiting. Could not assign traj ligand receptor " << endl;  
         return 0; } 

   int freq = atoi(sfreq.c_str()), beg = atoi(sbeg.c_str()), end = atoi(send.c_str());

// Find num frames and num ligand atoms
   ifstream input;
   int nframe = 0, tatom = 0;
   input.open(sinput.c_str(), ios::in);
     while (getline (input, strfr)) {
       if (strfr.find("END") != std::string::npos) { nframe++; }
         else { tatom++; }
       }
   int natom = tatom / nframe;
   cout << " >> Trajectory has " << nframe << " frames \n >> Ligand has " << natom << " atoms \n";
   input.close();  

   if (end > nframe) { end = nframe; cout << "!! Warning: End frame > than num frames in trajectory. Setting end = nframes.\n"; }

// Read ligand pqr for element types, LJ radii, and charges
   string element[natom];
   double ligcharge[natom], sigma[natom];
   int e = 0, q = 0, s = 0;
   ifstream ligand;
   ligand.open(ligin.c_str(), ios::in);
   while (getline (ligand, ligline)) {

    if (ligline.find("END") != std::string::npos) { cout<<"Hit END in pqr file. Done reading ligand pqr\n"; break; }
    element[e] = ligline.substr(13, 1).c_str();  // Atom element and charge store in order they appear in input pqr
    ligcharge[q] = atof(ligline.substr(61, 7).c_str());
    sigma[s] = atof(ligline.substr(69, 6).c_str());
    sigma[s] *= 2;
    e++; q++; s++;  // These counters may not be neccessary....
   }
   ligand.close(); 

// Assign LJ epsilon parameter by element
   string H ("H"), C ("C"), N ("N"), O ("O"), P ("P"), S ("S"), A ("A"), T ("T");   
   double epsilon[natom];
   int eps = 0;  // Why is this needed??
   for (int i = 0; i < natom; i++) {
     if (element[i].compare(H) == 0) { epsilon[i] = 0.0157; }
     else if (element[i].compare(C) == 0) { epsilon[i] = 0.0860; }
     else if (element[i].compare(N) == 0) { epsilon[i] = 0.1700; }
     else if (element[i].compare(O) == 0) { epsilon[i] = 0.2100; }
     else if (element[i].compare(P) == 0) { epsilon[i] = 0.2000; }
     else if (element[i].compare(S) == 0) { epsilon[i] = 0.2500; }
     else if (element[i].compare(T) == 0) { epsilon[i] = 0.5000; }
    }
   //cout << "The 10th element in your ligand is " << element[9] << " and its epsilon/sigma are " << epsilon[9] << "/" << sigma[9] << endl;

// Store ligand trajectory coords for each frame 
   string line;
   double rx[natom][nframe], ry[natom][nframe], rz[natom][nframe]; // 2D arrays hold natom*nframe coordinates

   ifstream traj;
   traj.open(sinput.c_str(), ios::in);
   int i = 0, j = 0;

   while (i < nframe) {
     getline (traj, line);
     
        if (line.substr(0, 4)=="ATOM") {  
            double x = atof(line.substr(30, 8).c_str());
            double y = atof(line.substr(38, 8).c_str());
            double z = atof(line.substr(46, 8).c_str());
            rx[j][i] = x; ry[j][i] = y; rz[j][i] = z;
            j++;
           }
        else if (line.substr(0, 3)=="CRY") { cout << "Ignoring crystal information\n"; continue; }
        else if (line.substr(0, 3)=="END") { i++; j=0; } // Frame is increased and atom num is reset to 0
    }
   traj.close();

// Store receptor coordinates from input pqr
   cout << "Reading receptor coordinates\n";
   int nrecatom = 0, ri = 0;
   ifstream getrecatoms;
   getrecatoms.open(recin.c_str(), ios::in);
   while (getline (getrecatoms, graline)) { 
      if (graline.substr(0, 4)=="ATOM"){ nrecatom++; } 
   }
   getrecatoms.close();
   
   double recx[nrecatom], recy[nrecatom], recz[nrecatom], reccharge[nrecatom], recepsilon[nrecatom], recsigma[nrecatom];
   string recelement[nrecatom];
   ifstream rec;
   rec.open(recin.c_str(), ios::in);

   while (getline (rec, recline)) {   
     if (recline.substr(0, 4)=="ATOM") {
        double rec_x = atof(recline.substr(31, 9).c_str());
        double rec_y = atof(recline.substr(41, 9).c_str());
        double rec_z = atof(recline.substr(51, 9).c_str());
        double rec_q = atof(recline.substr(61, 7).c_str());
        //cout << "getting recelement\n";
        recelement[ri] = recline.substr(13, 1).c_str();
        //cout<<"getting recsig\n";        
        recsigma[ri] = atof(recline.substr(69, 6).c_str());
        recsigma[ri] *= 2;
        recx[ri] = rec_x; recy[ri] = rec_y; recz[ri] = rec_z; // Can I directly assign the array values in the atof line?
        reccharge[ri] = rec_q;
        ri++;
       }

     else if (recline.find("END") != std::string::npos) {
       cout << "Done reading receptor pqr\n"; }
     else { "Error reading receptor pqr\n"; }
    }
   rec.close(); 
   cout << "assing rec LJ parm\n";
   for (int i = 0; i < nrecatom; i++) {
     if (recelement[i].compare(H) == 0) { recepsilon[i] = 0.0157; }
     else if (recelement[i].compare(C) == 0) { recepsilon[i] = 0.0860; }
     else if (recelement[i].compare(N) == 0) { recepsilon[i] = 0.1700; }
     else if (recelement[i].compare(O) == 0) { recepsilon[i] = 0.2100; }
     else if (recelement[i].compare(P) == 0) { recepsilon[i] = 0.2000; }
     else if (recelement[i].compare(S) == 0) { recepsilon[i] = 0.2500; }
     else if (recelement[i].compare(A) == 0) { recepsilon[i] = 0.7720; }  
     else if (recelement[i].compare(T) == 0) { recepsilon[i] = 0.5000; } 
    } 

// Calculating electrostatic and Lennard-Jones potentials for all ligand-receptor atom pairs    
   float Eelec, Evdw, elecfr[nframe], vdwfr[nframe];
   float Eelectot = 0, Evdwtot = 0, elec_cut = atof(selec_cut.c_str()), lj_cut = atof(slj_cut.c_str()), elec_trajtot=0, vdw_trajtot=0; 
    float iconc = atof(siconc.c_str());
   //double ions = 1. * iconc * 1e-27 * 6.02e23; // Default ionic conc= 0.002M, end converts to ions/A^3
    float Debye_length = sqrt( (78.5 * 8.85e-12 * 1.38e-23 * 298) / (1*1*2* 1.6e-19*1.6e-19 * iconc*1000 *6.02e23) );
    Debye_length *= 1e10; // Converts to Angstroms
    float scaleLJ = atof(sscaleLJ.c_str());
    int count = 0;
    float dist;
    float kappa = 1/Debye_length;
    cout << "Debye length = " << Debye_length << " Angstroms\n";

   cout << ">> Total ligand-receptor Eelec and Evdw for each frame\n           Eelec   |   Evdw\n";
   for (int frame = beg; frame < end; frame+=freq) {
     cout << "Frame " << frame << ": ";
     count++;
     for (int lai = 0; lai < natom; lai++) { // remind: lig coords stored as [atomnum][framenum]
      //cout<<"Showing E for "<<lai+1<<"th lig atom in frame "<<frame+1<<endl;
      for (int rai = 0; rai < nrecatom; rai++) {
        double dist = sqrt(pow(rx[lai][frame]-recx[rai], 2) + pow(ry[lai][frame]-recy[rai], 2) + pow(rz[lai][frame]-recz[rai], 2));
        if (dist <= elec_cut) {
           Eelec = ( 332.062 * ligcharge[lai]* reccharge[rai]* exp(-dist*kappa) )/ (78.5*dist);
           Eelectot += Eelec;
          } 
        if (dist <= lj_cut) {
           //double dist6 = pow(dist, 6);
           //double dist12 = dist6*dist6;
           //double sigma6 = pow(sigma[lai], 6); 
           //double sigma12 = sigma6*sigma6;
           //double Evdw = 4*epsilon[lai] * (sigma12/dist12 - sigma6/dist6);
           double epsij = sqrt(epsilon[lai]*recepsilon[rai]);
           double sigmaij = ((sigma[lai]+recsigma[rai]) / 2); //sigma is internuclear dist of atoms i & j at which LJ potential = 0
           Evdw = scaleLJ * 4 * epsij * (pow((sigmaij/dist), 12) - pow((sigmaij/dist), 6));
           Evdwtot += Evdw;
          }   
         }
        }
       elecfr[count-1] = Eelectot; vdwfr[count-1] = Evdwtot;      
       cout << Eelectot << "  |  " << Evdwtot << "\n";
       elec_trajtot += Eelectot; vdw_trajtot += Evdwtot;
       Eelectot = 0; Evdwtot = 0;  
       //if (frame == end-1) { break; }
       //cout << " * Beginning calculation for frame " << frame+2 << "\n";
     }

  double elecavg = elec_trajtot/count, vdwavg = vdw_trajtot/count, delec2, dvdw2, sumdelec2, sumdvdw2; 
  for (int i=0; i<count; i++) {
      delec2 = pow((elecfr[i]-elecavg), 2);
      dvdw2 =  pow((vdwfr[i]-vdwavg), 2);
      sumdelec2 += delec2; sumdvdw2 += dvdw2;  
      }
  cout << ">> " << count << " frames analyzed\n";
  double sdelec = sqrt(sumdelec2/(count-1)), sdvdw = sqrt(sumdvdw2/(count-1)); // Sample st. dev....
  cout << "Averages: " << elecavg << " +/- " << sdelec <<"  |  " << vdwavg <<" +/- " << sdvdw << endl;

  return 0;
}
      

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
using namespace std;

// Read user input
int main(int argc, char* argv[]) {
   string sdcd, strfr, ligin, ligline, recin, recline, graline, selec_cut, siconc, slj_cut, sscaleLJ, sbeg, sfreq, send, snrep;
   if (argc == 12) { sdcd = argv[1]; ligin = argv[2]; recin = argv[3]; selec_cut=argv[4]; siconc = argv[5]; slj_cut=argv[6]; sscaleLJ=argv[7];
                     sbeg=argv[8]; sfreq=argv[9]; send=argv[10]; snrep=argv[11]; }  
   else { cout << "Usage: ./E3.0 traj.dcd ligand.pqr receptor.pqr elec_cutoff ionic_conc (M) vdw_cutoff LJ_scale_coef (0 to 1) begin freq end n_replicates" << endl;
          cout << "Incorrect number of input arguments. Proceeding with default settings (cutoff=12A, iconc=0.0M, scaleLJ=1.0)\n";
                   selec_cut="12"; slj_cut="12"; siconc="0.0"; sscaleLJ="1.0";
           cout << " ! Exiting. Could not assign traj ligand receptor " << endl;  
         return 0; } 

   int beg = atoi(sbeg.c_str()); 
   int end = atoi(send.c_str());
   int nrep = atoi(snrep.c_str());

// Read DCD trajectory
   ifstream dcd;
   dcd.open(sdcd.c_str(), ios::in|ios::binary|ios::ate);
   long length;
   length = dcd.tellg();

   cout << "**Size of file is " << length/(1024*1024) << " MB.\n";
   char * buffer = new char[length];
   dcd.seekg(0, ios::beg);
   long headlen = 0;
   int tmp, natom, ifreq = atoi(sfreq.c_str());
   double freq = atof(sfreq.c_str());
   dcd.read((char*)&tmp, sizeof(int));
   headlen = headlen + 8 + tmp;
   dcd.seekg(headlen);
   dcd.read((char*)&tmp, sizeof(int));
   headlen = headlen + 8 + tmp;
   dcd.seekg(headlen);
   dcd.read((char*)&tmp, sizeof(int));
   headlen = headlen + 8 + tmp;
   dcd.read((char*)&natom, sizeof(int));

   int  nligatom = natom/nrep;
   long framelen = natom*12 + 80;
   long nframe = ((length-headlen)/framelen);
   cout << natom << " atoms in dcd file\n" << nligatom << " atoms in ligand\n" << nframe << " frames in dcd file\n";
   if (end > nframe) { end = nframe; }

   double cx, cy, cz; // Unit cell dimensions
   float coords[natom][3][nframe/ifreq];  
   // Coordinates stored as: frame 1 -> all x, all y, all z (ordered by atom index #), frame 2 -> all x, all y, all z, frame 3... etc
   cout << "** Reading ligand coordinates...\n";
    for (int i=0;i<nframe/ifreq;i++) {
       dcd.seekg(long(headlen) + framelen*i + 4);
       dcd.read((char*)&cx, sizeof(double));
       dcd.seekg(long(headlen) + framelen*i + 20);
       dcd.read((char*)&cy, sizeof(double));
       dcd.seekg(long(headlen) + framelen*i + 44);
       dcd.read((char*)&cz, sizeof(double));
       dcd.seekg(long(headlen) + framelen*i + 56);
       dcd.read((char*)coords, sizeof(float)); 
        for (int j=0;j<3;j++) {
         for (int k=0;k<natom;k++) {
         dcd.read((char*)&coords[k][j][i], sizeof(float));
        }
        dcd.read((char*)&tmp, sizeof(float));
        dcd.read((char*)&tmp, sizeof(float));
      }
     }
   dcd.close();
   delete[] buffer;

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
       
        recelement[ri] = recline.substr(13, 1).c_str();         
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
   double Eelec, Evdw, elecfr[nframe], vdwfr[nframe];
   double Eelectot = 0, Evdwtot = 0, elec_cut = atof(selec_cut.c_str()), lj_cut = atof(slj_cut.c_str()), elec_trajtot=0, vdw_trajtot=0; 
    float iconc = atof(siconc.c_str());
   cout << "iconc is " << iconc << endl;
   //double ions = 1. * iconc * 1e-27 * 6.02e23; // Default ionic conc= 0.002M, end converts to ions/A^3
    float Debye_length = sqrt( (78.5 * 8.85e-12 * 1.38e-23 * 298) / (1*1*2* 1.6e-19*1.6e-19 * iconc*1000 *6.02e23) ); 
  
    Debye_length *= 1e10; // Convert it into Angstroms
    float scaleLJ = atof(sscaleLJ.c_str());
    int count = 0;
    float dist;
    float kappa = 1/Debye_length;
    cout << "Debye length = " << Debye_length << " Angstroms\n"; 
    cout << ">> Total ligand-receptor Eelec and Evdw for each frame\n           Eelec   |   Evdw\n";
   for (int rep=0; rep<nrep; rep++) {
   cout << ">> Checking replicate " << rep+1 << " for adsorbed frames\n";
    for (int frame = beg; frame < end; frame+=freq) {
    float mindist = 100000, mincheck_dist = 100000;
  
     for (int b=0; b<nligatom; b++) { 
       float check_dist = coords[b+(rep*nligatom)][2][frame] - 15.4;
       if (check_dist < mincheck_dist) { mincheck_dist = check_dist; }
       } 
      if (mincheck_dist < 8) { 
     cout << "Frame " << frame << ": ";
     count++;
     for (int lai = 0; lai < nligatom; lai++) {
   
      for (int rai = 0; rai < nrecatom; rai++) {
        float dx2 = pow(coords[lai+(rep*nligatom)][0][frame]-recx[rai], 2);
        float dy2 = pow(coords[lai+(rep*nligatom)][1][frame]-recy[rai], 2);
        float dz2 = pow(coords[lai+(rep*nligatom)][2][frame]-recz[rai], 2);
             dist = sqrt(dx2 + dy2 + dz2);
             if (dist < mindist) { mindist = dist; }
        if (dist <= elec_cut) {
           Eelec = ( 332.062 * ligcharge[lai]* reccharge[rai]* exp(-dist*kappa) )/ (78.5*dist);
           Eelectot += Eelec;
          } 
        if (dist <= lj_cut) {
           double epsij = sqrt(epsilon[lai]*recepsilon[rai]);
           double sigmaij = ((sigma[lai]+recsigma[rai]) / 2); 
           Evdw = scaleLJ * 4 * epsij * (pow((sigmaij/dist), 12) - pow((sigmaij/dist), 6));
           Evdwtot += Evdw;
          }   
         }
        }   
       elecfr[count-1] = Eelectot; vdwfr[count-1] = Evdwtot;      
       cout << Eelectot << "  |  " << Evdwtot << "\nDistance: " << mindist << endl;
       elec_trajtot += Eelectot; vdw_trajtot += Evdwtot;
       Eelectot = 0; Evdwtot = 0; 
      }
      else { cout << "Not ads in frame " << frame << endl; continue; } 
     }
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
      

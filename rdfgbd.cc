#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

int main(int argc, char* argv[]) {

   string sdcd, srecpqr, infile, sbegframe, sendframe;
   if (argc == 6) { sdcd = argv[1]; srecpqr = argv[2]; infile = argv[3]; sbegframe = argv[4]; sendframe = argv[5];  }
      else { cout << ">> Please provide dcd trajectory, receptor.pqr, rdf input file, and first & last frame to be analyzed\n>> Usage: ./exec traj.dcd receptor.pqr rdf.in first_frame last_frame\n"; return 0; }

/// Read an input file with atom numbers of ref pts and prx atoms, etc. ///   
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

   int nrep = atoi(lines[1].c_str());
   int freq = atoi(lines[3].c_str());
   int nligatoms = atoi(lines[6].c_str());
   int ligrefatom = atoi(lines[8].c_str());
   int dcdfreq = atoi(lines[10].c_str());
   int ts = atoi(lines[12].c_str());   
   string indivrdfs = lines[14];
   cout << ">> Reading " << nrep << " ligand replicates in trajectory\n>> Reading every " << freq << " frames\n";

   inf.close();
   
/// Get number of atoms in receptor ///
   ifstream recpqr; recpqr.open(srecpqr.c_str(), ios::in); 
   string recline;
   int nrecatoms = 0;
   while (recpqr.is_open()) {
     getline (recpqr, recline);
     if (recpqr.eof()) { break; }
     if (recline.find("ATOM") != std::string::npos) { nrecatoms++; }
   }
   recpqr.close();
   cout << ">> " << nrecatoms << " atoms in receptor\n";

/// Read and store stationary receptor coordinates from pqr ///
   double rec[3][nrecatoms];
   ifstream recpqr2; recpqr2.open(srecpqr.c_str(), ios::in);

   for (int i=0; i<nrecatoms; i++) {
    getline (recpqr2, recline);
    rec[0][i] = atof(recline.substr(33, 7).c_str());
    rec[1][i] = atof(recline.substr(43, 7).c_str());
    rec[2][i] = atof(recline.substr(53, 7).c_str());
    }
    recpqr2.close();

/// Reading DCD trajectory ///
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
   float coords[natoms][3][nframe];
   for (int i=0;i<nframe;i++) {
     dcd.seekg(long(headlen) + framelen*i + 4); 
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
         }
      dcd.read((char*)&tmp, sizeof(float));
      dcd.read((char*)&tmp, sizeof(float));
      }
     }

   dcd.close();
   delete[] buffer;
   int bins[nrep][32]; // Create our bins
   for (int z=0; z<nrep; z++) { 
      for (int y=0; y<32; y++) { bins[z][y] = 0; } // Set all bin values to 0 first
    }
   int nmeasured = 0;
   int begframe = atoi(sbegframe.c_str());
   int endframe = atoi(sendframe.c_str());
   int assocframes[nrep];
   for (int i=0; i<nrep; i++) { assocframes[i] = 0; }

   if (endframe > nframe) { endframe = nframe; }   

   for (int j=begframe; j<endframe; j+=freq) {
     double totdist = 0;
     int nassoc = 0;
     //cout << "---------------------------\n** Reading frame " << j+1 << endl << "---------------------------\n";   
      //for (int k = 0; k<nrep; k++) {
      for (int k = 0; k<nrep; k++) {
       double mindist = 10000;
       nmeasured++;
        for (int m=0; m<nligatoms; m++) {
          for (int l=0; l<nrecatoms; l++) { 
            double dx2 = pow((rec[0][l]-coords[(k*nligatoms)+m][0][j]), 2);  // New way reads all lig atoms too, not just a COM ref atom...
            double dy2 = pow((rec[1][l]-coords[(k*nligatoms)+m][1][j]), 2);  // ...I think this is better for rigid body sims where its hard for lig...
            double dz2 = pow((rec[2][l]-coords[(k*nligatoms)+m][2][j]), 2);  // ...COM to get very close to receptor

            double dist = sqrt(dx2 + dy2 + dz2);
            if (dist < mindist) { mindist = dist; }        
          }
         }
     if (mindist <= 5) { 
        if (assocframes[k] == 0) { bins[k][0]++; nassoc++; assocframes[k] = j+1; } // Only record first frame of association
        else { nassoc++; bins[k][0]++; }
        }
     else if (mindist >  5 && mindist <= 10)  { bins[k][1] ++; }
     else if (mindist > 10 && mindist <= 15)  { bins[k][2] ++; }
     else if (mindist > 15 && mindist <= 20)  { bins[k][3] ++; }
     else if (mindist > 20 && mindist <= 25)  { bins[k][4] ++; }
     else if (mindist > 25 && mindist <= 30)  { bins[k][5] ++; }
     else if (mindist > 30 && mindist <= 35)  { bins[k][6] ++; }
     else if (mindist > 35 && mindist <= 40)  { bins[k][7] ++; }
     else if (mindist > 40 && mindist <= 45)  { bins[k][8] ++; }
     else if (mindist > 45 && mindist <= 50)  { bins[k][9] ++; }
     else if (mindist > 50 && mindist <= 55)  { bins[k][10] ++;}
     else if (mindist > 55 && mindist <= 60)  { bins[k][11] ++; }
     else if (mindist > 60 && mindist <= 65)  { bins[k][12] ++; }
     else if (mindist > 65 && mindist <= 70)  { bins[k][13] ++; }
     else if (mindist > 70 && mindist <= 75)  { bins[k][14] ++; }
     else if (mindist > 75 && mindist <= 80)  { bins[k][15] ++; }
     else if (mindist > 80 && mindist <= 85)  { bins[k][16] ++; }
     else if (mindist > 85 && mindist <= 90)  { bins[k][17] ++; }
     else if (mindist > 90 && mindist <= 95)  { bins[k][18] ++; }
     else if (mindist > 95 && mindist <= 100) { bins[k][19] ++; }
     else if (mindist > 100 && mindist <= 105){ bins[k][20] ++; }
     else if (mindist > 105 && mindist <= 110){ bins[k][21] ++; }
     else if (mindist > 110 && mindist <= 115){ bins[k][22] ++; }
     else if (mindist > 115 && mindist <= 120){ bins[k][23] ++; }
     else if (mindist > 120 && mindist <= 125){ bins[k][24] ++; }
     else if (mindist > 125 && mindist <= 130){ bins[k][25] ++; }
     else if (mindist > 130 && mindist <= 135){ bins[k][26] ++; }
     else if (mindist > 135 && mindist <= 140){ bins[k][27] ++; }
     else if (mindist > 140 && mindist <= 145){ bins[k][28] ++; }
     else if (mindist > 145 && mindist <= 150){ bins[k][29] ++; }
     else if (mindist > 150 && mindist <= 155){ bins[k][30] ++; }
     else if (mindist > 155 && mindist <= 160){ bins[k][31] ++; }
     totdist += mindist;
     //cout << "REP " << k+1 << " is " << mindist << "A away from DNA in frame " << j+1 << endl;
     }
     double avgdist = totdist/nrep;
   //cout<<">> Average ligand distance from receptor in frame "<<j+1<<" is "<<avgdist<<" A\n>> "<<nassoc<<" ligands are associated in frame "<<j+1<<endl;
   }

   if (indivrdfs == "YES") {
   for (int i=0; i<nrep; i++) { cout << "MOL " << i+1 << " rdf: " << endl;
       for (int j=0; j<32; j++) { cout << j*5 << " to " << (j+1)*5 << ": " << bins[i][j] << endl; }
     }
   }

   // Calculate 22 bins that hold sum of all MOLs in each bin
   int totbin[32];
   for (int t=0; t<32; t++) { totbin[t] = 0; } // Set all bins to 0 first
   for (int i=0; i<32; i++) { 
     for (int j=0; j<nrep; j++) { totbin[i] += bins[j][i]; }
    cout << ">> Bin " << i+1 << " total: " << totbin[i] << endl;
   }
   cout << "Total MOL distances measured: " << nmeasured << endl;
   double t_tot = 0, unassoc = 0;
   for (int i=0; i<nrep; i++) { 
      if (assocframes[i] == 0) { unassoc++; }
      else { cout << ">> Replicate " << i+1 << " associated in frame " << assocframes[i] << " in " << (assocframes[i]*ts*dcdfreq*freq)/1E6 << " ps\n"; } 
      t_tot += ((assocframes[i]*ts*dcdfreq*freq)/1E6);
      }

   double t_avg = t_tot/(nrep-unassoc);
   cout << ">> " << nrep-unassoc << " ligands associated to the receptor.\n>> Average association time: " << t_avg << " ps\n";

return 0;
}

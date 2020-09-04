#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

int main(int argc, char* argv[]) {

   string sdcd, srecpqr, infile, sbegframe, sendframe, sbinstart, sbinwidth;
   if (argc == 8) { sdcd = argv[1]; srecpqr = argv[2]; infile = argv[3]; sbegframe = argv[4]; sendframe = argv[5]; sbinstart = argv[6]; sbinwidth = argv[7];  }
      else { cout << ">> Please provide dcd trajectory, receptor.pqr, rdf input file, first & last frame to be analyzed, bin start (A), bin width (A)\n>> Usage: ./exec traj.dcd receptor.pqr rdf.in first_frame last_frame bin_start bin_width\n"; return 0; }

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
   //int ligrefatom = atoi(lines[8].c_str());
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
   float bins[nrep][30]; // Create our bins
   for (int z=0; z<nrep; z++) { 
      for (int y=0; y<30; y++) { bins[z][y] = 0; } // Set all bin values to 0 first
    }
   int nmeasured = 0; 
   float framesread = 0.0;
   double binstart = atof(sbinstart.c_str());
   double binwidth = atof(sbinwidth.c_str());
   cout << ">> Bin start radius: " << binstart << " A   Bin width: " << binwidth << " A\n";
   int begframe = atoi(sbegframe.c_str());
   int endframe = atoi(sendframe.c_str());
   int assocframes[nrep];
   for (int i=0; i<nrep; i++) { assocframes[i] = 0; }

   if (endframe > nframe) { endframe = nframe; }   

   for (int j=begframe; j<endframe; j+=freq) {
     double totdist = 0;
     int nassoc = 0;
     framesread++;
      for (int k = 0; k<nrep; k++) {
       float mindist = 10000., mindz = 10000., booldz = 0.;
       nmeasured++;
        for (int m=0; m<nligatoms; m++) {
          
            float dx2_1 = pow((60.7 - coords[(k*nligatoms)+m][0][j]), 2); 
            float dy2_1 = pow((64.0 - coords[(k*nligatoms)+m][1][j]), 2);
            float dx2_2 = pow((60.7 - coords[(k*nligatoms)+m][0][j]), 2);
            float dy2_2 = pow((59.0 - coords[(k*nligatoms)+m][1][j]), 2);
            float dx2_3 = pow((60.7 - coords[(k*nligatoms)+m][0][j]), 2);
            float dy2_3 = pow((54.0 - coords[(k*nligatoms)+m][1][j]), 2);  
            float xydist1 = sqrt(dx2_1 + dy2_1);
            float xydist2 = sqrt(dx2_2 + dy2_2);
            float xydist3 = sqrt(dx2_3 + dy2_3);
            if (xydist1 < mindist) { mindist = xydist1; }
            if (xydist2 < mindist) { mindist = xydist2; }
            if (xydist3 < mindist) { mindist = xydist3; }  // mindist ends up being the lowest of the 3 xy reference distances
            if (coords[(k*nligatoms)+m][2][j] > 98.73 || coords[(k*nligatoms)+m][2][j] < 25.16) { 
                booldz = 1;
                float dz;
                if( abs(coords[(k*nligatoms)+m][2][j] - 98.73) < abs(coords[(k*nligatoms)+m][2][j]-25.16)) { dz = abs((coords[(k*nligatoms)+m][2][j]-98.73)); }
                else { dz = abs((coords[(k*nligatoms)+m][2][j] - 25.16)); }
                if (dz < mindz) { mindz = dz; }                 
              }
           }
    
     if (booldz == 1) {  // Check if ligand was above or below the receptor
     //cout << "Mol " << k+1 << " had a dz in frame " << j+1 << endl << "mindz = " << mindz << " minxy = " << mindist << endl; 
        if (mindz > mindist) {  // If ligand was above or below receptor, and if z dist > than xy dist, bin according to z dist
  
        if (mindz <= binwidth) {
        if (assocframes[k] == 0) { bins[k][0]++; nassoc++; assocframes[k] = j+1; } // Only record first frame of association
        else { nassoc++; bins[k][0]++; }
        }
        else if (mindz > binwidth && mindz <= 2*binwidth)    { bins[k][1] ++; }
        else if (mindz > 2*binwidth && mindz <= 3*binwidth)  { bins[k][2] ++; }
        else if (mindz > 3*binwidth && mindz <= 4*binwidth)  { bins[k][3] ++; }
        else if (mindz > 4*binwidth && mindz <= 5*binwidth)  { bins[k][4] ++; }
        else if (mindz > 5*binwidth && mindz <= 6*binwidth)  { bins[k][5] ++; }
        else if (mindz > 6*binwidth && mindz <= 7*binwidth)  { bins[k][6] ++; }
        else if (mindz > 7*binwidth && mindz <= 8*binwidth)  { bins[k][7] ++; }
        else if (mindz > 8*binwidth && mindz <= 9*binwidth)  { bins[k][8] ++; }
        else if (mindz > 9*binwidth && mindz <= 10*binwidth)  { bins[k][9] ++; }
        else if (mindz > 10*binwidth && mindz <= 11*binwidth)  { bins[k][10] ++; }
        else if (mindz > 11*binwidth && mindz <= 12*binwidth)  { bins[k][11] ++; }
        else if (mindz > 12*binwidth && mindz <= 13*binwidth)  { bins[k][12] ++; }
        else if (mindz > 13*binwidth && mindz <= 14*binwidth)  { bins[k][13] ++; }
        else if (mindz > 14*binwidth && mindz <= 15*binwidth)  { bins[k][14] ++; }
        else if (mindz > 15*binwidth && mindz <= 16*binwidth)  { bins[k][15] ++; }
        else if (mindz > 16*binwidth && mindz <= 17*binwidth)  { bins[k][16] ++; }
        else if (mindz > 17*binwidth && mindz <= 18*binwidth)  { bins[k][17] ++; }
        else if (mindz > 18*binwidth && mindz <= 19*binwidth)  { bins[k][18] ++; }
        else if (mindz > 19*binwidth && mindz <= 20*binwidth)  { bins[k][19] ++; }
        else if (mindz > 20*binwidth && mindz <= 21*binwidth)  { bins[k][20] ++; }
        else if (mindz > 21*binwidth && mindz <= 22*binwidth)  { bins[k][21] ++; }
        else if (mindz > 22*binwidth && mindz <= 23*binwidth)  { bins[k][22] ++; }
        else if (mindz > 23*binwidth && mindz <= 24*binwidth)  { bins[k][23] ++; }
        else if (mindz > 24*binwidth && mindz <= 25*binwidth)  { bins[k][24] ++; }
        else if (mindz > 25*binwidth && mindz <= 26*binwidth)  { bins[k][25] ++; }
        else if (mindz > 26*binwidth && mindz <= 27*binwidth)  { bins[k][26] ++; }
        else if (mindz > 27*binwidth && mindz <= 28*binwidth)  { bins[k][27] ++; }
        else if (mindz > 28*binwidth && mindz <= 29*binwidth)  { bins[k][28] ++; }
        else if (mindz > 29*binwidth && mindz <= 30*binwidth)  { bins[k][29] ++; } 
        }
      else {  // If lig was above or below, but xy dist is was greater than the z dist, bin according to xy dist
         if (mindist <= binstart) {
         if (assocframes[k] == 0) { bins[k][0]++; nassoc++; assocframes[k] = j+1; } // Only record first frame of association
         else { nassoc++; bins[k][0]++; }
         }
        else if (mindist > binstart && mindist <= binstart+(1*binwidth))                 { bins[k][1] ++; }
        else if (mindist > binstart+(1*binwidth) && mindist <= binstart+(2*binwidth))    { bins[k][2] ++; }
        else if (mindist > binstart+(2*binwidth) && mindist <= binstart+(3*binwidth))    { bins[k][3] ++; }
        else if (mindist > binstart+(3*binwidth) && mindist <= binstart+(4*binwidth))    { bins[k][4] ++; }
        else if (mindist > binstart+(4*binwidth) && mindist <= binstart+(5*binwidth))    { bins[k][5] ++; }
        else if (mindist > binstart+(5*binwidth) && mindist <= binstart+(6*binwidth))    { bins[k][6] ++; }
        else if (mindist > binstart+(6*binwidth) && mindist <= binstart+(7*binwidth))    { bins[k][7] ++; }
        else if (mindist > binstart+(7*binwidth) && mindist <= binstart+(8*binwidth))    { bins[k][8] ++; }
        else if (mindist > binstart+(8*binwidth) && mindist <= binstart+(9*binwidth))    { bins[k][9] ++; }
        else if (mindist > binstart+(9*binwidth) && mindist <= binstart+(10*binwidth))   { bins[k][10] ++; } 
        else if (mindist > binstart+(10*binwidth) && mindist <= binstart+(11*binwidth))  { bins[k][11] ++; }
        else if (mindist > binstart+(11*binwidth) && mindist <= binstart+(12*binwidth))  { bins[k][12] ++; }
        else if (mindist > binstart+(12*binwidth) && mindist <= binstart+(13*binwidth))  { bins[k][13] ++; }
        else if (mindist > binstart+(13*binwidth) && mindist <= binstart+(14*binwidth))  { bins[k][14] ++; }
        else if (mindist > binstart+(14*binwidth) && mindist <= binstart+(15*binwidth))  { bins[k][15] ++; }
        else if (mindist > binstart+(15*binwidth) && mindist <= binstart+(16*binwidth))  { bins[k][16] ++; }
        else if (mindist > binstart+(16*binwidth) && mindist <= binstart+(17*binwidth))  { bins[k][17] ++; }
        else if (mindist > binstart+(17*binwidth) && mindist <= binstart+(18*binwidth))  { bins[k][18] ++; }
        else if (mindist > binstart+(18*binwidth) && mindist <= binstart+(19*binwidth))  { bins[k][19] ++; }
        else if (mindist > binstart+(19*binwidth) && mindist <= binstart+(20*binwidth))  { bins[k][20] ++; }
        else if (mindist > binstart+(20*binwidth) && mindist <= binstart+(21*binwidth))  { bins[k][21] ++; }
        else if (mindist > binstart+(21*binwidth) && mindist <= binstart+(22*binwidth))  { bins[k][22] ++; }
        else if (mindist > binstart+(22*binwidth) && mindist <= binstart+(23*binwidth))  { bins[k][23] ++; }
        else if (mindist > binstart+(23*binwidth) && mindist <= binstart+(24*binwidth))  { bins[k][24] ++; }
        else if (mindist > binstart+(24*binwidth) && mindist <= binstart+(25*binwidth))  { bins[k][25] ++; }
        else if (mindist > binstart+(25*binwidth) && mindist <= binstart+(26*binwidth))  { bins[k][26] ++; }
        else if (mindist > binstart+(26*binwidth) && mindist <= binstart+(27*binwidth))  { bins[k][27] ++; }
        else if (mindist > binstart+(27*binwidth) && mindist <= binstart+(28*binwidth))  { bins[k][28] ++; }
        else if (mindist > binstart+(28*binwidth))                                       { bins[k][29] ++; }
       }
     }
     
     else {  // Else if ligand wasn't above or below receptor, bin according to xy distance
     //cout << "Mol " << k+1 << " even w receptor in frame " << j+1 << endl << "mindz = " << mindz << " minxy = " << mindist << endl;   
     if (mindist <= binstart) {
          if (assocframes[k] == 0) { bins[k][0]++; nassoc++; assocframes[k] = j+1; } // Only record first frame of association
          else { nassoc++; bins[k][0]++; }
          }
        else if (mindist > binstart && mindist <= binstart+(1*binwidth))                 { bins[k][1] ++; }
        else if (mindist > binstart+(1*binwidth) && mindist <= binstart+(2*binwidth))    { bins[k][2] ++; }
        else if (mindist > binstart+(2*binwidth) && mindist <= binstart+(3*binwidth))    { bins[k][3] ++; }
        else if (mindist > binstart+(3*binwidth) && mindist <= binstart+(4*binwidth))    { bins[k][4] ++; }
        else if (mindist > binstart+(4*binwidth) && mindist <= binstart+(5*binwidth))    { bins[k][5] ++; }
        else if (mindist > binstart+(5*binwidth) && mindist <= binstart+(6*binwidth))    { bins[k][6] ++; }
        else if (mindist > binstart+(6*binwidth) && mindist <= binstart+(7*binwidth))    { bins[k][7] ++; }
        else if (mindist > binstart+(7*binwidth) && mindist <= binstart+(8*binwidth))    { bins[k][8] ++; }
        else if (mindist > binstart+(8*binwidth) && mindist <= binstart+(9*binwidth))    { bins[k][9] ++; }
        else if (mindist > binstart+(9*binwidth) && mindist <= binstart+(10*binwidth))   { bins[k][10] ++; }
        else if (mindist > binstart+(10*binwidth) && mindist <= binstart+(11*binwidth))  { bins[k][11] ++; }
        else if (mindist > binstart+(11*binwidth) && mindist <= binstart+(12*binwidth))  { bins[k][12] ++; }
        else if (mindist > binstart+(12*binwidth) && mindist <= binstart+(13*binwidth))  { bins[k][13] ++; }
        else if (mindist > binstart+(13*binwidth) && mindist <= binstart+(14*binwidth))  { bins[k][14] ++; }
        else if (mindist > binstart+(14*binwidth) && mindist <= binstart+(15*binwidth))  { bins[k][15] ++; }
        else if (mindist > binstart+(15*binwidth) && mindist <= binstart+(16*binwidth))  { bins[k][16] ++; }
        else if (mindist > binstart+(16*binwidth) && mindist <= binstart+(17*binwidth))  { bins[k][17] ++; }
        else if (mindist > binstart+(17*binwidth) && mindist <= binstart+(18*binwidth))  { bins[k][18] ++; }
        else if (mindist > binstart+(18*binwidth) && mindist <= binstart+(19*binwidth))  { bins[k][19] ++; }
        else if (mindist > binstart+(19*binwidth) && mindist <= binstart+(20*binwidth))  { bins[k][20] ++; }
        else if (mindist > binstart+(20*binwidth) && mindist <= binstart+(21*binwidth))  { bins[k][21] ++; }
        else if (mindist > binstart+(21*binwidth) && mindist <= binstart+(22*binwidth))  { bins[k][22] ++; }
        else if (mindist > binstart+(22*binwidth) && mindist <= binstart+(23*binwidth))  { bins[k][23] ++; }
        else if (mindist > binstart+(23*binwidth) && mindist <= binstart+(24*binwidth))  { bins[k][24] ++; }
        else if (mindist > binstart+(24*binwidth) && mindist <= binstart+(25*binwidth))  { bins[k][25] ++; }
        else if (mindist > binstart+(25*binwidth) && mindist <= binstart+(26*binwidth))  { bins[k][26] ++; }
        else if (mindist > binstart+(26*binwidth) && mindist <= binstart+(27*binwidth))  { bins[k][27] ++; }
        else if (mindist > binstart+(27*binwidth) && mindist <= binstart+(28*binwidth))  { bins[k][28] ++; }
        else if (mindist > binstart+(28*binwidth))                                       { bins[k][29] ++; } 
      }
     }
     //double avgdist = totdist/nrep;
   //cout<<">> Average ligand distance from receptor in frame "<<j+1<<" is "<<avgdist<<" A\n>> "<<nassoc<<" ligands are associated in frame "<<j+1<<endl;
   }

   if (indivrdfs == "YES") {
   for (int i=0; i<nrep; i++) { cout << "MOL " << i+1 << " rdf: " << endl;
       for (int j=0; j<30; j++) { cout << j*5 << " to " << (j+1)*5 << ": " << bins[i][j] << endl; }
     }
   }

   // Calculate 22 bins that hold sum of all MOLs in each bin
   float totbin[30];
   for (int t=0; t<30; t++) { totbin[t] = 0; } // Set all bins to 0 first
   for (int i=0; i<30; i++) { 
     for (int j=0; j<nrep; j++) { totbin[i] += bins[j][i]; }
    cout << ">> Bin " << i+1 << ": Total: " << totbin[i] << " ||  Average: " << totbin[i]/framesread <<  endl;
   }
   cout << "Total MOL distances measured: " << nmeasured << endl <<"Frames analyzed: " << framesread << endl;
   double t_tot = 0, unassoc = 0;
   for (int i=0; i<nrep; i++) { 
      if (assocframes[i] == 0) { unassoc++; }
      else { cout << ">> Replicate " << i+1 << " associated in frame " << assocframes[i] << " in " << (assocframes[i]*ts*dcdfreq*freq)/1000. << " ps\n"; } 
      t_tot += ((assocframes[i]*ts*dcdfreq*freq)/1000.);
      }

   double t_avg = t_tot/(nrep-unassoc);
   cout << ">> " << nrep-unassoc << " ligands associated to the receptor.\n>> Average association time: " << t_avg << " ps\n";

return 0;
}

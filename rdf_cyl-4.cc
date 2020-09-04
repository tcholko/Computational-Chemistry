#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

int main(int argc, char* argv[]) {

   string sdcd, infile, sbegframe, sendframe, sbinstart, sbinwidth;
   if (argc == 7) { sdcd = argv[1]; infile = argv[2]; sbegframe = argv[3]; sendframe = argv[4]; sbinstart = argv[5]; sbinwidth = argv[6]; }
      else { cout << ">> Please provide dcd trajectory, rdf input file, first frame, last frame, bin start, bin width\n>> Usage: ./exec traj.dcd rdf.in first_frame last_frame bin_start bin_width\n"; return 0; }
   cout << "User input: " << sdcd << " " << infile << " begin: " << sbegframe << " end: " << sendframe << " binstart: " << sbinstart << " bin width: " << sbinwidth << endl; 
///// Read an input file with atom numbers of ref pts and prx atoms, etc. /////   
   ifstream inf;
   inf.open(infile.c_str(), ios::in);
   string infline;
   int i=0; 
   double ref[10];// Number of ref pts is 10 for now

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

   int refbeg = atoi(lines[6].c_str());
   int refend = atoi(lines[8].c_str());
   int nligatoms = atoi(lines[10].c_str());
   int refnatom = (refend-refbeg)+1;
   cout << "refbeg= " << refbeg<<" refend= "<< refend<<"refnatom= "<<refnatom<<endl;
   int nmol = atoi(lines[1].c_str());
   int freq = atoi(lines[3].c_str());   
   cout << ">> " << nmol << " molecules in rdf\n>> Reading every " << freq << " frames\n";
   double xref1 = atof(lines[12].c_str());
   double xref2 = atof(lines[13].c_str());
   double xref3 = atof(lines[14].c_str());
   double yref1 = atof(lines[16].c_str());
   double yref2 = atof(lines[17].c_str());
   double yref3 = atof(lines[18].c_str());
   double topzref = atof(lines[20].c_str());
   double botzref = atof(lines[21].c_str());
   double mol[nmol];
   for (int j=0; j<nmol; j++) { 
     mol[j] = atof(lines[j+24].c_str());  // Reads in RDF molecule  atom numbers
     //cout << "mol " <<j+1<< " atom is " << mol[j] << endl; 
    }

   inf.close();

///// Reading DCD trajectory ////// 
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

   cout << ">> " << natoms << " atoms in dcd file\n>> " << nframe << " frames in dcd file\n";
  
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
   float bins[nmol][30]; // Create our bins
   for (int z=0; z<nmol; z++) { 
      for (int y=0; y<30; y++) { bins[z][y] = 0; } // Set all bin values to 0 first
    }
   int nmeasured = 0, analyzed = 0;
   double binstart = atof(sbinstart.c_str());
   double binwidth = atof(sbinwidth.c_str());
   int begframe = atoi(sbegframe.c_str());
   int endframe = atoi(sendframe.c_str());
   if (endframe > nframe) { endframe = nframe; }   

   int frame_accum_bin[endframe][30];
   int assocframes[nmol];
   for (int i=0; i<nmol; i++) { assocframes[i] = 0; }

   for (int j=begframe; j<endframe; j+=freq) {
     double nassoc = 0;
     double totdist = 0;
     //cout << "---------------------------\n** Reading frame " << j+1 << endl << "---------------------------\n";
     for (int k = 0; k<nmol; k++) { 
       float mindist = 10000., mindz = 10000., booldz = 0.;
       int molfirstatom = mol[k]-1;
       int ads_atoms = 0;   
       nmeasured++;
      for (int m=0; m<nligatoms; m++) {
 
        float dx2_1 = pow((xref1 - coords[molfirstatom+m][0][j]), 2);
        float dy2_1 = pow((yref1 - coords[molfirstatom+m][1][j]), 2);
        float dx2_2 = pow((xref2 - coords[molfirstatom+m][0][j]), 2);
        float dy2_2 = pow((yref2 - coords[molfirstatom+m][1][j]), 2);
        float dx2_3 = pow((xref3 - coords[molfirstatom+m][0][j]), 2);
        float dy2_3 = pow((yref3 - coords[molfirstatom+m][1][j]), 2);
        float xydist1 = sqrt(dx2_1 + dy2_1);
        float xydist2 = sqrt(dx2_2 + dy2_2);
        float xydist3 = sqrt(dx2_3 + dy2_3);
        if (xydist1 < mindist) { mindist = xydist1; }
        if (xydist2 < mindist) { mindist = xydist2; }
        if (xydist3 < mindist) { mindist = xydist3; }     

      if (coords[molfirstatom+m][2][j] > topzref || coords[molfirstatom+m][2][j] < botzref) { //119.42, 47.46
           booldz = 1;
           float dz;
           if( abs(coords[molfirstatom+m][2][j] - topzref) < abs(coords[molfirstatom+m][2][j]-botzref)) { dz = abs((coords[molfirstatom+m][2][j]-topzref)); }
            else { dz = abs((coords[molfirstatom+m][2][j] - botzref)); }
            if (dz < mindz) { mindz = dz; }
          }
      }

    if (booldz == 1) {  // Check if ligand was above or below the receptor

     if (mindz > mindist) {  // If ligand was above or below receptor, and if z dist > than xy dist, bin according to z dist
      if (mindz <= binwidth) {
        if (assocframes[k] == 0) { 
          for (int n=0; n<nligatoms; n++) { 
            for (int o=0; o<refnatom; o++) { 
              float ads_chk_x2 = pow((coords[refbeg-1+o][0][j] - coords[molfirstatom+n][0][j]), 2);
              float ads_chk_y2 = pow((coords[refbeg-1+o][1][j] - coords[molfirstatom+n][1][j]), 2);
              float ads_chk_z2 = pow((coords[refbeg-1+o][2][j] - coords[molfirstatom+n][2][j]), 2);
              float ads_chk    = sqrt(ads_chk_x2 + ads_chk_y2 + ads_chk_z2);
              if (ads_chk < 3.0) { ads_atoms++; }
              }
             }
           if ( ads_atoms < 2) {  bins[k][0]++; } 
           nassoc++; assocframes[k] = j+1; 
           } // Only record first frame of association
        else { 
            for (int n=0; n<nligatoms; n++) {
             for (int o=0; o<refnatom; o++) {
              float ads_chk_x2 = pow((coords[refbeg-1+o][0][j] - coords[molfirstatom+n][0][j]), 2);
              float ads_chk_y2 = pow((coords[refbeg-1+o][1][j] - coords[molfirstatom+n][1][j]), 2);
              float ads_chk_z2 = pow((coords[refbeg-1+o][2][j] - coords[molfirstatom+n][2][j]), 2);
              float ads_chk    = sqrt(ads_chk_x2 + ads_chk_y2 + ads_chk_z2);
              if (ads_chk < 3.0) { ads_atoms++; }
              }
             }
           if ( ads_atoms < 2) {  bins[k][0]++; }
           nassoc++; 
           }
       }
     else if (mindz > binwidth && mindz <= 2*binwidth)       { bins[k][1] ++; }
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
         if (assocframes[k] == 0) { 
           for (int n=0; n<nligatoms; n++) {
            for (int o=0; o<refnatom; o++) {
              float ads_chk_x2 = pow((coords[refbeg-1+o][0][j] - coords[molfirstatom+n][0][j]), 2);
              float ads_chk_y2 = pow((coords[refbeg-1+o][1][j] - coords[molfirstatom+n][1][j]), 2);
              float ads_chk_z2 = pow((coords[refbeg-1+o][2][j] - coords[molfirstatom+n][2][j]), 2);
              float ads_chk    = sqrt(ads_chk_x2 + ads_chk_y2 + ads_chk_z2);
              if (ads_chk < 3.0) { ads_atoms++; }
              }
             }
           if ( ads_atoms < 2) {  bins[k][0]++; }
           nassoc++; assocframes[k] = j+1; 
           }

         else { 
            for (int n=0; n<nligatoms; n++) {
            for (int o=0; o<refnatom; o++) {
              float ads_chk_x2 = pow((coords[refbeg-1+o][0][j] - coords[molfirstatom+n][0][j]), 2);
              float ads_chk_y2 = pow((coords[refbeg-1+o][1][j] - coords[molfirstatom+n][1][j]), 2);
              float ads_chk_z2 = pow((coords[refbeg-1+o][2][j] - coords[molfirstatom+n][2][j]), 2);
              float ads_chk    = sqrt(ads_chk_x2 + ads_chk_y2 + ads_chk_z2);
              if (ads_chk < 3.0) { ads_atoms++; }
              }
             }
           if ( ads_atoms < 2) {  bins[k][0]++; }
          nassoc++;  
          }
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
     if (mindist <= binstart) {
          if (assocframes[k] == 0) { 
           for (int n=0; n<nligatoms; n++) {
             for (int o=0; o<refnatom; o++) {
              float ads_chk_x2 = pow((coords[refbeg-1+o][0][j] - coords[molfirstatom+n][0][j]), 2);
              float ads_chk_y2 = pow((coords[refbeg-1+o][1][j] - coords[molfirstatom+n][1][j]), 2);
              float ads_chk_z2 = pow((coords[refbeg-1+o][2][j] - coords[molfirstatom+n][2][j]), 2);
              float ads_chk    = sqrt(ads_chk_x2 + ads_chk_y2 + ads_chk_z2);
              if (ads_chk < 3.0) { ads_atoms++; }
             }
            }
           if ( ads_atoms < 2) {  bins[k][0]++; }
         nassoc++; assocframes[k] = j+1; 
        }
        else { 
         for (int n=0; n<nligatoms; n++) {
            for (int o=0; o<refnatom; o++) {
              float ads_chk_x2 = pow((coords[refbeg-1+o][0][j] - coords[molfirstatom+n][0][j]), 2);
              float ads_chk_y2 = pow((coords[refbeg-1+o][1][j] - coords[molfirstatom+n][1][j]), 2);
              float ads_chk_z2 = pow((coords[refbeg-1+o][2][j] - coords[molfirstatom+n][2][j]), 2);
              float ads_chk    = sqrt(ads_chk_x2 + ads_chk_y2 + ads_chk_z2);
              if (ads_chk < 3.0) { ads_atoms++; }
              }
             }
           if ( ads_atoms < 2) { bins[k][0]++; }
        nassoc++; 
        }
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

            //totdist += mindist;
            //cout << "MOL " << k+1 << " is " << mindist << "A away from DNA in frame " << j+1 << endl;
     } // END mol k loop
  
   analyzed++;
     for (int w=0; w<30; w++) {
       for (int s=0; s<nmol; s++) { 
        frame_accum_bin[j][w] += bins[s][w]; 
       }
     cout << "frame " << j << " bin " << w << " accum = " << frame_accum_bin[j][w] << endl;
     }
          //double avgdist = totdist/nmol;
          //cout << ">> Average MOL distance from DNA in frame " <<j+1<< " is " << avgdist << " A\n";
    }  // END frame j loop

   //for (int i=0; i<nmol; i++) { cout << "MOL " << i+1 << " rdf: " << endl;
   //  for (int j=0; j<30; j++) { cout << j*5 << " to " << (j+1)*5 << ": " << bins[i][j] << endl; }
   //}
   // Calculate bins that hold sum of all MOLs in each bin
   float totbin[30];
   float binavg[30];
   for (int t=0; t<30; t++) { totbin[t] = 0; } // Set all bins to 0 first
   for (int i=0; i<30; i++) { 
     for (int j=0; j<nmol; j++) { totbin[i] += bins[j][i]; }
      binavg[i] = totbin[i]/analyzed;
     // cout << ">> Bin " << i+1 << " Total: " << totbin[i] << " | Average: " << binavg << endl;
     }  

   float frame_totbin[analyzed][30];

    for (int i=1; i<analyzed; i++) {
      for (int j=0; j<30; j++) { 
       frame_totbin[i][j] = frame_accum_bin[i][j] - frame_accum_bin[i-1][j]; 
       frame_totbin[0][j] = frame_accum_bin[0][j]; // In frame 1, frame tot = frame accumulation since its first frame
      }
     cout << "Frame " << i << " bin " << j << " has " << frame_totbin[i-1][j+1] << " mols\n"; 
    }
 
    for (int i=0; i<30; i++) {
     float sum_d2=0;
      for (int j=0; j<analyzed; j++) { 
      float frame_d2 = pow((frame_totbin[j][i] - binavg[i]), 2);
      sum_d2 += frame_d2;
     }
      float SD = sqrt(sum_d2/(analyzed-1));
      cout << ">> Bin " << i+1 << " Total: " << totbin[i] << " | Average: " << binavg[i] << "+/- " << SD << endl;
    }
    
   cout << "Total MOL distances measured: " << nmeasured << endl << "Frames analyzed: " << analyzed << endl;
   
return 0;
}

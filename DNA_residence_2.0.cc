#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

int main(int argc, char* argv[]) {

   string sdcd, infile, sbegframe, sendframe;
   if (argc == 5) { sdcd = argv[1]; infile = argv[2]; sbegframe = argv[3]; sendframe = argv[4];  }
      else { cout << ">> Please provide dcd trajectory, input file, and first and last frame to be analyzed\n>> Usage: ./exec traj.dcd restime.in first_frame last_frame\n"; return 0; }
  
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

   int refbeg = atoi(lines[6].c_str());
   int refend = atoi(lines[8].c_str());
   int refnatom = (refend-refbeg)+1;
   cout << "refbeg= " << refbeg<<" refend= "<< refend<<"refnatom= "<<refnatom<<endl;

   int nmol = atoi(lines[1].c_str());  
   int freq = atoi(lines[3].c_str());
   int ts  = atoi(lines[11].c_str()); 
   int nligatoms = atoi(lines[13].c_str());  
   cout << ">> " << nmol << " molecules in rdf\n>> Reading every " << freq << " frames, timestep is " <<ts<<" ps\n";
/// TO DO - Make this read in frist lig atom #, then count forward # of lig atoms for each one ///
   double mol[nmol];
   for (int j=0; j<nmol; j++) { 
     mol[j] = atof(lines[j+16].c_str());  // Reads in ligand identifier atom numbers
     //cout << "mol " <<j+1<< " identifier atom is " << mol[j] << endl; 
    }

   inf.close();

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

   int begframe = atoi(sbegframe.c_str());
   int endframe = atoi(sendframe.c_str());
   if (endframe > nframe) { endframe = nframe; }   

   int moln = 1, nassoc = 0, assoc, dissoctimer, assoctimer;
   double dist, totrestime = 0.0;
   for (int k = 0; k<nmol; k++) {
     int beginlig = mol[k]-1;
     assoc = 0; dissoctimer = 0; assoctimer = 0;
     double restimer = 0.0;  // Reset timer to 0 for each new ligand
 
      for (int j=begframe; j<endframe; j+=freq) {   
        double mindist = 10000.0;
         for (int m=beginlig; m<(beginlig+nligatoms); m++) {
            for (int l=refbeg-1; l<refend; l++) {
              double dx2 = pow((coords[l][0][j]-coords[m][0][j]), 2);
              double dy2 = pow((coords[l][1][j]-coords[m][1][j]), 2);
              double dz2 = pow((coords[l][2][j]-coords[m][2][j]), 2);
    
              dist = sqrt(dx2 + dy2 + dz2);
              if (dist < mindist) { mindist = dist; }

              if (dist <= 4.0 && assoc == 1) { restimer += (ts*freq); dissoctimer = 0; break; }  // Break out of loop so we don't over-count
              if (dist <= 4.0 && assoctimer*freq*ts >= 10) { cout << ">> MOL "<<moln<<" associated in frame "<< j <<" to receptor atom "<< l << endl;    
                 assoc = 1; nassoc++; dissoctimer = 0; assoctimer = 0; restimer += (3*ts*freq); break; }  // Add res time for the 3 uncounted frames
              if (dist <= 4.0 && assoc == 0) { cout << ">> MOL "<< moln << " is within association cutoff in frame "<< j << ". ++ timer\n";
                 assoctimer++;  break; }
              }
           if (assoc == 1 || assoctimer > 0) { break; }  // Break loop over ligand too as soon as we find association 
           }

        if (mindist > 4.0 && assoc == 1) { dissoctimer++; }  // Count the time ligand is outside association threshold
        if (mindist > 4.0 && dissoctimer*freq*ts >= 10) { totrestime += restimer;
           cout<<">> MOL "<<moln<<" dissociated at frame "<< j <<". Mindist = "<<mindist<< endl <<">> Last residence time: " <<restimer<<" ps\n";
           restimer = 0.0, dissoctimer = 0; assoc = 0; }  // Reset timer to 0 upon dissociation
        if (mindist > 4.0 && assoc == 0) { assoctimer = 0; }  // Forces association timer to count only consecutive frames
      }
     cout << ">> End of MOL "<< moln << " trajectory. Current residence time: " << restimer << " ps\n";  // Tell res time if still assoc at end
     totrestime += restimer;
     moln++; 
   }

   cout <<">> " << nassoc <<" association events detected.\n >> Average residence time: " << totrestime/nassoc << " ps\n";

return 0;
}

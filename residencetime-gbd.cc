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

   double nrep = atof(lines[1].c_str());
   double freq = atof(lines[3].c_str());
   int nligatoms = atoi(lines[6].c_str());
   //double dcdfreq = atoi(lines[8].c_str());
   double ts = atof(lines[8].c_str());
   double resdist = atof(lines[10].c_str());   
 
  cout << ">> Reading " << nrep << " ligand replicates in trajectory\n>> Reading every " << freq << " frames\n>> Simulation timestep is " << ts << " picoseconds\n";

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

   int begframe = atoi(sbegframe.c_str());
   int endframe = atoi(sendframe.c_str());
   double assoc, nassoc = 0, assoctimer, dissoctimer;
   float totrestime = 0.0;
   if (endframe > nframe) { endframe = nframe; }      
  
      for (int k = 0; k<nrep; k++) {
       cout << ">> Analyzing replicate " << k+1 << " trajectory...\n";
       assoc=0; assoctimer=0; dissoctimer=0;
       double restimer = 0;
       for (int j=begframe; j<endframe; j+=freq) {
       
       double mindist = 10000;
        for (int m=0; m<nligatoms; m++) {
          for (int l=0; l<nrecatoms; l++) {
            double dx2 = pow((rec[0][l]-coords[(k*nligatoms)+m][0][j]), 2);
            double dy2 = pow((rec[1][l]-coords[(k*nligatoms)+m][1][j]), 2);  
            double dz2 = pow((rec[2][l]-coords[(k*nligatoms)+m][2][j]), 2);  

            double dist = sqrt(dx2 + dy2 + dz2);
            if (dist < mindist) { mindist = dist; }        
        
            if (dist <= resdist && assoc == 1) { restimer += (ts*freq); dissoctimer = 0; break; }  // Break out of loop so we don't over-count
            if (dist <= resdist && assoctimer*freq*ts >= 20) { cout << ">> Rep "<< k+1 <<" associated in frame "<< j <<" to receptor atom "<< l+1 << endl;
                assoc = 1; nassoc++; dissoctimer = 0; assoctimer = 0; restimer += (2*ts*freq); break; }  // Add res time for the 3 uncounted frames
            if (dist <= resdist && assoc == 0) { cout << ">> Rep "<< k+1 << " is within association cutoff in frame "<< j << ". ++ timer\n";
                 assoctimer++;  break; }
              }
           if (assoc == 1 || assoctimer > 0) { break; }  // Break loop over ligand too as soon as we find association
           }

        if (mindist > resdist && assoc == 1) { dissoctimer++; }  // Count the time ligand is outside association threshold
        if (mindist > resdist && dissoctimer*freq*ts >= 20) { totrestime += restimer;
           cout<<">> Replicate "<< k+1 <<" dissociated at frame "<< j <<". Mindist = "<<mindist<< endl <<">> Last residence time: " <<restimer<<" ps\n";
           restimer = 0.0, dissoctimer = 0; assoc = 0; }  // Reset timer to 0 upon dissociation
        if (mindist > resdist && assoc == 0) { assoctimer = 0; }  // Forces association timer to count only consecutive frames

     }

     cout << ">> End of replicate "<< k+1 << " trajectory. Current residence time: " << restimer << " ps\n";  // Tell res time if still assoc at end
     totrestime += restimer;
   }

   cout <<">> " << nassoc <<" association events detected.\n >> Average residence time: " << totrestime/nassoc << " ps\n";

return 0;
}

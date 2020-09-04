#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <cstdio>
#include <sstream>
using namespace std;

int main(int argc, char* argv[]) {
   string temp, inpqr, inpdb;
   if (argc == 3) { inpqr = argv[1]; inpdb = argv[2]; }
      else { cout << "Command line error. Usage: ./exec gbd_traj.pdb pdb_in.pdb\n"; return 0; }

  ifstream pqrcount;
  pqrcount.open(inpqr.c_str(), ios::in);
 
  int j=0, k=0, count=0; 

  while (pqrcount.is_open()) { 
    getline (pqrcount, temp); count++; 
    if (pqrcount.eof()) { break; }
   }

  pqrcount.close();

  string line1[count];
  string pdbinfo[count];
  
  ifstream pqr;
  pqr.open(inpqr.c_str(), ios::in);
  
  while (pqr.is_open() && j < count) {
     getline (pqr, line1[j]);
     j++;
    }
  pqr.close();

  ifstream pdb;
  pdb.open(inpdb.c_str(), ios::in);

  while (pdb.is_open() && k < count) {
     getline (pdb, pdbinfo[k]);
    k++;
    }  
  pdb.close();

   for (int i=0; i<count; i++) { 
     istringstream iss(pdbinfo[i]);
     string col1, col2, col3, col4, col5;
     iss >> col1 >> col2 >> col3 >> col4 >> col5;
            
     istringstream pqrcol(line1[i]);
     string pcol1, pcol2, pcol3, pcol4, pcol5, pcol6, pcol7, pcol8, pcol9, pcol10;
     pqrcol >> pcol1 >> pcol2 >> pcol3 >> pcol4 >> pcol5 >> pcol6 >> pcol7 >> pcol8 >> pcol9 >> pcol10;   
     
          if (col2.length() == 1) { cout << col1<<"      "<< col2; } 
     else if (col2.length() == 2) { cout << col1<<"     "<< col2; }
     else if (col2.length() == 3) { cout << col1<<"    "<< col2; }
     else if (col2.length() == 4) { cout << col1<<"   "<< col2; }   

          if (col3.length() == 1) { cout <<"  "<< col3 <<"   "; }
     else if (col3.length() == 2) { cout <<"  "<< col3 <<"  "; }
     else if (col3.length() == 3) { cout <<"  "<< col3 <<" ";}
     else if (col3.length() == 4) { cout <<" "<< col3 <<" ";}

          if (col5.length() == 1) { cout << col4 <<"     "<< col5;}
     else if (col5.length() == 2) { cout << col4 <<"    "<< col5;}
     else if (col5.length() == 3) { cout << col4 <<"   "<< col5;}
     else if (col5.length() == 4) { cout << col4 <<"  "<< col5;}

          //if (pcol6.length()  < 10) { cout <<"      " << pcol6 <<"   "<< pcol7 <<"  "<<pcol8 <<"  "<< pcol9 <<"  "<< pcol10 << endl; break; }
          if (pcol6.length() == 14) { cout <<"      " << pcol6; }
     else if (pcol6.length() == 15) { cout <<"     " << pcol6; } 
     else if (pcol6.length() == 16) { cout <<"    " << pcol6; }

          if (pcol7.length() == 5) { cout <<"   "<< pcol7 <<"  "<<pcol8 <<"  "<< pcol9 <<"  "<< pcol10 << endl; }
     else if (pcol7.length() == 6) { cout <<"  "<< pcol7 <<"  "<<pcol8 <<"  "<< pcol9 <<"  "<< pcol10 << endl; }
     else if (pcol7.length() == 7) { cout <<" "<< pcol7 <<"  "<<pcol8 <<"  "<< pcol9 <<"  "<< pcol10 << endl; }
   }     
  return 0;
}
  

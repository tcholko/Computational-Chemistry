#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <cstdio>
using namespace std;

int main(int argc, char* argv[]) {
   string inf, line1, line2;
   if (argc == 2) { inf = argv[1]; }
      else { cout << "Command line error. Usage: ./exec in_filename\n"; return 0; }

   ifstream input;
   input.open(inf.c_str(), ios::in);
  int i = 1, count = 0;
  char delim('X');
   
  while (input.is_open()) {
   getline (input, line1, delim);
   getline (input, line2);
   
   cout << line1 << " " << i << " " << line2 << endl;
   count ++;
   if ( count % 18 == 0) { i++; }
   }
 
    


   
  // printf("ATOM      1  Na+ Na+     1    %8.3f%8.3f%8.3f", x, y, z);
  // printf("  0.00  0.00\n"); 
  // }

   return 0;
}
  

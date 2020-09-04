#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
using namespace std;

// Read input file
int main(int argc, char* argv[]) {
   string sinput, dx, dy, line;
   if (argc == 4) { sinput = argv[1]; dx = argv[2]; dy = argv[3]; } 
      else { cout << "Usage: ./executable input_file x_translation y_translation\n"; return 0; }

   double ddx = atof(dx.c_str());
   double ddy = atof(dy.c_str());

   ifstream input;
   input.open(sinput.c_str(), ios::in);

     while (getline (input, line)) {

      double x = atof(line.substr(0, 7).c_str());
      double y = atof(line.substr(8, 6).c_str());
      double z = atof(line.substr(15, 7).c_str());
      int atomnum = atoi(line.substr(23, 3).c_str());
      x += ddx;
      y += ddy;

      cout << x << " " << y << " " << z << " " << atomnum << " " << "7.0 OR ";
     }

   input.close();

return 0;
}


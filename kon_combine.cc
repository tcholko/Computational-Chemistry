#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
using namespace std;

int main(int argc, char* argv[]) {
   ifstream input, input2, input3;
   string sinput, sinput2, sinput3;
   string binders;
   if (argc == 5) { sinput = argv[1]; sinput2 = argv[2]; sinput3 = argv[3]; binders = argv[4]; }
    else { cout << "Please specify your input files and number of binders \nUsage: ./kon_combine input_file1 input_file2 input_file3 #_binders\n";
     return 0; 
    }

   int bind = atoi(binders.c_str());
   char delimiter('=');
   string garbage, stime[bind];
   double time[bind], time2[bind], time3[bind], tavg = 0, ttot = 0;
// Input 1
   int i = 0;
   input.open(sinput.c_str(), ios::in);
 
   while (input.is_open() && i < bind) {
     getline (input, garbage, delimiter);
     getline (input, stime[i]);
     time[i] = atof(stime[i].c_str());
     i++;
    }
    input.close();
// Input 2
    string stime2[bind];
    input2.open(sinput2.c_str(), ios::in);
    int j = 0;

   while (input2.is_open() && j < bind) {
     getline (input2, garbage, delimiter);
     getline (input2, stime2[j]);
     time2[j] = atof(stime2[j].c_str());
     j++;
    }
    input2.close();
// Input 3
    string stime3[bind];
    input3.open(sinput3.c_str(), ios::in);
    int h = 0;

   while (input3.is_open() && h < bind) {
     getline (input3, garbage, delimiter);
     getline (input3, stime3[h]);
     time3[h] = atof(stime3[h].c_str());
     h++;
    }
    input3.close();


   for (int k=0; k<bind; k++) {
     cout << time[k] << "\n" << time2[k] << "\n" << time3[k] <<"\n";
    }

return 0;
}

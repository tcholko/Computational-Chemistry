#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
using namespace std;

int main(int argc, char* argv[]) {
   ifstream input, input2, input3, input4, input5;
   string sinput, sinput2, sinput3, sinput4, sinput5, smax, holder;
   if (argc == 7) { sinput = argv[1]; sinput2 = argv[2]; sinput3 = argv[3]; sinput4 = argv[4]; sinput5 = argv[5]; smax = argv[6]; }
    else { cout << "Please specify your input files\nUsage: ./kon_combine input_file1 input_file2 input_file3... max_binders\n";
     return 0; 
    } 

   int i = 0;
   input.open(sinput.c_str(), ios::in);
 
   while (getline (input, holder)) { i++; }
   int nline1 = i, count = 0;
   string line1[nline1];
   cout << "file has " << i << " lines\n";
   input.clear();   
   input.seekg(0);   // These two lines seek back to beginning of file

   while (input.is_open() && count < nline1) {
    getline (input, line1[count]);
    count++;
   }
   input.close();
   //for (int j=0; j<nline; j++) { cout << line[j] << endl; }

   int j = 0;
   input2.open(sinput2.c_str(), ios::in);

   while (getline (input2, holder)) { j++; }
   int nline2 = j, count2 = 0;
   string line2[nline2];
   cout << "file has " << j << " lines\n";
   input2.clear();
   input2.seekg(0);   // These two lines seek back to beginning of file

   while (input2.is_open() && count2 < nline2) {
    getline (input2, line2[count2]);
    count2++;
   }
   input2.close();

   int l = 0;
   input3.open(sinput3.c_str(), ios::in);

   while (getline (input3, holder)) { l++; }
   int nline3 = l, count3 = 0;
   string line3[nline3];
   cout << "file has " << l << " lines\n";
   input3.clear();
   input3.seekg(0);   // These two lines seek back to beginning of file

   while (input3.is_open() && count3 < nline3) {
    getline (input3, line3[count3]);
    count3++;
   }
   input3.close();

   int m = 0;
   input4.open(sinput4.c_str(), ios::in);

   while (getline (input4, holder)) { m++; }
   int nline4 = m, count4 = 0;
   string line4[nline4];
   cout << "file has " << m << " lines\n";
   input4.clear();
   input4.seekg(0);   // These two lines seek back to beginning of file

   while (input4.is_open() && count4 < nline4) {
    getline (input4, line4[count4]);
    count4++;
   }
   input4.close();

   int n = 0;
   input5.open(sinput5.c_str(), ios::in);

   while (getline (input5, holder)) { n++; }
   int nline5 = n, count5 = 0;
   string line5[nline5];
   cout << "file has " << n << " lines\n";
   input5.clear();
   input5.seekg(0);   // These two lines seek back to beginning of file

   while (input5.is_open() && count5 < nline5) {
    getline (input5, line5[count5]);
    count5++;
   }
   input5.close();


   int maxlines = atoi(smax.c_str());

   for (int k=0; k<maxlines; k++) {
     if (k < nline1) { cout << line1[k] << endl; } 
     if (k < nline2) { cout << line2[k] << endl; }
     if (k < nline3) { cout << line3[k] << endl; }
     if (k < nline4) { cout << line4[k] << endl; }
     if (k < nline5) { cout << line5[k] << endl; }
   }
return 0;
}

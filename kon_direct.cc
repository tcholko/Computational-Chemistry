#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
using namespace std;

int main(int argc, char* argv[]) {
   ifstream input;
   string sinput, binders, svol, swindow;
   if (argc == 5) { sinput = argv[1]; binders = argv[2]; svol = argv[3]; swindow = argv[4]; }
    else { cout << "Specify your input file, # binders, system volume, and # k and t values to include in calculation window\nUsage: ./kon input_file #_binders vol(A^3) window\n";
     return 0; 
    }

   int bind = atoi(binders.c_str()), window = atoi(swindow.c_str());
   char delimiter('=');
   string line, stime[bind];
   input.open(sinput.c_str(), ios::in);
   double tavg = 0, ttot = 0;
   int i = 0, o = 0;
   double vol = atof(svol.c_str()), times[bind], kons[bind], tavgs[bind];   

   cout << " > assoc. time (ps)  |  t_avg  |  k_on \n -------------------------------------\n -------------------------------------\n";

   while (input.is_open() && i < bind) {
     getline (input, line);
    /* cout << garbage << endl;
     getline (input, stime[i]);
     cout << stime[i] << endl;
     double time = atof(stime[i].c_str());
     cout << time << endl; */
    if (line.find("event") != std::string::npos) {
     times[i] = atof(line.substr(23, 11).c_str()); 
     ttot += times[i];
     tavg = ttot/(i+1);
     tavgs[o] = tavg;
     o++;
     kons[i] = (vol * 6.02E23 * 1E-27 * 1E12)/tavg;
     cout << " > " << times[i] << "  |  " << tavg << "  |  " << kons[i] << "\n";
     i++;
     }
    else { continue; }
   }
   input.close();

   double sumd2 = 0, sd, relsd;
   for (int k=0; k<bind; k++) {
     double d = tavg - times[k];
     double d2 = d*d;
     sumd2 += d2;  
     sd = sqrt(sumd2/bind);
     relsd = (sd*100)/tavg;
      
    }
   cout << " >> Average association time = " << tavg << " +/- " << sd << " ps (" << relsd << "%)\n";
   
   double sumk = 0, sumdk2 = 0, sd_k, relsd_k, avgk, sumtavg = 0, sumdtavg2 = 0, sd_tavg, relsd_tavg, avgtavg;
   for (int j=0; j<bind; j++) { sumk += kons[j]; sumtavg += tavgs[j]; }
   avgk = sumk/(bind);
   avgtavg = sumtavg/bind;
   for (int l=bind-window; l<bind; l++) {
     double dk = avgk - kons[l];
     double dk2 = dk*dk;
     sumdk2 += dk2;
     double dtavg = avgtavg - tavgs[l];
     double dtavg2 = dtavg*dtavg;
     sumdtavg2 += dtavg2;
    }
   sd_k = sqrt(sumdk2/bind);
   relsd_k = (sd_k*100)/avgk;
   sd_tavg = sqrt(sumdtavg2/bind);
   relsd_tavg = (sd_tavg*100)/avgtavg;
   cout << ">> Standard deviation k = " << sd_k << " (" << relsd_k << "%)\n>> Standard deviation t_avg = " << sd_tavg << " (" << relsd_tavg << "%)" << endl;

return 0;
}

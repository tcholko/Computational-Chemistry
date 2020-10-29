import sys
import math 

if len(sys.argv) < 7:
   print("Usage: program.py tstep(ps) tint(frames over which to measure D) trest(frames between time origins) com_atom natoms nframes")

tstep = float(sys.argv[1]); tint = int(sys.argv[2]); trest = int(sys.argv[3]) 
com_atom = int(sys.argv[4]); natoms = int(sys.argv[5]); nframe = int(sys.argv[6]); 
nrest = int((nframe-tint)/trest)

#nframe = 100
#natoms= 32
#com_atom = 2
#tstep = 10 #time between each frame (ps)
xc = []
yc = []
zc = []


with open("prx_diff_10ps.pdb", "r") as file:
       for line in file:
          if line.startswith("ATOM"):
           # for i in range(nframe):
             # for j in range(natom):
                xc.append(line[30:37])
                yc.append(line[38:45])
                zc.append(line[46:53])
tdisp2 = 0
#inrest = int(nrest)
print(str(nrest) + " restarts will be used\n")
#trest = 20 #frames
#inrest = 4
#tint = 40 #frames
# Make multiple time origins
for i in range(nrest):
    x_disp = float(xc[ (tint*natoms) +(i*trest*natoms) + com_atom]) - float(xc[ (i*trest*natoms) + com_atom])
    y_disp = float(yc[ (tint*natoms) +(i*trest*natoms) + com_atom]) - float(yc[ (i*trest*natoms) + com_atom])
    z_disp = float(zc[ (tint*natoms) +(i*trest*natoms) + com_atom]) - float(zc[ (i*trest*natoms) + com_atom])
    disp2 = (x_disp**2 + y_disp**2 + z_disp**2)
    tdisp2 += disp2
    print(str(i))
avgdisp2 = tdisp2/nrest # Average over different time origins
print("Avereage Displacement in interval = " + str(round(avgdisp2, 8)) + " cm^2") 

D = avgdisp2/(2 * 3 *tint*tstep) 
print("D = " + str(round(D*1e12*1e-16, 8)) + " cm^2/s")

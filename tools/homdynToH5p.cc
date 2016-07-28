/*
  Read Homdyn Hxxx.DAT and convert it to the H5Part data
  format.

  g++  -I/Users/adelmann/install/hdf5-1.6.5/hdf5/include -I/Users/adelmann/svnwork/H5Part/src  -g -c homdynToH5p.cc
  g++ -o homdynToH5p homdynToH5p.o -L/Users/adelmann/svnwork/H5Part/src -lH5Part -L/Users/adelmann/install/hdf5-1.6.5/hdf5/lib -lhdf5 -lz   -lm


  Usage: homdynToH5p [-f newFilename]

  Reads HBUNCH.OUT and writes the data to HBUNCH.h5 or newFilename.h5 



*/

#include "H5Part.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc,char *argv[]){

  const int nCol = 25;
  const int nHeader = 1;
  int nLines  = 0;

  H5PartFile *file;

  double data[nCol][10000];
  string headers[nCol];
  string units[nCol];

  string fnStr("HBUNCH.OUT");

  ifstream in;    

  /*
  Open and read HOMDYN File:
  */

  in.open(fnStr.c_str());	

  /* 
    over read possible header
  */

  for (int l=0;l<nHeader;l++) {
    for (int c=0;c<nCol;c++) {	
     in >> headers[c];
     cout << c << " - - " << headers[c] << endl;
    }
   }
  headers[0] = string("SPOS"); // H5Root needs this name 
   /* 
     read in file data
   */

  while (1) {    
   for (int c=0;c<nCol;c++)
     in >> data[c][nLines];
   if (!in.good()) break;      
     nLines++;   
  }   
  
  in.close();

  cout  << "In HBUNCH.OUT found " << nLines << " lines " << endl;

  file=H5PartOpenFile("HBUNCH.h5",H5PART_WRITE);
  if(!file) {
    perror("File open failed:  exiting!");
    exit(0);
  }

  H5PartWriteFileAttribString(file,"File Description", "This file contains HOMDYN HBUNCH.OUT data");

  H5PartWriteFileAttribString(file,"tUnit","s");
  H5PartWriteFileAttribString(file,"xUnit","m");
  H5PartWriteFileAttribString(file,"yUnit","m");
  H5PartWriteFileAttribString(file,"zUnit","m");
  H5PartWriteFileAttribString(file,"pxUnit","#beta#gamma");
  H5PartWriteFileAttribString(file,"pyUnit","#beta#gamma");
  H5PartWriteFileAttribString(file,"pzUnit","#beta#gamma");
  H5PartWriteFileAttribString(file,"idUnit","1");
  H5PartWriteFileAttribString(file,"SPOSUnit","m");
  H5PartWriteFileAttribString(file,"TIMEUnit","s");
  H5PartWriteFileAttribString(file,"#gammaUnit","1");
  H5PartWriteFileAttribString(file,"ENERGYUnit","MeV");
  H5PartWriteFileAttribString(file,"#varepsilonUnit","m rad");
  H5PartWriteFileAttribString(file,"#varepsilonrUnit","m rad");

  H5PartWriteFileAttribString(file,"#varepsilonr-geomUnit","m rad");
  H5PartWriteFileAttribString(file,"RMSXUnit","m");
  H5PartWriteFileAttribString(file,"RMSRUnit","m");
  H5PartWriteFileAttribString(file,"RMSPUnit","#beta#gamma");

  H5PartWriteFileAttribString(file,"maxdEUnit","MeV");
  H5PartWriteFileAttribString(file,"max#phiUnit","deg");

  H5PartWriteFileAttribString(file,"phizUnit","deg");
  H5PartWriteFileAttribString(file,"enezUnit","keV");
  
/*
1       Z_[m] 
2	sigma_r_[mm] 
3	sigma_z_[mm] 
4	I_[A] 
5	sigma_x_[mm]         
6	enx_[um] 
7	sigma_y_[mm] 
8	eny_[um] 
9	T_[MeV] 
10	dg/g_[%] 
11	DE_[MeV]          
12	elz_[KeVmm] 	
13	Ez_[MV/m] 
14	Bz_[T] Bx_[G] By_[G] 	
15	Qgrad_[T/m] 
16	ByWig     
17	BHOR_[T] 
18	Time_[nsec] 
19	beta 
20	R/gL 
21	Lplas_[m] 
22	Zeq_[m] 
23	EWsteady
*/

  for (int t=0;t<nLines;t++)	{
    H5PartSetStep(file,t); /* must set the current timestep in file */

    double d3[3];
    double dummy = 0.0;
    int rc;
	
    rc = H5PartWriteStepAttrib(file,"SPOS", H5T_NATIVE_DOUBLE,&data[0][t],1);

    dummy = data[1][t]/1000.0;
    rc = H5PartWriteStepAttrib(file,"RMSR", H5T_NATIVE_DOUBLE,&dummy,1);
    rc = H5PartWriteStepAttrib(file,"TIME", H5T_NATIVE_DOUBLE,&data[17][t],1);
    rc = H5PartWriteStepAttrib(file,"ENERGY", H5T_NATIVE_DOUBLE,&data[8][t],1);

    d3[0] = data[4][t]/1000.0;
    d3[1] = data[6][t]/1000.0;
    d3[2] = 0.0;
    rc = H5PartWriteStepAttrib(file,"RMSX", H5T_NATIVE_DOUBLE,&d3,3);

    d3[0] = data[5][t]/1000000.0;
    d3[1] = data[7][t]/1000000.0;
    d3[2] = 0.0;
    rc = H5PartWriteStepAttrib(file,"#varepsilon", H5T_NATIVE_DOUBLE,&d3,3);


    d3[0] = 0.0;
    d3[1] = 0.0;
    d3[2] = 0.0;

    rc = H5PartWriteStepAttrib(file,"#gamma", H5T_NATIVE_DOUBLE,&dummy,1);
    rc = H5PartWriteStepAttrib(file,"#varepsilonr", H5T_NATIVE_DOUBLE,&d3,3);
    rc = H5PartWriteStepAttrib(file,"#varepsilonr-geom", H5T_NATIVE_DOUBLE,&d3,3);
    rc = H5PartWriteStepAttrib(file,"RMSP", H5T_NATIVE_DOUBLE,&d3,3);

    rc = H5PartWriteStepAttrib(file,"maxdE", H5T_NATIVE_DOUBLE,&dummy,1);
    rc = H5PartWriteStepAttrib(file,"max#phi", H5T_NATIVE_DOUBLE,&dummy,1);
    rc = H5PartWriteStepAttrib(file,"phiz", H5T_NATIVE_DOUBLE,&dummy,1);
    rc = H5PartWriteStepAttrib(file,"enez", H5T_NATIVE_DOUBLE,&dummy,1);
  }
  H5PartCloseFile(file);
} 


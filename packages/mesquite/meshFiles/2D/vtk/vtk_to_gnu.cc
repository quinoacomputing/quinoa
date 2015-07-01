//
// converts a vtk file to gnuplot format

#include <iostream.h>  
#include <fstream.h>  


int main()
{
  int m,n, pts_mn;
  cout << " Input: m,n for the input gnuplot file \n";
  cin >> m;
  cin >> n;
  pts_mn = (m+1)*(n+1);

  ifstream in_file("vtkdata");

  ofstream ou_file("gnuplotdata");

  // read the vtk file and store data
  const int length = 128;
  char xxx[length];

  // skip over the first 5 lines
  in_file.getline( xxx, length, '\n' );
  in_file.getline( xxx, length, '\n' );
  in_file.getline( xxx, length, '\n' );
  in_file.getline( xxx, length, '\n' );
  in_file.getline( xxx, length, '\n' );

  char coords[pts_mn][length];
   for(int i = 0; i < pts_mn; ++i)
   {
     in_file.getline( coords[i], length, '\n' );
   }

   for (int i=0; i<m+1; i++) {
     for (int j=0; j<n+1; j++) {
       ou_file << coords[j+i*(n+1)] << "\n";
     }
     ou_file << " \n";
   }

   for (int j=0; j<n+1; j++) {
     for (int i=0; i<m+1; i++) {
       ou_file << coords[j+i*(n+1)] << "\n";
     }
     ou_file << " \n";
   }

}

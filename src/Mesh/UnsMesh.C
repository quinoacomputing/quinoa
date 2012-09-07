//******************************************************************************
/*!
  \file      src/Mesh/UnsMesh.C
  \author    J. Bakosi
  \date      Fri 07 Sep 2012 12:37:51 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Unstructured mesh class definition
  \details   Unstructured mesh class definition
*/
//******************************************************************************

#include <UnsMesh.h>

using namespace Quinoa;

UnsMesh::UnsMesh()
//******************************************************************************
//  Constructor
//! \author J. Bakosi
//******************************************************************************
{
}

UnsMesh::~UnsMesh()
//******************************************************************************
//  Destructor
//! \author J. Bakosi
//******************************************************************************
{
}

// #include <cstdlib>
// #include <cstdio>
// #include <cmath>
// #include <string.h>
// #include "Mesh.h"
// #include "Macros.h"
// #include "Const.h"
// 
// static int nelem_g;
// static int *pg, *etags;
// 
// static void readmeshfile(int *npoin, int *nbpoin, int *nelem, double **coord,
//                          int **binpoel, int **inpoel, int **bptags, int **betags)
// // -----------------------------------------------------------------------------
// // Routine: readmeshfile - Read in mesh from file
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // reads meshfile, also puts some header info into the file specified by
// // OUT_MESH_FILENAME, where the reordered (and finally used) mesh will be
// // written, arguments:
// //   o npoin: the total the number of points (nodes) read in from the meshfile,
// //   o nbpoin: the total number of boundary points/elements read in from the
// //             meshfile,
// //   o nelem: the total number of domain (triangular) elements read in from the
// //            meshfile,
// //   o coord[npoin*NDIMN]: matrix that stores x,y coordinates of all points,
// //   o binpoel[nelem_g*NBNODE]: matrix that stores the global element
// //                              connectivity indices of boundary elements,
// //   o inpoel[nelem_g*NNODE]: matrix that stores the global element connectivity
// //                            indices of domain elements,
// //   o bptags[nbpoin*3]: matrix that stores the boundary tags assigned to
// //                       boundary points,
// //   o betags[nelem_g*2]: matrix that stores the boundary tags assigned to
// //                        boundary elements.
// //  the output arrays are allocated dynamically in this function, so the caller
// //  needs to take care of freeing them
// // -----------------------------------------------------------------------------
// {
//   FILE *imesh, *omesh, *opos3d, *igeo;
//   char str[STRLEN];
//   int meshversion, meshtype, meshdoublesize, e, p, etype, ntags, itmp, num, i,
//       q;
//   double dtmp;
// 
//   // start writing omesh...
//   if (!(omesh = fopen(OUT_MESH_FILENAME,"w")))
//     ERR("cannot open OUT_MESH_FILENAME!");
// 
//   // start reading mesh data in...
//   if (!(imesh = fopen(INP_MESH_FILENAME,"r")))
//     ERR("cannot open INP_MESH_FILENAME!");
// 
//   fgets( str, STRLEN, imesh );			// str << "$MeshFormat\n"
//   fprintf( omesh, "%s", str );			// copy to omesh
// 
//   // 2 0 8
//   fscanf( imesh, "%d %d %d\n", &meshversion, &meshtype, &meshdoublesize );
//   // copy to omesh
//   fprintf( omesh, "%d %d %d\n", meshversion, meshtype, meshdoublesize );
//   if ( (meshversion != 2) ||
//        (meshtype != 0) ||
//        (meshdoublesize != sizeof(double)) )
//     ERR("wrong meshformat");
// 
//   fgets( str, STRLEN, imesh );			// str << "$EndMeshFormat\n"
//   fprintf( omesh, "%s", str );			// copy to omesh
// 
//   fgets( str, STRLEN, imesh );			// str << "$Nodes\n"
//   fprintf( omesh, "%s", str );			// copy to omesh
// 
//   fscanf( imesh, "%d\n", npoin );		// read in number of points
//   fprintf( omesh, "%d\n", *npoin );		// copy to omesh
//   
//   // read in points...
//   // the node numbering and the order can be arbitrary in imesh and probably
//   // doesn't start with 0, so save the indices in Pg...
//   // array for storing point indices...
//   if ( !(pg = (int*)malloc((*npoin)*sizeof(int))) )
//     ERR("Can't allocate memory!");
// 
//   // array for storing x and y coordinates of vertices,
//   if ( !(*coord = (double*)malloc((*npoin)*NDIMN*sizeof(double))) )
//     ERR("Can't allocate memory!");
// 
//   for ( p = 0; p < *npoin; p++ ) {
//     // read in coordinates of point p and throw away z coordinate...
//     fscanf( imesh, "%d %lf %lf %lf\n", pg+p, *coord+p*NDIMN+0,
//                                        *coord+p*NDIMN+1, &dtmp );
//   }
// 
//   // check if numbering starts from 0...
//   for ( itmp = p = 0; p < *npoin; p++ ) if (pg[p] == 0) itmp = 1;
//   if (!itmp)
//     ERR("point numbering should start with 0, create a Point(0) in"
//         " the gmsh .geo file");
// 
//   fgets( str, STRLEN, imesh );	// str << "$EndNodes\n"
//   fgets( str, STRLEN, imesh );	// str << "$Elements\n"
//   // read in Nelem_g = number of boundary (line) elements +
//   //                   number of domain elements (triangles)
//   // (Nelem_g: total number of elements from gmsh)
//   fscanf( imesh, "%d\n", &nelem_g );
// 
//   // allocate maximum possible size for Etags, Betags, Binpoel, Inpoel,
//   // (the real sizes are certainly smaller than Nelem_g, however, we
//   // don't know this yet...
//   if ( !(etags = (int*)malloc(nelem_g*2*sizeof(int))) )
//     ERR("Can't allocate memory!");
//   if ( !(*betags = (int*)malloc(nelem_g*2*sizeof(int))) )
//     ERR("Can't allocate memory!");
//   if ( !(*binpoel = (int*)malloc(nelem_g*NBNODE*sizeof(int))) )
//     ERR("Can't allocate memory!");
//   if ( !(*inpoel = (int*)malloc(nelem_g*NNODE*sizeof(int))) )
//     ERR("Can't allocate memory!");
// 
//   // read in all Nelem_g elements, also count up total number of boundary
//   // elements/points (Nbpoin) and total number of domain elements (Nelem)
//   // (Nelem_g = Nbpoin + Nelem)...
//   for ( e = *nelem = *nbpoin = 0; e < nelem_g; e++ ) {
//     // element index, element type, number of tags
//     fscanf( imesh, "%d %d %d\n", &itmp, &etype, &ntags );
//     if ( ntags > MAXNTAGS ) ERR("ntags > MAXNTAGS in meshfile");
//     switch ( etype )
//     {
//       case 1: // line (boundary element with 2 nodes)...
// 	      // store physical entity index in betags
// 	      if ( ntags > 0 ) fscanf( imesh, "%d", *betags+(*nbpoin)*2+0 );
// 	      // store elementary geometrical entity index in betags
// 	      if ( ntags > 1 ) fscanf( imesh, "%d", *betags+(*nbpoin)*2+1 );
// 	      // throw away mesh partition index (if any)
// 	      if ( ntags > 2 ) fscanf( imesh, "%d", &itmp );
// 	      // store boundary element connectivity indices in binpoel...
//               fscanf( imesh, "%d %d\n", *binpoel+(*nbpoin)*NBNODE+0,
//                                         *binpoel+(*nbpoin)*NBNODE+1 );
// 	      // increase number of boundary elements/points...
// 	      (*nbpoin)++;
// 	      break;
// 
//       case 2: // triangle (surface element with 3 nodes)
// 	      // store physical entity index in etags
// 	      if ( ntags > 0 ) fscanf( imesh, "%d", etags+(*nelem)*2+0 );
// 	      // store elementary geometrical entity index in etags
// 	      if ( ntags > 1 ) fscanf( imesh, "%d", etags+(*nelem)*2+1 );
// 	      // throw away mesh partition index (if any)
// 	      if ( ntags > 2 ) fscanf( imesh, "%d", &itmp );
// 	      // store three node indices in inpoel...
//     	      fscanf( imesh, "%d %d %d\n", *inpoel+(*nelem)*NNODE+0,
//                                            *inpoel+(*nelem)*NNODE+1,
// 					   *inpoel+(*nelem)*NNODE+2 );
// 	      // increase number of domain elements...
// 	      (*nelem)++;
// 	      break;
// 
//       default: ERR("unknown element type");
//     }
//   }
// 
//   fclose( imesh );	// finish reading mesh data
// 
//   // generate boundary tags assigned to points instead of elements...
//   // array for containing the above: the first column is the global indices
//   // of the boundary points, the other two are the physical and elementary
//   // geometrical tags
//   if ( !(*bptags = (int*)malloc((*nbpoin)*3*sizeof(int))) )
//     ERR("Can't allocate memory!");
//   
//   for ( e = num = 0; e < *nbpoin; e++ )
//     for ( p = 0; p < NBNODE; p++ ) {
//       itmp = (*binpoel)[e*NBNODE+0];	// get next point index
//       
//       for ( i = q = 0; i < num; i++ )   // see if we have encountered it already
//          if ( (*bptags)[i*3+0] == itmp ) { q = 1; break; }
// 
//       if ( q == 0 ) {    // add new only if we haven't encountered it yet
//         (*bptags)[num*3+0] = itmp;		// copy global point index
//         (*bptags)[num*3+1] = (*betags)[e*2+0];	// copy physical index
//         (*bptags)[num*3+2] = (*betags)[e*2+1];	// copy geometrical index
//         num++;
//       }
//     }
// 
//   // start writing postprocess-file (specified by OUT_POSTPROCESS_3D_FILENAME)
//   // by copying file specified by INP_GEO_FILENAME in it...
//   if (!(opos3d = fopen(OUT_POSTPROCESS_3D_FILENAME,"w")))
//     ERR("cannot open OUT_POSTPROCESS_3D_FILENAME!");
// 
//   if (!(igeo = fopen(INP_GEO_FILENAME,"r")))
//     ERR("cannot open INP_GEO_FILENAME!");
//   
//   while ( fgets(str,STRLEN,igeo) ) fprintf( opos3d, "%s", str );
// 
//   fclose( igeo );	// finish reading mesh data
//   fclose( opos3d );	// close opos3d for now (cont'd during postprocess)
//   fclose( omesh );	// close omesh for now (cont'd in writemeshfile)
// }
// 
// static void writemeshfile(int npoin, int nelem, int nbpoin, double *coord,
//                           int *inpoel)
// // -----------------------------------------------------------------------------
// // Routine: writemeshfile - Write mesh to file
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   int p, e, cnt;
//   FILE *omesh;
// 
//   printf(" * outputing preprocessed mesh...\n");
//   fflush(stdout);
// 
//   // continue writing omesh...
//   if (!(omesh = fopen(OUT_MESH_FILENAME,"a")))
//     ERR("cannot open OUT_MESH_FILENAME!");
//   
//   // coordinates of point p...
//   for ( p = 0; p < npoin; p++ )
//     fprintf( omesh, "%d %g %g 0\n", pg[p], coord[p*NDIMN+0], coord[p*NDIMN+1] );
// 
//   fprintf( omesh, "$EndNodes\n" );
//   fprintf( omesh, "$Elements\n" );
//   fprintf( omesh, "%d\n", nelem_g-nbpoin );
// 
//   cnt = 0;
//   /*for ( e = 0; e < nbpoin; e++ )	// boundary elements...
//   {
//     // (gmsh) element index, element type = 1 (line), MAXNTAGS...
//     fprintf( omesh, "%d 1 %d ", cnt, MAXNTAGS );
// 
//     // tags...
//     // output physical entity index...
//     fprintf( omesh, "%d ", betags[e*2+0] );
//     // output elementary geometrical entity index...
//     fprintf( omesh, "%d ", betags[e*2+1] );
//     // output partition index...
//     fprintf( omesh, "0 " );
// 
//     // two boundary element numbers...
//     fprintf( omesh, "%d %d\n", binpoel[e*NBNODE+0], binpoel[e*NBNODE+1] );
// 
//     cnt++;	// increase (gmsh) element index counter
//   }*/
//   
//   for ( e = 0; e < nelem; e++ )		// triangle elements...
//   {
//     // (gmsh) element index, element type = 2 (triangle), MAXNTAGS...
//     fprintf( omesh, "%d 2 %d ", cnt, MAXNTAGS );
//    
//     // tags...
//     // output physical entity index...
//     fprintf( omesh, "%d ", etags[e*2+0] );
//     // output elementary geometrical entity index...
//     fprintf( omesh, "%d ", etags[e*2+1] );
//     // output partition index...
//     fprintf( omesh, "0 " );
//     
//     // three point numbers...
//     fprintf( omesh, "%d %d %d\n", inpoel[e*NNODE+0],
//                                   inpoel[e*NNODE+1],
//                                   inpoel[e*NNODE+2] );
// 
//     cnt++;	// increase (gmsh) element index counter
//   }
//   
//   fprintf( omesh, "$EndElements\n" );
// 
//   fclose( omesh );		// finish writing omesh
// }
// 
// static void gensup( int npoin, int nelem, int *inpoel,
// 	            int **esup1, int **esup2, int **psup1, int **psup2 )
// // -----------------------------------------------------------------------------
// // Routine: gensup - Generate derived mesh data structures
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // generates derived data structures for:
// //   o elements surrounding points (linked lists: esup1, esup2),
// //   o points surroinding points (linked lists: psup1, psup2),
// // input arguments:
// //   o npoin: number of points in mesh,
// //   o nelem: number of elements in mesh,
// //   o inpoel: domain-element connectivity,
// // output arguments:
// //   o esup1, esup2, psup1, psup2
// // output arrays are allocated dynamically here, caller should take
// // care of freeing them
// // -----------------------------------------------------------------------------
// {
//   int e, p, i, j, q, size;
//   int *lpoin;
// 
//   // calculate elements surrounding points: esup1, esup2...
//   // one of the linked lists for storing elements surrounding points: esup2
//   if ( !(*esup2 = (int*)malloc((npoin+1)*sizeof(int))) )
//     ERR("Can't allocate memory!");
//   
//   // calculate elements surrounding points: esup1, esup2...
//   // element pass 1: count number of elements connected to each point...
//   for ( i = 0; i < npoin+1; i++ )       // initialize esup2
//     (*esup2)[i] = 0;
//   for ( e = 0; e < nelem; e++ ) {
//     (*esup2)[inpoel[e*NNODE+0]+1]++;
//     (*esup2)[inpoel[e*NNODE+1]+1]++;
//     (*esup2)[inpoel[e*NNODE+2]+1]++;
//   }
// 
//   // storage/reshuffling pass 1: update storage counter and store,
//   // also find out the maximum size of esup1 (size)...
//   size = (*esup2)[0]+1;
//   for ( i = 1; i < npoin+1; i++ ) {
//      (*esup2)[i] += (*esup2)[i-1];
//      if ( (*esup2)[i]+1 > size ) size = (*esup2)[i]+1;
//   }
//   
//   // now we know size, so allocate the other one of the linked lists for
//   // storing elements surrounding points: esup1...
//   if ( !(*esup1 = (int*)malloc(size*sizeof(int))) )
//     ERR("Can't allocate memory!");
// 
//   // store the elements in esup1...
//   for ( e = 0; e < nelem; e++ ) {
//     i = inpoel[e*NNODE+0];
//     j = (*esup2)[i]+1;
//     (*esup2)[i] = j;
//     (*esup1)[j] = e;
// 
//     i = inpoel[e*NNODE+1];
//     j = (*esup2)[i]+1;
//     (*esup2)[i] = j;
//     (*esup1)[j] = e;
// 
//     i = inpoel[e*NNODE+2];
//     j = (*esup2)[i]+1;
//     (*esup2)[i] = j;
//     (*esup1)[j] = e;
//   }
// 
//   // storage/reshuffling pass 2...
//   for ( i = npoin; i > 0; i-- ) {
//     (*esup2)[i] = (*esup2)[i-1];
//   }
//   (*esup2)[0] = 0;
// 
//   // calculate points surrounding points: psup1, psup2...
//   // one of the linked lists for storing points surrounding points: psup2
//   if ( !(*psup2 = (int*)malloc((npoin+1)*sizeof(int))) )
//     ERR("Can't allocate memory!");
// 
//   // temporary array, used only in this function...
//   if ( !(lpoin = (int*)malloc(npoin*sizeof(int))) )
//     ERR("Can't allocate memory!");
// 
//   // determine size...
//   for ( i = 0; i < npoin; i++ )
//     lpoin[i] = 0;
// 
//   (*psup2)[0] = 0;
//   size = 0;
//   for ( p = 0; p < npoin; p++ ) {
//     for ( i = (*esup2)[p]+1; i <= (*esup2)[p+1]; i++ ) {
//       e = (*esup1)[i];
// 
//       q = inpoel[e*NNODE+0];
//       if ( (q != p) && (lpoin[q] != (p+1)) )
//       {	size++; lpoin[q] = p+1; }
// 
//       q = inpoel[e*NNODE+1];
//       if ( (q != p) && (lpoin[q] != (p+1)) )
//       { size++; lpoin[q] = p+1; }
// 
//       q = inpoel[e*NNODE+2];
//       if ( (q != p) && (lpoin[q] != (p+1)) )
//       { size++; lpoin[q] = p+1; }
//     }
//     (*psup2)[p+1] = size;
//   }
// 
//   // now we know Mpsup, so allocate the other one of the linked lists for
//   // storing points surrounding points: psup1
//   if ( !(*psup1 = (int*)malloc((size+1)*sizeof(int))) )
//     ERR("Can't allocate memory!");
// 
//   // ...and rerun the loop above to include "psup1[j] = q" as well...
//   for ( i = 0; i < npoin; i++ )
//     lpoin[i] = 0;
// 
//   (*psup2)[0] = 0;
//   j = 0;
//   for ( p = 0; p < npoin; p++ ) {
//     for ( i = (*esup2)[p]+1; i <= (*esup2)[p+1]; i++ ) {
//       e = (*esup1)[i];
// 	
//       q = inpoel[e*NNODE+0];
//       if ( (q != p) && (lpoin[q] != (p+1)) )
//       { j++; (*psup1)[j] = q; lpoin[q] = p+1; }
//       
//       q = inpoel[e*NNODE+1];
//       if ( (q != p) && (lpoin[q] != (p+1)) )
//       { j++; (*psup1)[j] = q; lpoin[q] = p+1; }
//       
//       q = inpoel[e*NNODE+2];
//       if ( (q != p) && (lpoin[q] != (p+1)) )
//       { j++; (*psup1)[j] = q; lpoin[q] = p+1; }
//     }
//     (*psup2)[p+1] = j;
//   }
// 
//   free( lpoin );
// }
// 
// static void regensup( int npoin, int nelem, int *inpoel,
// 	              int *esup1, int *esup2, int *psup1, int *psup2 )
// // -----------------------------------------------------------------------------
// // Routine: regensup - Regenerate derived mesh data structures
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // re-generates derived data structures for:
// //   o elements surrounding points (linked lists: esup1, esup2),
// //   o points surroinding points (linked lists: psup1, psup2),
// // input arguments:
// //   o npoin: number of points in mesh,
// //   o nelem: number of elements in mesh,
// //   o inpoel: domain-element connectivity,
// // output arguments:
// //   o esup1, esup2, psup1, psup2
// // -----------------------------------------------------------------------------
// {
//   int e, p, i, j, q;
//   int *lpoin;
// 
//   // recalculate elements surrounding points: esup1, esup2...
//   // element pass 1: count number of elements connected to each point...
//   for ( i = 0; i < npoin+1; i++ )       // initialize esup2
//     esup2[i] = 0;
//   for ( e = 0; e < nelem; e++ ) {
//     esup2[inpoel[e*NNODE+0]+1]++;
//     esup2[inpoel[e*NNODE+1]+1]++;
//     esup2[inpoel[e*NNODE+2]+1]++;
//   }
// 
//   // storage/reshuffling pass 1: update storage counter and store
//   for ( i = 1; i < npoin+1; i++ )
//     esup2[i] += esup2[i-1];
//   
//   // store the elements in esup1...
//   for ( e = 0; e < nelem; e++ ) {
//     i = inpoel[e*NNODE+0];
//     j = esup2[i]+1;
//     esup2[i] = j;
//     esup1[j] = e;
// 
//     i = inpoel[e*NNODE+1];
//     j = esup2[i]+1;
//     esup2[i] = j;
//     esup1[j] = e;
// 
//     i = inpoel[e*NNODE+2];
//     j = esup2[i]+1;
//     esup2[i] = j;
//     esup1[j] = e;
//   }
// 
//   // storage/reshuffling pass 2...
//   for ( i = npoin; i > 0; i-- )
//     esup2[i] = esup2[i-1];
//   esup2[0] = 0;
// 
//   // calculate points surrounding points: psup1, psup2...
//   // temporary array, used only in this function...
//   if ( !(lpoin = (int*)malloc(npoin*sizeof(int))) )
//     ERR("Can't allocate memory!");
// 
//   for ( i = 0; i < npoin; i++ )
//     lpoin[i] = 0;
// 
//   psup2[0] = 0;
//   j = 0;
//   for ( p = 0; p < npoin; p++ ) {
//     for ( i = esup2[p]+1; i <= esup2[p+1]; i++ ) {
//       e = esup1[i];
// 	
//       q = inpoel[e*NNODE+0];
//       if ( (q != p) && (lpoin[q] != (p+1)) )
//       { j++; psup1[j] = q; lpoin[q] = p+1; }
//       
//       q = inpoel[e*NNODE+1];
//       if ( (q != p) && (lpoin[q] != (p+1)) )
//       { j++; psup1[j] = q; lpoin[q] = p+1; }
//       
//       q = inpoel[e*NNODE+2];
//       if ( (q != p) && (lpoin[q] != (p+1)) )
//       { j++; psup1[j] = q; lpoin[q] = p+1; }
//     }
//     psup2[p+1] = j;
//   }
// 
//   free( lpoin );
// }
// 
// static void genbpsup( int nbpoin, int *binpoel, int *psup1, int *psup2,
//                       int **bpsup1, int **bpsup2, int **bpg )
// // -----------------------------------------------------------------------------
// // Routine: genbpsup - Generate derived boundary data structures
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // generates derived data structures for:
// //   o boundary points surrounding points (linked lists: bpsup1, bpsup2),
// //   o boundary->global index mapping array (bpg)
// // input arguments:
// //   o nbpoin: number of boundary points,
// //   o binpoel: boundary-element connectivity
// //   o psup1, psup2: linked lists storing points surrounding points,
// // output arguments:
// //   o bpsup1, bpsup2, bpg
// // -----------------------------------------------------------------------------
// {
//   int i, j, q, size;
// 
//   // array to store boundary->global index mapping array
//   if ( !(*bpg = (int*)malloc(nbpoin*sizeof(int))) )
//     ERR("Can't allocate memory!");
//   
//   // loop over boundary elements (equal amount of boundary points)
//   for ( size=1, i=0; i < nbpoin; i++ ) {
//     // generate boundary->global index mapping array
//     (*bpg)[i] = binpoel[i*NBNODE+0];
//     // count up total number of points surrounding boundary points
//     size += psup2[(*bpg)[i]+1]-psup2[(*bpg)[i]];
//   }
// 
//   // linked lists to store points surrounding boundary points
//   if ( !(*bpsup1 = (int*)malloc(size*sizeof(int))) )
//     ERR("Can't allocate memory!");
//   if ( !(*bpsup2 = (int*)malloc((nbpoin+1)*sizeof(int))) )
//     ERR("Can't allocate memory!");
// 
//   // generate points surrounding boundary points
//   (*bpsup1)[0] = (*bpsup2)[0] = 0;
//   for ( q=1, i=0; i < nbpoin; i++ ) { // loop over boundary points
//     // loop over points surrounding boundary point i
//     for ( j = psup2[(*bpg)[i]]+1; j <= psup2[(*bpg)[i]+1]; j++ )
//       // copy out indices of points surrounding boundary point i
//       (*bpsup1)[q++] = psup1[j];
//     // put in ending index of boundary point surrounding point indices
//     (*bpsup2)[i+1] = q-1;
//   }
// }
// 
// static void genesupel( int nelem, int *inpoel, int *esup1, int *esup2,
//                        int **esupel1, int **esupel2 )
// // -----------------------------------------------------------------------------
// // Routine: genesupel - Elements surrounding points of elements
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // generates derived data structure for unstructured Eulerian mesh for storing
// // the elements surrounding points of elements (linked lists: esupel1, esupel2)
// // -----------------------------------------------------------------------------
// {
//   int e, i, j, k, l, A, B, C, mesupel, mresupel, haveit;
//   int *resupel, *nresupel;
// 
//   // find out maximum size of "redundant" elements surrounding points of
//   // elements (resupel): mresupel
//   for ( mresupel = e = 0; e < nelem; e++ ) {
//     // get node-information of element e
//     A = inpoel[e*NNODE+0];
//     B = inpoel[e*NNODE+1];
//     C = inpoel[e*NNODE+2];
// 
//     // count up "redundant" number of elements surrounding all three nodes
//     // of element e
//     j = esup2[A+1]-esup2[A] + esup2[B+1]-esup2[B] + esup2[C+1]-esup2[C];
//     if ( j > mresupel ) mresupel = j;	// find out the biggest one
//   }
// 
//   // help array for storing indices of elements surrounding nodes of an element
//   // in a "redundant" way (resupel) and in a "non-redundant" way (nresupel),
//   // (both of size mresupel for worst case of nresupel)
//   if ( !(resupel = (int*)malloc(mresupel*sizeof(int))) )
//     ERR("Can't allocate memory!");
//   if ( !(nresupel = (int*)malloc(mresupel*sizeof(int))) )
//     ERR("Can't allocate memory!");
// 
//   // find out size of esupel1 (mesupel)
//   for ( mesupel = e = 0; e < nelem; e++ ) {
//     // get node-information of element e
//     A = inpoel[e*NNODE+0];
//     B = inpoel[e*NNODE+1];
//     C = inpoel[e*NNODE+2];
// 
//     // collect indices of elements surrounding all three nodes (redundant)
//     j = 0;
//     for ( i = esup2[A]+1; i <= esup2[A+1]; i++ ) resupel[j++] = esup1[i];
//     for ( i = esup2[B]+1; i <= esup2[B+1]; i++ ) resupel[j++] = esup1[i];
//     for ( i = esup2[C]+1; i <= esup2[C+1]; i++ ) resupel[j++] = esup1[i];
// 
//     // get rid of redundancy: resupel -> nresupel
//     nresupel[0] = resupel[0];	// the first one is certainly needed
//     for ( l = i = 1; i < j; i++ ) {
//        // check whether nresupel has resupel[i] already
//        for ( haveit = k = 0; k < l; k++ )
//          if (nresupel[k] == resupel[i]) {
//            haveit = 1;
//            k = l;
//          }
// 
//        if ( !haveit )	// put resupel[i] into nresupel if new
// 	 nresupel[l++] = resupel[i];
//     }
// 
//     // add up number of "non-redundant" elements surrounding nodes of elements
//     // -1, because we don't store element e itself (but it's in resupel)
//     mesupel += l-1;
//   }
// 
//   // linked lists for storing elements surrounding points of elements:
//   // esupel1, esupel2
//   if ( !(*esupel1 = (int*)malloc(mesupel*sizeof(int))) )
//     ERR("Can't allocate memory!");
//   if ( !(*esupel2 = (int*)malloc((nelem+1)*sizeof(int))) )
//     ERR("Can't allocate memory!");
// 
//   // finally, generate linked lists
//   (*esupel2)[0] = -1;
//   for ( e = 0; e < nelem; e++ ) {
//     // get node-information of element e
//     A = inpoel[e*NNODE+0];
//     B = inpoel[e*NNODE+1];
//     C = inpoel[e*NNODE+2];
// 
//     // collect indices of elements surrounding all three nodes (redundant)
//     j = 0;
//     for ( i = esup2[A]+1; i <= esup2[A+1]; i++ ) resupel[j++] = esup1[i];
//     for ( i = esup2[B]+1; i <= esup2[B+1]; i++ ) resupel[j++] = esup1[i];
//     for ( i = esup2[C]+1; i <= esup2[C+1]; i++ ) resupel[j++] = esup1[i];
// 
//     // get rid of redundancy: resupel -> nresupel
//     nresupel[0] = resupel[0];	// the first one is certainly needed
//     for ( l = i = 1; i < j; i++ ) {
//        // check whether nresupel has resupel[i] already
//        for ( haveit = k = 0; k < l; k++ )
//          if ( nresupel[k] == resupel[i] ) {
//            haveit = 1;
//            k = l;
//          }
// 
//        if ( !haveit )	// put resupel[i] into nresupel if new
// 	 nresupel[l++] = resupel[i];
//     }
// 
//     // store nresupel of element e in global linked lists (esupel1, esupel2),
//     // except element e itself
//     (*esupel2)[e+1] = (*esupel2)[e] + l - 1;    // -1: don't store element e
//     for ( k = i = 0; i < l; i++ )
//       if (nresupel[i] != e)
//         (*esupel1)[(*esupel2)[e]+(k++)+1] = nresupel[i];
//   }
// 
//   free( resupel );
//   free( nresupel );
// }
// 
// static void outputmeshinfo( int nelem, int npoin )
// // -----------------------------------------------------------------------------
// // Routine: outputmeshinfo - output some information about filenames and number
// //                           of elemente/points to stdout
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   printf( "----------------------------------------------------------------\n" );
//   printf( "INP_GEO_FILENAME = %s\n",
//            INP_GEO_FILENAME );
//   printf( "INP_MESH_FILENAME = %s\n",
//            INP_MESH_FILENAME );
//   printf( "OUT_MESH_FILENAME = %s\n",
//            OUT_MESH_FILENAME );
//   printf( "OUT_POSTPROCESS_3D_FILENAME = %s\n",
//            OUT_POSTPROCESS_3D_FILENAME );
//   printf( "OUT_POSTPROCESS_3DTAV_FILENAME = %s\n",
//            OUT_POSTPROCESS_3DTAV_FILENAME );
//   printf( "OUT_POSTPROCESS_OUTFLOW_FILENAME = %s\n",
//            OUT_POSTPROCESS_OUTFLOW_FILENAME );
//   printf( "OUT_POSTPROCESS_TOUTFLOW_FILENAME = %s\n",
//            OUT_POSTPROCESS_TOUTFLOW_FILENAME );
//   printf( "OUT_POSTPROCESS_TIME_FILENAME = %s\n",
//            OUT_POSTPROCESS_TIME_FILENAME );
//   printf( "OUT_POSTPROCESS_SURF_FILENAME = %s\n",
//            OUT_POSTPROCESS_SURF_FILENAME );
//   printf( "OUT_POSTPROCESS_PARTICLES_FILENAME = %s\n",
//            OUT_POSTPROCESS_PARTICLES_FILENAME );
//   printf( "PDF_BASE_FILENAME = %s\n",
//            PDF_BASE_FILENAME );
//   printf( "TPDF_BASE_FILENAME = %s\n",
//            TPDF_BASE_FILENAME );
//   printf( "STRIPE_BASE_FILENAME = %s\n",
//            STRIPE_BASE_FILENAME );
//   printf( "TSTRIPE_BASE_FILENAME = %s\n",
//            TSTRIPE_BASE_FILENAME );
//   printf( "DOWNSTREAM_EVOLUTION_FILENAME = %s\n",
//            DOWNSTREAM_EVOLUTION_FILENAME );
//   printf( "TDOWNSTREAM_EVOLUTION_FILENAME = %s\n",
//            TDOWNSTREAM_EVOLUTION_FILENAME );
//   printf( "STREAMWISE_CENTERLINE_FILENAME = %s\n",
//            STREAMWISE_CENTERLINE_FILENAME );
//   printf( "TSTREAMWISE_CENTERLINE_FILENAME = %s\n",
//            TSTREAMWISE_CENTERLINE_FILENAME );
//   printf( "WALL_FILENAME = %s\n",
//            WALL_FILENAME );
//   printf( "TWALL_FILENAME = %s\n",
//            TWALL_FILENAME );
//   printf( "AZ_BASE_FILENAME = %s\n",
//            AZ_BASE_FILENAME );
//   printf( "TAZ_BASE_FILENAME = %s\n",
//            TAZ_BASE_FILENAME );
//   #ifdef SAVERESTART
//   printf( "RESTART_FILENAME = %s\n",
//            RESTART_FILENAME );
//   #endif
//   printf( "npoin = %d\n", npoin );
//   printf( "nelem = %d\n", nelem );
//   printf("-----------------------------------------------------------------\n" );
//   fflush(stdout);
// }
// 
// static void correctordering(int npoin, int nelem, int nbpoin, int *inpoel,
//                             int *binpoel)
// // -----------------------------------------------------------------------------
// // Routine: correctordering - Correct mesh ordering
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // corrects the ordering, gmsh doesn't necessarily sort the numbering of the
// // nodes and/or elements, this is not usual, however, I have encountered it, so
// // I correct it here (if they're in the right order, nothing happens), new
// // versions of gmsh have a feature to order at save, so this should not be
// // necessary anymore, but it's fast, done only once in the beginning and doesn't
// // hurt
// // input arguments:
// //   o npoin: total number points in mesh,
// //   o nelem: total number of domain-elements in mesh,
// //   o nbpoin: total number of boundary-points/elements in mesh,
// // output arguments:
// //   o inpoel[nelem*NNODE]: the updated domain-element connectivity,
// //   o binpoel[nbpoin*NBNODE]: the updated boundary-element connectivity,
// // -----------------------------------------------------------------------------
// {
//   int maxpg, p, e;
//   int *ipg;
// 
//   printf( " * correcting the ordering of mesh...\n" );
//   fflush(stdout);
// 
//   // get maximum extent of pg...
//   for ( maxpg = p = 0; p < npoin; p++ )
//     if (maxpg < pg[p])
//       maxpg = pg[p];
//     
//   // temporary array for inverse mapping of pg
//   if ( !(ipg = (int*)malloc((maxpg+1)*sizeof(int))) )
//     ERR("Can't allocate memory!");
//    
//   // create inverse mapping vector and correct Pg...
//   for ( p = 0; p < npoin; p++ ) { 
//     ipg[pg[p]] = p;	// set inverse mapping element
//     pg[p] = p;		// the new order will be simply increasing order from 0
//   }
// 
//   // remap element connectivity...
//   for ( e = 0; e < nelem; e++ ) {
//     inpoel[e*NNODE+0] = ipg[inpoel[e*NNODE+0]];
//     inpoel[e*NNODE+1] = ipg[inpoel[e*NNODE+1]];
//     inpoel[e*NNODE+2] = ipg[inpoel[e*NNODE+2]];
//   }
//   
//   // remap boundary element connectivity, global boundary point indices...
//   for ( e = 0; e < nbpoin; e++ ) {
//     binpoel[e*NBNODE+0] = ipg[binpoel[e*NBNODE+0]];
//     binpoel[e*NBNODE+1] = ipg[binpoel[e*NBNODE+1]];
//   }
// 
//   free( ipg );
// }
// 
// static void renumber( int npoin, int nelem, int nbpoin, int *psup1, int *psup2,
//                       int *inpoel, int *binpoel, double *coord, int *bptags )
// // -----------------------------------------------------------------------------
// // Routine: renumber - Mesh renumber for better locality
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // renumber points for better locality with the advancing front technique
// // input arguments:
// //   o npoin: total number points in mesh,
// //   o nelem: total number of domain-elements in mesh,
// //   o nbpoin: the total number of boundary points/elements read in from the
// //             meshfile,
// //   o psup1, psup2: linked lists storing points surrounding points,
// // output arguments:
// //   o inpoel[nelem*NNODE]: the updated domain-element connectivity,
// //   o binpoel[nbpoin*NBNODE]: the updated boundary-element connectivity,
// //   o coord[npoin*NDIMN]: the update matrix that stores x,y coordinates of all
// //                         points,
// //   o bptags[nbpoin*3]: matrix that stores the boundary tags assigned to
// //                       boundary points.
// // -----------------------------------------------------------------------------
// {
//   int i, j, e, p, q; 
//   int num, cnt;
//   int *lpoin, *hpoin, *kpoin, *mapvec;
//   double x1, x2, y1, y2;
//   
//   printf( " * renumbering mesh...\n" );
//   fflush(stdout);
// 
//   // allocate memory
//   if ( !(lpoin = (int*)calloc(npoin, sizeof(int))) )
//     ERR( "Cannot allocate memory for lpoin");
//   if ( !(hpoin = (int*)calloc(npoin, sizeof(int))) )
//     ERR( "Cannot allocate memory for hpoin");
//   if ( !(kpoin = (int*)calloc(npoin, sizeof(int))) )
//     ERR( "Cannot allocate memory for kpoin");
//   if ( !(mapvec = (int*)calloc(npoin, sizeof(int))) )
//     ERR( "Cannot allocate memory for mapvec");
// 
//   // construct mapping with advancing front
//   hpoin[0]=0;  lpoin[0]=1;  mapvec[0]=0;  num=1;  //set first point
//   for ( i = 1; i < npoin; i++ ) hpoin[i] = -1;
//   while ( num < npoin ) {
//     cnt = 0;
//     i = 0;
//     for ( j = 0; j < npoin; j++ ) kpoin[j] = -1;
//     while ( (p = hpoin[i]) != -1 ) {
//       i++;
//       for ( j = psup2[p]+1; j <= psup2[p+1]; j++ ) {
//         q = psup1[j];
//         // loop only through those points that have not been counted yet
//         if ( lpoin[q] != 1 ) {
//           mapvec[q] = num++;
// 	  kpoin[cnt] = q;	// register the point as just been counted
// 	  lpoin[q] = 1;		// register the point as counted
// 	  cnt++;
//         }
//       }
//     }
//     memcpy( hpoin, kpoin, npoin*sizeof(int) );
//   }
// 
//   // apply mapping to boundary element conectivity and boundary point tag
//   // indices after renumbering
//   for ( e = 0; e < nbpoin; e++ ) {
//     binpoel[e*NBNODE+0] = mapvec[binpoel[e*NBNODE+0]];
//     binpoel[e*NBNODE+1] = mapvec[binpoel[e*NBNODE+1]];
//     bptags[e*3+0] = mapvec[bptags[e*3+0]];
//   }
// 
//   // apply mapping to element conectivity after renumbering
//   for ( e = 0; e < nelem; e++ ) {
//     inpoel[e*NNODE+0] = mapvec[inpoel[e*NNODE+0]];
//     inpoel[e*NNODE+1] = mapvec[inpoel[e*NNODE+1]];
//     inpoel[e*NNODE+2] = mapvec[inpoel[e*NNODE+2]];
//   }
// 
//   //  apply mapping to coordinates after renumbering
//   memset( lpoin, 0, npoin*sizeof(int) );	// clean lpoin
//   for (j = 0; j < npoin; j++) {
//     p = j;
//     if  (lpoin[p] == 0) {
//       x1 = coord[p*NDIMN+0];
//       y1 = coord[p*NDIMN+1];
//       while (lpoin[p] == 0) {
// 	x2 = coord[mapvec[p]*NDIMN+0];
// 	y2 = coord[mapvec[p]*NDIMN+1];
// 	coord[mapvec[p]*NDIMN+0] = x1;
// 	coord[mapvec[p]*NDIMN+1] = y1;
// 	x1 = x2; y1 = y2;
// 	lpoin[p] = 1;
// 	p = mapvec[p];
//       }
//     }
//   }
//   
//   // free memory
//   free( lpoin );
//   free( hpoin );
//   free( kpoin );
//   free( mapvec );
// }
// 
// static void precompute_element_data(int nelem, int npoin, int *inpoel,
//                                     double *coord, int *esup1, int *esup2,
//                                     double **Ae, double **dNx, double **dNy,
//                                     double **dete, double **sqrtAp,
//                                     double *minsqrtAp )
// // -----------------------------------------------------------------------------
// // Routine: precompute_element_data - Precompute element data
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // precalculates:
// //   * element sizes,
// //   * shapefunction derivatives,
// //   * determinant of elements
// // input arguments:
// //   o nelem: total number of domain-elements in mesh,
// //   o npoin: total number points in mesh,
// //   o inpoel[nelem*NNODE]: element connectivity,
// //   o coord[npoin*NDIMN]: point coordinates,
// //   o esup1[],esup2[]: linked lists of elements surrounding points
// // output arguments:
// //   o Ae[nelemv]: the element sizes,
// //   o dNx[nelem*NNODE]: x derivatives of nodal shapefunctions,
// //   o dNy[nelem*NNODE]: y derivatives of nodal shapefunctions,
// //   o dete[nelem]: determinant of element,
// //   o sqrtAp[npoin]: average of the square roots of the elements surrounding
// //                    each point.
// //   o minsqrtAp: the smallest characteristic element size
// // assigns characteristic element size for each point,
// // finds smallest characteristic element size
// // output arrays allocate here, caller should take care of freeing them
// // -----------------------------------------------------------------------------
// {
//   int i, n, e, eN, p, A, B, C, NA, NB, NC;
//   double xba, yca, xca, yba, Ae2;
//   double x[NNODE], y[NNODE];
// 
//   // for storing areas of elements, shapefunction derivatives, determinants
//   if ( !(*Ae = (double*)malloc(nelem*sizeof(double))) )
//     ERR("Can't allocate memory!");
//   if ( !(*dNx = (double*)malloc(nelem*NNODE*sizeof(double))) )
//     ERR("Can't allocate memory!");
//   if ( !(*dNy = (double*)malloc(nelem*NNODE*sizeof(double))) )
//     ERR("Can't allocate memory!");
//   if ( !(*dete = (double*)malloc(nelem*sizeof(double))) )
//     ERR("Can't allocate memory!");
//   
//   for ( e = 0; e < nelem; e++ ) {	// loop over all elements
//     eN = e*NNODE;
//     A = inpoel[eN+0];
//     B = inpoel[eN+1];
//     C = inpoel[eN+2];
//     NA = NDIMN*A;
//     NB = NDIMN*B;
//     NC = NDIMN*C;
//     xba = coord[NB+0] - coord[NA+0];
//     yca = coord[NC+1] - coord[NA+1];
//     xca = coord[NC+0] - coord[NA+0];
//     yba = coord[NB+1] - coord[NA+1];
//     Ae2 = xba*yca - xca*yba;
//     (*Ae)[e] = Ae2/2.0;			// size of element e
// 
//     // shapefunction derivatives...
//     (*dNx)[eN+0] = (yba - yca)/Ae2;	// NAx
//     (*dNx)[eN+1] = yca/Ae2;		// NBx
//     (*dNx)[eN+2] = -yba/Ae2;		// NCx
//     (*dNy)[eN+0] = (xca - xba)/Ae2;	// NAy
//     (*dNy)[eN+1] = -xca/Ae2;		// NBy
//     (*dNy)[eN+2] = xba/Ae2;		// NCy
// 
//     // determinant for element e
//     x[0] = coord[NA+0];  y[0] = coord[NA+1];
//     x[1] = coord[NB+0];  y[1] = coord[NB+1];
//     x[2] = coord[NC+0];  y[2] = coord[NC+1];
//     (*dete)[e] = x[0]*(y[1]-y[2]) + x[1]*(y[2]-y[0]) + x[2]*(y[0]-y[1]);
//   }
// 
//   // assign a characteristic element size to each point...
//   // for storing the characteristic element size
//   // (square-root of average element areas surrounding points)
//   if ( !(*sqrtAp = (double*)malloc(npoin*sizeof(double))) )
//     ERR("Can't allocate memory!");
//   
//   for ( p = 0; p < npoin; p++ )	{	// loop over all owned points
//     // loop over all elements surrounding point p
//     for ( (*sqrtAp)[p]=0.0, n=0, i=esup2[p]+1; i <= esup2[p+1]; n++, i++ )
//       (*sqrtAp)[p] += (*Ae)[esup1[i]];	// add up element sizes
// 
//     // calculate square-root of average element size
//     (*sqrtAp)[p] = sqrt((*sqrtAp)[p]/n);
//   }
//   
//   // find smallest sqrtAp (approx. smallest element size)
//   for ( (*minsqrtAp) = 1e+5, p = 0; p < npoin; p++ )
//      if ( (*sqrtAp)[p] < (*minsqrtAp) ) (*minsqrtAp) = (*sqrtAp)[p];
// }
// 
// static void genesuel(int nelem, int *inpoel, int *esupel1, int *esupel2,
//                      int **esuel)
// // -----------------------------------------------------------------------------
// // Routine: genesuel - Elements surrounding elements
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   int e, f, A, B, C;
// 
//   // array for storing the max three element indices surrounding elements
//   if ( !(*esuel = (int*)malloc(nelem*NNODE*sizeof(int))) )
//     ERR("Can't allocate memory!");
// 
//   // initialize esuel with -1 indicating no element on that side
//   for ( e = 0; e < nelem*NNODE; e++ )
//     (*esuel)[e] = -1;
// 
//   for ( e = 0; e < nelem; e++ )	{		// loop over all elements
//     A = inpoel[e*NNODE+0];			// get node indices
//     B = inpoel[e*NNODE+1];
//     C = inpoel[e*NNODE+2];
//     // loop over elements surround nodes of element e
//     for ( f = esupel2[e]+1; f <= esupel2[e+1]; f++ ) {
//       if ( ((A==inpoel[esupel1[f]*NNODE+0]) ||
//             (A==inpoel[esupel1[f]*NNODE+1]) ||
// 	    (A==inpoel[esupel1[f]*NNODE+2])) &&
//            ((B==inpoel[esupel1[f]*NNODE+0]) ||
//             (B==inpoel[esupel1[f]*NNODE+1]) ||
// 	    (B==inpoel[esupel1[f]*NNODE+2])) )
//          (*esuel)[e*NNODE+2] = esupel1[f];
//       
//       if ( ((B==inpoel[esupel1[f]*NNODE+0]) ||
//             (B==inpoel[esupel1[f]*NNODE+1]) ||
// 	    (B==inpoel[esupel1[f]*NNODE+2])) &&
//            ((C==inpoel[esupel1[f]*NNODE+0]) ||
//             (C==inpoel[esupel1[f]*NNODE+1]) ||
// 	    (C==inpoel[esupel1[f]*NNODE+2])) )
//          (*esuel)[e*NNODE+0] = esupel1[f];
//       
//       if ( ((C==inpoel[esupel1[f]*NNODE+0]) ||
//             (C==inpoel[esupel1[f]*NNODE+1]) ||
// 	    (C==inpoel[esupel1[f]*NNODE+2])) &&
//            ((A==inpoel[esupel1[f]*NNODE+0]) ||
//             (A==inpoel[esupel1[f]*NNODE+1]) ||
// 	    (A==inpoel[esupel1[f]*NNODE+2])) )
//          (*esuel)[e*NNODE+1] = esupel1[f];
//     }
//   }
// }
// 
// static void free_temp_arrays( void )
// // -----------------------------------------------------------------------------
// // Routine: free_temp_arrays - Free temporary arrays needed for mesh preprocessing
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   free( pg );
//   free( etags );
// }
// 
// void prepmsh(int *npoin, int *nbpoin, int *nelem, double **coord, int **bpg,
//              int **binpoel, int **inpoel, int **esup1, int **esup2, int **psup1,
//              int **psup2, int **bpsup1, int **bpsup2, int **esupel1,
//              int **esupel2, int **esuel, int **bptags, int **betags,
//              double **Ae, double **dNx, double **dNy, double **dete,
//              double **sqrtAp, double *minsqrtAp)
// // -----------------------------------------------------------------------------
// // Routine: prepmsh - Main mesh preprocessing function
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   // read meshfile
//   readmeshfile( npoin, nbpoin, nelem,				// <- modified
//                 coord, binpoel, inpoel, bptags, betags );	// <- allocated
// 
//   // output some information to stdout
//   outputmeshinfo( *nelem, *npoin );
//     
//   // correct the ordering
//   correctordering( *npoin, *nelem, *nbpoin,
//                    *inpoel, *binpoel );				// <- modified
// 
//   // generate derived data structures for elements and points surrounding points
//   gensup( *npoin, *nelem, *inpoel,
//           esup1, esup2, psup1, psup2 );				// <- allocated
// 
//   // renumber for better data locality
//   renumber( *npoin, *nelem, *nbpoin, *psup1, *psup2,
//             *inpoel, *binpoel, *coord, *bptags );		// <- modified
//   
//   // regenerate derived data structures for elements and points surrounding
//   // points after renumbering
//   regensup( *npoin, *nelem, *inpoel,
//             *esup1, *esup2, *psup1, *psup2 );			// <- modified
// 
//   // generate derived data structures for boundary points
//   genbpsup( *nbpoin, *binpoel, *psup1, *psup2,
//             bpsup1, bpsup2, bpg );				// <- allocated
// 
//   // generate derived data structure to store elements surrounding points of
//   // elements
//   genesupel( *nelem, *inpoel, *esup1, *esup2,
//              esupel1, esupel2 );				// <- allocated
// 
//   // generate derived data structure to store elements surrounding elements
//   genesuel( *nelem, *inpoel, *esupel1, *esupel2,
//             esuel );						// <- allocated
// 
//   // write out meshfile
//   writemeshfile( *npoin, *nelem, *nbpoin, *coord, *inpoel );
// 
//   // precompute element related data, shapefunction derivatives, etc.
//   precompute_element_data( *nelem, *npoin, *inpoel, *coord, *esup1, *esup2,
//                            Ae, dNx, dNy, dete, sqrtAp,		// <- allocated
// 			   minsqrtAp );				// <- modified
// 
//   // free temporary arrays needed for mesh preprocessing
//   free_temp_arrays();
// }

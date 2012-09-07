//******************************************************************************
/*!
  \file      src/Mesh/GmshReader.C
  \author    J. Bakosi
  \date      Fri 07 Sep 2012 03:55:50 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh reader class definition
  \details   Gmsh mesh reader class definition
*/
//******************************************************************************

#include <fstream>

#include <GmshReader.h>
#include <GmshException.h>

using namespace Quinoa;

void
GmshReader::open()
//******************************************************************************
//  Open Gmsh mesh file
//! \author J. Bakosi
//******************************************************************************
{
  m_mesh.open(m_filename, ifstream::in);
  if (!m_mesh.good()) throw GmshException(FATAL, FAILED_OPEN);
}

void
GmshReader::close()
//******************************************************************************
//  Close Gmsh mesh file
//! \author J. Bakosi
//******************************************************************************
{
  m_mesh.close();
  if (m_mesh.fail()) throw GmshException(WARNING, FAILED_CLOSE);
}

void
GmshReader::read(Mesh* mesh)
//******************************************************************************
//  Read Gmsh mesh from file
//! \author J. Bakosi
//******************************************************************************
{
}

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


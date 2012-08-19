// -----------------------------------------------------------------------------
// \file    src/Base/Const.h
// \author  jbakosi
// \date    The Aug 14 9:32:00 2012
// \brief   Global constants
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------

#include <cmath>
#include "Macros.h"

// flow and geometrical constants
// most of these must be the same as in the .geo file
const double RE = 3900.0;               // Reynolds number
const double DP = -1.0;                 // imposed mean pressure gradient
const double sqrtREp2 = sqrt(RE/2);

#ifdef WALLFUNCTIONS
const double D = 1.0+0.1;               // diameter of cylinder (1.0+YP)
// distance from walls where particle-boundary conditions are imposed
const double YP = 0.1*D;
#else
const double D = 1.0;                   // diameter of cylinder (1.0+YP)
const double YP = 0.0;  // elliptic relaxation: boundary layer fully resolved
#endif

// turbulence model constants
#ifdef WALLFUNCTIONS                    // needed only for wall-functions
const double C0 = 3.5;
const double CMU = 0.09;
const double KAPPA = 0.41;
const double E = 8.5;
#else                                   // needed only for elliptic relaxation
const double C1 = 1.85;
const double C2 = 0.63;
const double CT = 6.0;
const double CL = 0.134;
const double CXI = 1.0;
const double CETA = 72.0;
const double CV = 1.4;
const double GAMMA5 = 0.1;
#endif

const double C3 = 5.0;
const double C4 = 0.25;
const double COM1 = 0.5;                // Dreeben: 0.44
const double COM2 = 0.73;               // Dreeben: 0.9
const double Ct = 0.2;                  // micromixing model constant

#ifdef WALLFUNCTIONS
const double BETA = -2.0/(0.75*C0 + 0.5 + C3 + COM2 - COM1);
#endif

// numerical model constants
const int NPEL = 50;            // numer of particles per element (initially)
const int MINNPEL = 5;          // minimum number of particles per element
const int MAXTS = 147483647;    // maximum number of timesteps to simulate
const double RETIME = 0.0;      // time to start the scalar release
const double AVTIME = 60.0;     // time to start time-averaging
const double MAXTIME = 100.0;   // maximum time to simulate
const int MAX_IT = 1500;        // maximum number of iterations in linear solve
const double STOP_TOL = 1.0e-3; // stopping tolerance in linear solve
const double CFL = 0.8;         // CFL stability coefficient (Courant number)
const double Cp = 2.0e-3*CFL;   // pressure stability coefficient
const int NBI = 50;             // number of bins of sample space for outputing pdfs

#ifdef PROJECTION
const int CNBI = 5;             // number of bins of sample space for velocity conditioning
#else                           //  in VCIEM (using Fox's projection)
const int BINd = 5;             // number of bins/dimension of sample space for velocity conditioning
const int CNBI = BINd*BINd*BINd;//  in VCIEM (without Fox's projection) (DO NOT CHANGE CNBI, only BINd)
const int MINNPBI = 2;          // minimum numer of particles per conditioning bin
#endif

const double SD = YP+1.0e-3;    // a small distance, so that sampling locations are just inside the domain
const double XS = 3.0;          // release source location (x,y)
const double YS = 0.0+SD;
const double SS = 1.0;          // source strength
const double DS = D;            // diameter of source
const int ST2DUMP = 100;        // number of timesteps to be elapsed to dump data
const int NNODE = 3;            // number of nodes per element for scalar (3 = triangle)
const int NBNODE = 2;           // number of nodes per boundary element (2 = line)
const int MAXNTAGS = 3;         // maximum number of ntags belonging to each element in meshfile
const int NDIMN = 2;            // 2 spatial dimensions
const int RDOF = 9;             // degrees of freedom of elliptic relaxation tensor
const int U2DOF = 6;            // number of Reynolds stress tensor components retained
const int NURPP = 1;            // number of uniform random numbers per particle in table

#ifndef WALLFUNCTIONS
const int NGRPP = 6;            // number of Gaussian random numbers per particle in table (ell. relax)
#else
const int NGRPP = 4;            // number of Gaussian random numbers per particle in table (wall-func)
#endif // WALLFUNCTIONS

#ifdef ME
const int STAGES = 2;
const double ALPHA[] = {1/2,1}; // constants for 2-stage modified Euler timestepping
#else
const int STAGES = 1;           // Euler-Maruyama
#endif

const int STRLEN = 256;         // maximum length of c-style strings
const int MAXNTHREADS = 256;    // maximum number of threads
const double SMIN = -4.0;       // most negative value of sample space for outputing pdfs
const double SMAX = +4.0;       // most positive value of sample space for outputing pdfs
const double BINSIZE = (SMAX-SMIN)/NBI; // size of bin of sample space for outputing pdfs
const double EPS = 1.0e-14;     // machine zero
const double PI = 3.14159265358979323846;// pi
const double BOUND = 0.01;      // lower bound for quantities getting too close to zero

#ifndef WALLFUNCTIONS
const double CT2 = CT*CT;
const double RE3 = RE*RE*RE;
#endif

const double DSH = DS/2;
const double DSHS = DSH*DSH;
const int MAXDL = 11;           // number of downstream locations to output scross-stream statistics
const double DL[MAXDL] = { 1.06,// downstream locations to output cross-stream statistics
                           1.54,
                           2.02,
                           3.0,
                           4.0,
                           5.0,
                           6.0,
                           7.0,
                           8.0,
                          10.0,
                          16.0 };
const int MAXPL = 2;            // number of locations to output scalar pdfs at
const double PL[] = {  4.0, 0.0,// locations to output scalar pdfs at, pairs of (x/d, y/d)
		       8.0, 0.0 };
const int MAXAZ = 3;            // number of azimuthal locations to output statistics at
const double AZ[] = { 90.0001, 120.0, 150.0 }; // azimuthal locations to output statistics at (angle in degrees)


// filenames
// inputfile: original .geo file (not used directly, just copied over to postprocess file)
const char INP_GEO_FILENAME[] = "cylinder_g.geo";
// inputfile: original .msh file saved by gmsh (mesh version 2)
const char INP_MESH_FILENAME[] = "cylinder_g.msh";
// outputfile: renumbered and corrected (finally used mesh)
const char OUT_MESH_FILENAME[] = "cylinder.msh";
// outputfile: postprocess results file with time-averaged 3d functions at last timestep
const char OUT_POSTPROCESS_3DTAV_FILENAME[] = "cylinder_3dtav.pos";
// outputfile: postprocess results file with gmsh views of 3d functions at a given timestep
const char OUT_POSTPROCESS_3D_FILENAME[] = "cylinder_3d.pos";
// outputfile: profiles at outflow
const char OUT_POSTPROCESS_OUTFLOW_FILENAME[] = "dist";
// outputfile: time-averaged profiles at outflow
const char OUT_POSTPROCESS_TOUTFLOW_FILENAME[] = "dist.tav";
// outputfile: time-evolution of domain-average statistics
const char OUT_POSTPROCESS_TIME_FILENAME[] = "dist.t";
// outputfile: time-averaged statistics on the cylinder surface (lift, drag, etc.)
const char OUT_POSTPROCESS_SURF_FILENAME[] = "dist.surf";
// outputfile: time evolution of particle positions
const char OUT_POSTPROCESS_PARTICLES_FILENAME[] = "cylinder_particles.pos";
// restartfile for full restart: particle positions,velocity,frequency,scalar
const char RESTART_FILENAME[] = "cylinder.restart";
// outputfile: base filename for pdfs of scalar in selected locations
const char PDF_BASE_FILENAME[] = "pdf";
// outputfile: base filename for time-averged pdfs of scalar in selected locations
const char TPDF_BASE_FILENAME[] = "pdf.tav";
// outputfiles: base filename for scalar in points of selected stripes
const char STRIPE_BASE_FILENAME[] = "sdist";
// outputfiles: base filename for time-averaged scalar in points of selected stripes
const char TSTRIPE_BASE_FILENAME[] = "sdist.tav";
// outputfile: downstream evolution of different quantities
const char DOWNSTREAM_EVOLUTION_FILENAME[] = "dev";
// outputfile: downstream evolution of different time-averaged quantities
const char TDOWNSTREAM_EVOLUTION_FILENAME[] = "dev.tav";
// outputfile: streamwise centerline profiles
const char STREAMWISE_CENTERLINE_FILENAME[] = "dist.cen";
// outputfile: time-averaged streamwise centerline profiles
const char TSTREAMWISE_CENTERLINE_FILENAME[] = "dist.cen.tav";
// outputfile: wall profiles
const char WALL_FILENAME[] = "dist.wall";
// outputfile: time-averaged wall profiles
const char TWALL_FILENAME[] = "dist.wall.tav";
// outputfile: base filename for azimuthal profiles
const char AZ_BASE_FILENAME[] = "dist.az";
// outputfile: base filename for time-averaged azimuthal profiles
const char TAZ_BASE_FILENAME[] = "dist.az.tav";


//////////////////////////////
// labels in postprocess files
//////////////////////////////
const char VIEWNAME_MEANVELOCITY_U[] = "<U>";	// name of view for x mean velocity
const char VIEWNAME_MEANVELOCITY_V[] = "<V>";	// name of view for y mean velocity
const char VIEWNAME_REYNOLDS11[] = "<uu>";	// name of view for Reynolds-stress component <uu>
const char VIEWNAME_REYNOLDS22[] = "<vv>";	// name of view for Reynolds-stress component <vv>
const char VIEWNAME_REYNOLDS33[] = "<ww>";	// name of view for Reynolds-stress component <ww>
const char VIEWNAME_REYNOLDS12[] = "<uv>";	// name of view for Reynolds-stress component <uv>
const char VIEWNAME_REYNOLDS13[] = "<uw>";	// name of view for Reynolds-stress component <uw>
const char VIEWNAME_REYNOLDS23[] = "<vw>";	// name of view for Reynolds-stress component <vw>
const char VIEWNAME_SKEWNESS11[] = "S(U)";	// name of view for skewness of streamwise velocity
const char VIEWNAME_SKEWNESS22[] = "S(V)";	// name of view for skewness of cross-stream velocity
const char VIEWNAME_FLATNESS11[] = "F(U)";	// name of view for flatness of streamwise velocity
const char VIEWNAME_FLATNESS22[] = "F(V)";	// name of view for flatness of cross-stream velocity
const char VIEWNAME_DISSIPATION[] = "eps";	// name of view for rate of dissipation of turbulent kinetic energy
const char VIEWNAME_TKE[] = "tke";		// name of view for turbulent kinetic energy
const char VIEWNAME_MEAN_PRESSURE[] = "<P>";	// name of view for mean pressure
const char VIEWNAME_MEANVELOCITY[] = "<Ui>";	// name of view for mean velocity vectorfield
const char VIEWNAME_VORTICITY[] = "<Wz>";	// name of view for mean spanwise vorticity
#ifndef WALLFUNCTIONS
const char VIEWNAME_WIGGLY_P11[] = "<wp11>";	// name of view for wiggly-p11
const char VIEWNAME_WIGGLY_P22[] = "<wp22>";	// name of view for wiggly-p22
const char VIEWNAME_WIGGLY_P33[] = "<wp33>";	// name of view for wiggly-p33
const char VIEWNAME_WIGGLY_P12[] = "<wp12>";	// name of view for wiggly-p12
const char VIEWNAME_WIGGLY_P13[] = "<wp13>";	// name of view for wiggly-p13
const char VIEWNAME_WIGGLY_P23[] = "<wp23>";	// name of view for wiggly-p23
const char VIEWNAME_WIGGLY_P21[] = "<wp21>";	// name of view for wiggly-p21
const char VIEWNAME_WIGGLY_P31[] = "<wp31>";	// name of view for wiggly-p31
const char VIEWNAME_WIGGLY_P32[] = "<wp32>";	// name of view for wiggly-p32
#endif
const char VIEWNAME_MEAN_SCALAR[] = "<C>";	// name of view for mean of scalar
const char VIEWNAME_VARIANCE_SCALAR[] = "<c2>";	// name of view for variance of scalar
const char VIEWNAME_SKEWNESS_SCALAR[] = "<c3>";	// name of view for skewness of scalar
const char VIEWNAME_KURTOSIS_SCALAR[] = "<c4>";	// name of view for kurtosis of scalar
const char VIEWNAME_MEAN_TM[] = "<tm>";		// name of view for mean micromixing timescale

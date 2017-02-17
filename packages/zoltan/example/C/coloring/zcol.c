/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */
/***************************************************************
  File   : zcoldrive.c
  Date   :
  Author : Umit V. Catalyurek
  Descr  : based on Zoltan simpleGraph example, a driver
           for coloring codes that generates 5, 7 or 9
           point stencil on the fly and colors them.
 ***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "zoltan.h"

/*
#define _DEBUG
#define ZOLTANV31
*/

#define pID(i, j)  ((i)*uz->procC + (j))
#define gID(i, j)  ((i)*uz->meshC + (j))
#define lID(i, j)  (((i)-uz->_sr)*(uz->_ec-uz->_sc) + ((j)-uz->_sc))

#define getR(id)   ((id) / uz->meshC)
#define getC(id)   ((id) % uz->meshC)

#define pIDfRC(r, c)    pID(r/(uz->meshR/uz->procR), c/(uz->meshC/uz->procC))
#define pIDfGID(gid)    pIDfRC(getR(gid), getC(gid))

#define lIDfGID(gid)    lID(getR(gid), getC(gid))
#define gIDfLID(lid)    gID(uz->_sr+(lid/(uz->_ec-uz->_sc)), uz->_sc+(lid%(uz->_ec-uz->_sc)))

static int *xadj=NULL, *adj=NULL, *adjTemp=NULL;

struct edge {
    int c, r;
};

struct edge neig[8] = {
    {0, -1}, {1, 0}, {0, 1}, {-1, 0},
    {-1, -1}, {1, 1}, {1, -1}, {-1, 1}
};

typedef struct {
    int        stencil; /* for 5, 7 and 9-point, 4, 6 and 8 edges will be created, respectively, 
                           hence stencil value will be 4, 6 or 8 */

    int        procR, procC; /* proc mesh is R x C  */
    int        meshR, meshC; /* input mesh is R x C */
    int        _sr, _er, _sc, _ec; /* proc's start (inclusive) and end (exclusive) rows/cols of input mesh */
    int        myRank, myR, myC;
    int        numProcs;
    double     beta;
    int        numredge;
    int        *redgeto; /* if beta>0; contains end point of my random edges;
                            redgeto[lID(i,j)] = -1 means no random edge,
                                                    >= 0 points to end point */    
} UZData;

extern unsigned int Zoltan_Rand_InRange(unsigned int *, unsigned int);

double u_wseconds(void)
{
 struct timeval tp;

 gettimeofday(&tp, NULL);
 
 return (double) tp.tv_sec + (double) tp.tv_usec / 1000000.0;
}


static int get_number_of_objects(void *data, int *ierr)
{
    UZData *uz = (UZData *) data;
    
    *ierr = ZOLTAN_OK;
    
    return (uz->_er-uz->_sr) * (uz->_ec-uz->_sc);
}

static void get_object_list(void *data, int sizeGID, int sizeLID,
                            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                            int wgt_dim, float *obj_wgts, int *ierr)
{
    UZData *uz = (UZData *) data;
    int i, j, next=0;

  /* By setting parameters, I previously gave the Zoltan library 
   * the length of my global IDs, local IDs, and object weights.
   */

    if (sizeGID != 1 || sizeLID !=1 ){  /* My global/local IDs are 1 integer */
        *ierr = ZOLTAN_FATAL;
        return;
    }

    for (i=uz->_sr; i<uz->_er; ++i)
        for (j=uz->_sc; j<uz->_ec; ++j) {
            globalID[next] = gID(i, j);   /* application wide global ID */
            localID[next] = lID(i, j);    /* process specific local ID  */

            /* sanity check */
            if (lID(i, j)!=next) {
                fprintf(stderr, "umit something is wrong lID(%d, %d) is %d but next is %d\n", i, j, lID(i, j), next);
                exit(1);
            }
                
#if 0            
            if (wgt_dim>0)
                obj_wgts[next] = uz->stencil;  /* weight */
#endif
            ++next;
        }

  *ierr = ZOLTAN_OK;
}


/* let's keep this in the code for now; may be we'll use it later */
#if 0
/* 
 **************************************************************
 * Prototype: ZOLTAN_NUM_GEOM_FN
 * Return the dimension of a vertex, for geometric methods
 **************************************************************
 */
static int get_num_geometry(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return 2;
}

/* 
 **************************************************************
 * Prototype: ZOLTAN_GEOM_MULTI_FN
 * Return the coordinates of my objects (vertices), for
 * geometric methods.
 **************************************************************
 */
static void get_geometry_list(void *data, int sizeGID, int sizeLID,
                              int num_obj,
                              ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                              int num_dim, double *geom_vec, int *ierr)
{
    UZData *uz = (UZData *) data;
    int i, j, next=0;
    
    if ( (sizeGID != 1) || (sizeLID != 1) || (num_dim != 2)){
        *ierr = ZOLTAN_FATAL;
        return;
    }

    for (i=uz->_sr; i<uz->_er; ++i)
        for (j=uz->_sc; j<uz->_ec; ++j) {
            geom_vec[next++] = (double) i;
            geom_vec[next++] = (double) j;
        }

    *ierr = ZOLTAN_OK;
}
#endif

static int filledges(int lid, int *adjID, int *adjProc, UZData *uz)
{
    int *nextID = adjID, *ptr;
    int *nextProc = adjProc;
    int gid = gIDfLID(lid);
    int r = getR(gid), r2;
    int c = getC(gid), c2, e, j;

    /* printf(" %d (%d,%d) : ", id, r, c); */
    for (e=0; e<uz->stencil; ++e) {
        r2 = r+neig[e].r;
        c2 = c+neig[e].c;
            
        if (r2>=0 && r2<uz->meshR &&
            c2>=0 && c2<uz->meshC) {
            *nextID++ = gID(r2, c2);
            if (adjProc)
                *nextProc++ = pID(r2/(uz->meshR/uz->procR), c2/(uz->meshC/uz->procC));
                /* printf(" %d (%d, %d) [%d] ", *(nextID-1), r2, c2, *(nextProc-1)); */
            }
            
    }
    if (uz->redgeto && (uz->redgeto[lid]>=0)) {
        for (ptr=adjID; ptr<nextID && *ptr!=uz->redgeto[lid]; ++ptr);
        if (ptr>=nextID) {
            r2 = getR(uz->redgeto[lid]);
            c2 = getC(uz->redgeto[lid]);
            
            *nextID++ = uz->redgeto[lid];
            if (adjProc)
                *nextProc++ = pID(r2/(uz->meshR/uz->procR), c2/(uz->meshC/uz->procC));
            /* printf(" %d (%d, %d) [%d] ", *(nextID-1), r2, c2, *(nextProc-1)); */
        }
    }

    if (xadj) {
        for (j=xadj[lid]; j<xadj[lid+1]; ++j) {
            for (ptr=adjID; ptr<nextID && *ptr!=uz->redgeto[lid]; ++ptr);
            if (ptr>=nextID) {
                *nextID++ = adj[j];
                if (adjProc)
                    *nextProc++ = pIDfGID(adj[j]);
            }
        }
    }

    /* printf("\n"); */
    return nextID-adjID;
}
    

static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
                               int num_obj,
                               ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                               int *numEdges, int *ierr)
{
    UZData *uz = (UZData *) data;
    int i, lid;
    
    if ( (sizeGID != 1) || (sizeLID != 1)){
        *ierr = ZOLTAN_FATAL;
        return;
    }
    
    for (i=0;  i < num_obj ; ++i){
        lid = lIDfGID(globalID[i]);
        numEdges[i] = filledges(lid, adjTemp, NULL, uz);
    }
    
    *ierr = ZOLTAN_OK;
}

static void get_edge_list(void *data, int sizeGID, int sizeLID,
        int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int *num_edges, 
        ZOLTAN_ID_PTR nborGID, int *nborProc,
        int wgt_dim, float *ewgts, int *ierr)
{
    UZData *uz = (UZData *) data;
    int i;
    ZOLTAN_ID_PTR nextID; 
    int *nextProc;

    if ( (sizeGID != 1) || (sizeLID != 1) ||
         (wgt_dim != 0)) {      /* we are not using edge weights */
        *ierr = ZOLTAN_FATAL;
        return;
    }
    
    nextID = nborGID;
    nextProc = nborProc;

    /* printf("querying %d vertices\n", num_obj); */
    for (i=0;  i < num_obj ; ++i) {
        int lid = lIDfGID(globalID[i]), ne;

        ne = filledges(lid, (int *) nextID, nextProc, uz);
        nextID += ne;
        nextProc += ne;
    }
    *ierr = ZOLTAN_OK;
}

void computeProcMesh(UZData *uz)
{
    int i = (int) sqrt((double)uz->numProcs+0.1);
    while (uz->numProcs % i)
        --i;
    uz->procR = i;
    uz->procC = uz->numProcs / i;
}

void printStats(char *st, int val,  int myRank, int nProc)
{
    int *allval=NULL, minv, maxv, totv, i;

    if (!myRank)
        allval = (int *) malloc(sizeof(int)*nProc);

    MPI_Gather(&val, 1, MPI_INT, allval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!myRank) {
        minv = maxv = totv = allval[0];
        for (i=1; i<nProc; ++i) {
            if (allval[i]>maxv)
                maxv = allval[i];
            if (allval[i]<minv)
                minv = allval[i];
            totv += allval[i];
        }

        printf("%-15s  Min: %10d  Max: %10d  Avg: %10.1lf  Tot: %10d\n", st, minv, maxv, (double) totv/(double)nProc, totv);
        free(allval);
    }
}

#ifdef _DEBUG
void saveColor(char *exename, UZData *uz, int *gid_list, int *color, int ngids)
{
    char fname[2048];
    FILE *fp;
    int i;

    sprintf(fname, "%s.%dx%d.%d-point.%d.Of.%dx%d.txt", exename,
            uz->meshR, uz->meshC, uz->stencil+1, uz->myRank, uz->procR, uz->procC);

    if (!uz->myRank)
        printf("Saving the colors in file '%s'\n", fname);
    fp = fopen(fname, "w");
    for (i=0; i<ngids; ++i) 
        fprintf(fp, "%d %d\n", gid_list[i], color[i]);
    fclose(fp);
            
}
#endif


int main(int argc, char *argv[])
{
    int rc, i, ngids, maxcol, ncolors;
    float ver;
    struct Zoltan_Struct *zz=NULL;
#ifdef ZOLTANV31
    int numGidEntries, numLidEntries;
#else
    ZOLTAN_GRAPH_EVAL graph;
#endif
    int *color;
    ZOLTAN_ID_PTR gid_list;
    UZData guz, *uz=&guz;
    int msg_tag = 9999;
    int nlvtx, next, maxdeg=0;
    double times[9]={0.,0.,0.,0.,0.,0.,0.,0.}; /* Used for timing measurements */
    double gtimes[9]={0.,0.,0.,0.,0.,0.,0.,0.}; /* Used for timing measurements */
    
    /******************************************************************
     ** Initialize MPI and Zoltan
     ******************************************************************/

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &uz->myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &uz->numProcs);

    MPI_Barrier(MPI_COMM_WORLD);
    times[0] = u_wseconds();
    
    rc = Zoltan_Initialize(argc, argv, &ver);
    if (rc != ZOLTAN_OK){
        fprintf(stderr, "Sorry Zoltan initialize failed...\n");
        goto End;
    }
    zz = Zoltan_Create(MPI_COMM_WORLD);

    if (argc<3 && !uz->myRank) {
        fprintf(stderr, "usage: %s [meshR] [meshC] [X-point stencil] [procR] [procC] [ws-beta] [<ZoltanParam>=<Val>] ...\n\n", argv[0]);
        fprintf(stderr, "ws-beta: is the probablity of adding an edge to a vertex to generate Watts-Strogatz graphs\n");
        fprintf(stderr, "Valid values for Stencil are 5, 7 and 9\n");
        fprintf(stderr, "Zoltan Coloring Parameters and values are\n");
        fprintf(stderr, "\tDISTANCE        : 1 or 2\n");
        fprintf(stderr, "\tSUPERSTEP_SIZE  : suggested >= 100\n"); 
        fprintf(stderr, "\tCOMM_PATTERN    : S or A\n");
        fprintf(stderr, "\tCOLOR_ORDER     : I, B, U\n");
        fprintf(stderr, "\tCOLORING_METHOD : F (for now)\n");
        fprintf(stderr, "\n");
    }

    uz->procR = uz->procC = 0;
    uz->meshR = uz->meshC = 1024;
    uz->stencil = 9;

    if (argc>1)
        uz->meshR = atoi(argv[1]);
    if (argc>2)
        uz->meshC = atoi(argv[2]);
    if (argc>3)
        uz->stencil = atoi(argv[3]);
    if (uz->stencil!=5 && uz->stencil!=7 && uz->stencil!=9) {
        fprintf(stderr, "\t invalid stencil value. Valid values are 5, 7 and 9. Assumed 9.\n");
        uz->stencil = 9;
    }
    --uz->stencil;

    if (argc>4)
        uz->procR = atoi(argv[4]);
    if (argc>5)
        uz->procC = atoi(argv[5]);
    if (uz->procR <= 0 || uz->procC <= 0)
        computeProcMesh(uz);
    
    if (uz->procR*uz->procC!=uz->numProcs) {
        fprintf(stderr, "#Procs=%d but requested %dx%d Proc Mesh Partitioning...\n", uz->numProcs, uz->procR, uz->procC);
        goto End;
    }

    if (argc>6)
        uz->beta = atof(argv[6]);
    else
        uz->beta = 0.0;
    
    /* compute which part of mesh I will compute */
    uz->myR = uz->myRank / uz->procC;
    uz->myC = uz->myRank % uz->procC;

    uz->_sr = uz->myR * (uz->meshR / uz->procR);
    uz->_er = (uz->myR+1) * (uz->meshR / uz->procR);
    if (uz->_er>uz->meshR)
        uz->_er = uz->meshR;
    uz->_sc = uz->myC * (uz->meshC / uz->procC);
    uz->_ec = (uz->myC+1) * (uz->meshC / uz->procC);
    if (uz->_ec>uz->meshC)
        uz->_ec = uz->meshC;


    if ( (uz->meshR % uz->procR) !=0 || (uz->meshC % uz->procC)!=0) {
        printf("Mesh dimensions are not divisible with proc mesh.\nRequested mesh is %dx%d and proc mesh is %d x %d\n", uz->meshR, uz->meshC, uz->procR, uz->procC);
        exit(1);
    }
    nlvtx= (uz->_er-uz->_sr) * (uz->_ec-uz->_sc);

    if (uz->myRank==0)
        printf("Running %s on %d x %d processor mesh, generating %d-point %d x %d mesh with beta=%.3lf\n", argv[0], uz->procR, uz->procC, uz->stencil+1, uz->meshR, uz->meshC, uz->beta);

    times[1] = u_wseconds();    
    uz->numredge = 0;
    uz->redgeto = NULL;
    if (uz->beta>0) { /* create random edges for WS graph */
        int ngvtx=uz->meshC*uz->meshR, trsh=(int) (uz->beta*100.0);
        int ierr=0;
        int *edges=NULL, *redges=NULL, *proclist=NULL, nedge;
        ZOLTAN_COMM_OBJ *plan;
            
        uz->redgeto = (int *) malloc(nlvtx*sizeof(int));
        for (i=0; i<nlvtx; ++i) {
            int rv = Zoltan_Rand_InRange(NULL, 100);
            if ( rv < trsh) {
                if ((uz->redgeto[i] = Zoltan_Rand_InRange(NULL,  ngvtx))==gIDfLID(i)) /* is it a self edge */
                    uz->redgeto[i] = -1;
                else 
                    ++uz->numredge;
            } else
                uz->redgeto[i] = -1;
        }

        edges = (int *) malloc(sizeof(int)*2*uz->numredge);
        proclist = (int *) malloc(sizeof(int)*uz->numredge);
        next = 0;
        for (i=0; i<nlvtx; ++i)
            if (uz->redgeto[i]>0) {
                edges[2*next] = uz->redgeto[i];
                edges[2*next+1] = gIDfLID(i);
                proclist[next] = pIDfGID(uz->redgeto[i]);
                ++next;
            }

        ierr = Zoltan_Comm_Create(&plan, uz->numredge, proclist, MPI_COMM_WORLD,
                                  msg_tag, &nedge);
        
        redges = (int *) malloc(sizeof(int)*2*nedge);
        
        --msg_tag;
        ierr |= Zoltan_Comm_Do(plan, msg_tag, (char *) edges, 2*sizeof(int),
                               (char *) redges);
        ierr |= Zoltan_Comm_Destroy(&plan);
        free(proclist);
        free(edges);
                
        if (ierr) {
            printf("error while communicating edges!\n");
            exit(1);
        }
        xadj = (int *) calloc(1+nlvtx, sizeof(int));
        adj = (int *) malloc(sizeof(int)*nedge);
        for (i=0; i<nedge; ++i)  {
            if (redges[2*i] < gID(uz->_sr, uz->_sc) || redges[2*i] >= gID(uz->_er, uz->_ec)) {
                printf("[%d/%d] received gid=%d doesn't blong to processor range [%d, %d) should go to proc %d\n", uz->myRank, uz->numProcs, redges[2*i],  gID(uz->_sr, uz->_sc), gID(uz->_er, uz->_ec), pIDfGID(redges[2*i]));
            }
                
            ++xadj[lIDfGID(redges[2*i])];
        }
        xadj[nlvtx] = nedge;
        maxdeg = xadj[0];
        for (i=1; i<nlvtx; ++i) {
            maxdeg = xadj[i]>maxdeg ? xadj[i] : maxdeg;
            xadj[i] += xadj[i-1];
        }
        
        for (i=0; i<nedge; ++i) {
            int u = lIDfGID(redges[2*i]);
            int v = redges[2*i+1];
            adj[--xadj[u]] = v;
        }
        free(redges);                    
    }
    maxdeg += uz->stencil+1;
    adjTemp = (int *) malloc(sizeof(int)*2*maxdeg);
    times[2] = u_wseconds();
    
    /*
      printf("My rank %d/%d at proc-mesh loc (%d, %d) generating [%d, %d) x [%d, %d) + %d random edges TotEdge=%d\n", uz->myRank, uz->numProcs, uz->myR, uz->myC, uz->_sr, uz->_er, uz->_sc, uz->_ec, uz->numredge, xadj[nlvtx]);  */
    printStats("Number of Vertices  ", nlvtx, uz->myRank, uz->numProcs);
    if (xadj)
        printStats("Number of Rand Edges", xadj[nlvtx], uz->myRank, uz->numProcs);

    
    /* General parameters */
#ifndef ZOLTANV31
    Zoltan_Set_Param(zz, "GRAPH_BUILD_TYPE", "FAST_NO_DUP");
#endif

    /* General parameters */
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "3");
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");


    /* coloring parameters */
    Zoltan_Set_Param(zz, "SUPERSTEP_SIZE", "500"); /* let's make S=500 default */
    for (i=7; i<argc; ++i) {
        char param[256], *eq;

        if (!uz->myRank)
            printf("processing argv[%d]='%s'\n", i, argv[i]);
        strncpy(param, argv[i], sizeof(param));
        eq = strchr(param, '=');
        if (!eq) {
            fprintf(stderr, "invalid argument '%s', Zoltan Paramters should be in the format <ZoltanParam>=<Val>\n", param);
            goto End;
        }
        *eq = 0;
        Zoltan_Set_Param(zz, param, eq+1);
    }


#if 0    
    /* Graph parameters */
    Zoltan_Set_Param(zz, "CHECK_GRAPH", "2");
#endif

    /* set call backs */
    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_objects, uz);
    Zoltan_Set_Obj_List_Fn(zz, get_object_list, uz);
    Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list, uz);
    Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list, uz);

#if 0    
#ifndef ZOLTANV31
    Zoltan_LB_Eval_Graph(zz, 0, &graph);

    if (!uz->myRank) {
        printf("EdgeCut   Min=%8.0f  Max=%8.0f  Sum=%8.0f\n", graph.cuts[EVAL_GLOBAL_MIN], graph.cuts[EVAL_GLOBAL_MAX], graph.cuts[EVAL_GLOBAL_SUM]);
        printf("#Vertices Min=%8.0f  Max=%8.0f  Sum=%8.0f imbal=%.2f\n", graph.nobj[EVAL_GLOBAL_MIN], graph.nobj[EVAL_GLOBAL_MAX], graph.nobj[EVAL_GLOBAL_SUM], graph.obj_imbalance);        
    }
#endif
#endif

    /* now color */
    ngids = get_number_of_objects(uz, &rc);

    gid_list = (ZOLTAN_ID_PTR) malloc(sizeof(ZOLTAN_ID_TYPE) * ngids);
#ifndef ZOLTANV31
    next = 0;
    for (i=uz->_sr; i<uz->_er; ++i) {
        int j;
        for (j=uz->_sc; j<uz->_ec; ++j) {
            gid_list[next++] = i*uz->meshC + j;
        }
    }    
#endif
    color = (int *) malloc(sizeof(int) * ngids);    

    MPI_Barrier(MPI_COMM_WORLD);
    times[3] = u_wseconds();    
#ifdef ZOLTANV31
    rc = Zoltan_Color(zz, /* input (all remaining fields are output) */
                      &numGidEntries,  /* Number of integers used for a global ID */
                      &numLidEntries,  /* Number of integers used for a local ID */
                      ngids,           /* #objects to color in this proc */
                      gid_list,        /* global ids of colored vertices */
                      NULL,            /* we ignore local ids */
                      color);          /* result color */    
#else    
    rc = Zoltan_Color(zz, /* input (all remaining fields are output) */
                      1,  /* Number of integers used for a global ID */
                      ngids,           /* #objects to color in this proc */
                      gid_list,        /* global ids of colored vertices */
                      color);          /* result color */
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    times[4] = u_wseconds();
    MPI_Reduce(times, gtimes, 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


    if (rc != ZOLTAN_OK) 
        fprintf(stderr, "Zoltan_Color failed with return code %d...\n", rc);

    for (maxcol=i=0; i<ngids; ++i)
        if (color[i] > maxcol)
            maxcol = color[i];
    MPI_Reduce(&maxcol, &ncolors, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    if (uz->myRank==0) {
        struct rusage usage;
                
        printf("%s setup             Proc-0: %8.2lf   Max: %8.2lf\n", argv[0], times[1]-times[0], gtimes[1]-gtimes[0]);
        printf("%s gen rand edges    Proc-0: %8.2lf   Max: %8.2lf\n", argv[0], times[2]-times[1], gtimes[2]-gtimes[1]);
        printf("%s set gids          Proc-0: %8.2lf   Max: %8.2lf\n", argv[0], times[3]-times[2], gtimes[3]-gtimes[2]);
        printf("%s Zoltan_Color call Proc-0: %8.2lf   Max: %8.2lf\n", argv[0], times[4]-times[3], gtimes[4]-gtimes[3]);
        printf("%s Coloring Time    : %.2lf   # Colors used : %d\n", argv[0], gtimes[4]-gtimes[0], ncolors);
        getrusage(RUSAGE_SELF, &usage);
        printf("%s maxrss=%ld minflt=%ld majflt=%ld nswap=%ld\n", argv[0], usage.ru_maxrss, usage.ru_minflt, usage.ru_majflt, usage.ru_nswap);                
    }
    
#ifdef _DEBUG
    saveColor(argv[0], uz, (int *) gid_list, color, ngids);
#endif

    /******************************************************************
     ** Clean up
     ******************************************************************/

    if (gid_list)
        free(gid_list);
    if (color)
        free(color);
    if (xadj)
        free(xadj);
    if (adj)
        free(adj);
    if (adjTemp)
        free(adjTemp);
    if (uz->redgeto)
        free(uz->redgeto);

End:    
    Zoltan_Destroy(&zz);
    MPI_Finalize();

    return 0;
}

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<unistd.h>
#include<time.h>
#include<metis.h>

#define NOS 6

/* data structure for complete cubed sphere, CS */
typedef struct {
    int dim;
    int nParts;
   int nVerts;         // number of verts, non-duplicate
   int nEls;           // number of elements
   double **my_coords; // vert coordinates
   int **my_nadj;      // vert adjacency
   unsigned nvattribs; // number of vert attributes per node
   int **vattribs;     // vert attributes
   int **my_e_n;       // element to vert table
   int *epart;         // element partitioning
   int *vpart;         // vert partitioning
   int *activeList;    // element to vert table
} CS;

/* data structure for individual sixths, that will make the cubed sphere */
typedef struct {
   int dim;
   unsigned sizes[3];
   unsigned elSizes[3];
   unsigned basis[3];

   // the angle limits of each grid
   double eta_limits[2];
   double xi_limits[2];
   double r_limits[2];  // radial min & max

   unsigned nVerts;
   unsigned nEls;
   unsigned nvattribs; // number of vert attributes per node
   int **vattribs;     // vertice attributes
   int *activeList;    // the sixth's element-node mapping - indexed with local elements resulting in global vert ids */
   int **my_e_n;       // the sixth's element-node mapping - indexed with local elements resulting in global vert ids
   int **my_nadj;      // the sixth's node-node mapping - indexed with local verts resulting in global vert ids
   double **my_coords; // the sixth's node xyz coords
   unsigned v_offset;  // global vert ordering offset
   unsigned e_offset;  // global element ordering offset
} Sixth;

/* ID_TYPE used to refer to local sixth ordering or global CS ordering of elements and verts */
typedef enum {
   GLOBAL,
   LOCAL
} ID_TYPE;

#define VA_BC 0 // variable attribute boundary condition

/* function definitions */
void partitionCS( CS *mesh );
void write_vtu( CS *mesh, char* filename ) ;
void write_ascii( CS *cs );
void write_sixth_vtu( Sixth* self, char* name );

void build_rotation_matrix( double rot_around_x, double rot_around_y, double *mat )
/*@
   builds a 3D rotation matrix. Used to rotate a build sixth into a new position
   input params - rot_around_x, rot_around_y - are in radians
@*/
{

   mat[0] = cos(rot_around_y);
   mat[1] = 0;
   mat[2] = -sin(rot_around_y);

   mat[3] = 0;
   mat[4] = cos(rot_around_x);
   mat[5] = -sin(rot_around_x);

   mat[6] = sin(rot_around_y);
   mat[7] = sin(rot_around_x);
   mat[8] = cos(rot_around_x)*cos(rot_around_y);
}


void Sixth_Setup(  Sixth* self, unsigned e_offset, unsigned v_offset, unsigned* inputSizes, double *xi_limits, double *eta_limits, double *r_limits )
{
   /*@
      Setup the Sixth data structure.
      Assume inputSizes is length 3
      @*/

   int **e_n=NULL;
   int **nadj=NULL;
   int *activeList=NULL;
   int **vattribs=NULL;
   double **coords=NULL;
   unsigned nEls, nVerts, v_i;
   int dim = self->dim;

   memcpy( self->sizes, inputSizes, 3*sizeof(unsigned) );

   // setup basis
   self->basis[0] = 1;
   self->basis[1] = self->sizes[0];
   if( dim == 3 ) {
       self->basis[2] = self->sizes[0]*self->sizes[1];
   } else {
       self->basis[2] = NAN;
   }

   // setup elment sizes
   self->elSizes[0] = self->sizes[0]-1;
   self->elSizes[1] = self->sizes[1]-1;
   if( dim == 3 ) {
       self->elSizes[2] = self->sizes[2]-1;
   } else {
       self->elSizes[2] = NAN;
   }

   // used to create globalIds
   self->v_offset = v_offset;
   self->e_offset = e_offset;

   self->nvattribs=1;

   if( dim == 3 ) {
       self->nEls = self->elSizes[0] * self->elSizes[1] * self->elSizes[2];
       self->nVerts = self->sizes[0]*self->sizes[1]*self->sizes[2];
   } else {
       self->nEls = self->elSizes[0] * self->elSizes[1];
       self->nVerts = self->sizes[0]*self->sizes[1];
   }

   nEls = self->nEls;
   nVerts = self->nVerts;

   memcpy( self->xi_limits, xi_limits, 2*sizeof(double) );
   memcpy( self->eta_limits, eta_limits, 2*sizeof(double) );
   memcpy( self->r_limits, r_limits, 2*sizeof(double) );

   /*
      build memory for local
       -node-node adj structure( nadj )
       -element-node structure( e_n )
       -node xyz structure( coords )
       -a list of active verts, initially all are active (activeList)
    */
   vattribs = malloc( nVerts * sizeof(int*) );
   activeList = malloc( nVerts *sizeof(int));
   nadj = malloc( nVerts *sizeof(int*));
   e_n = malloc( nEls *sizeof(int*));
   coords = malloc( nVerts *sizeof(double*));

   // build node adjacency and element-node tables, leave them as the 3D length
   for( v_i=0; v_i<nVerts; v_i++ ) {
      activeList[v_i]=1;
      nadj[v_i] = malloc( 6*sizeof(int) );

      // allocate & initialise attributes
      vattribs[v_i] = malloc( self->nvattribs*sizeof(int) );
      vattribs[v_i][VA_BC] = 0;

      memset( nadj[v_i], -1, 6*sizeof(int) ); // set to -1 to represent no nbr connection

      coords[v_i] = malloc( 3*sizeof(double) );
      memset( coords[v_i], 0, 3*sizeof(double) );
   }
   for( v_i=0; v_i<nEls; v_i++) {
      e_n[v_i] = malloc( 8*sizeof(int) );
      memset( e_n[v_i], -1, 8*sizeof(int) ); // set to -1, i.e. uninitialised
   }

   self->my_nadj = nadj;
   self->my_e_n = e_n;
   self->my_coords = coords;
   self->activeList = activeList;
   self->vattribs = vattribs;
}

void Sixth_FreeMem( Sixth *self )
{
   /*@
     A function to deallocate a sixth

     @*/

   int v_i, e_i;

   // free local node-node array
   for( v_i = 0 ; v_i<self->nVerts; v_i++ ) {
      free( self->my_nadj[v_i] );
      free( self->my_coords[v_i] );
      free( self->vattribs[v_i] );
   }
   free( self->my_nadj );
   self->my_nadj=NULL;
   free( self->vattribs );
   self->vattribs=NULL;
   free( self->my_coords );
   self->my_coords=NULL;
   free( self->activeList );
   self->activeList=NULL;

   // free local element-node array
   for( e_i = 0; e_i<self->nEls; e_i++  ) {
      free( self->my_e_n[e_i] );
   }
   free( self->my_e_n );
   self->my_e_n=NULL;
}

unsigned Sixth_eProject( Sixth *self, unsigned *ijk, ID_TYPE id_tag )
{
   /*@
     Given ijk parametrisation find the spherical id
     @*/

   unsigned elId;

   // sanity check input
   assert( ijk[0] < self->elSizes[0] );
   assert( ijk[1] < self->elSizes[1] );
   assert( ijk[2] < self->elSizes[2] );

   // then project
   elId = ijk[0] +
          ijk[1] * self->elSizes[0] +
          ijk[2] * (self->elSizes[0] * self->elSizes[1] );

   if( id_tag == GLOBAL )
      elId += self->e_offset;

   return elId;
}


unsigned Sixth_Project( Sixth *self, unsigned *ijk, ID_TYPE id_tag )
{
   /*@
     Given ijk parametrisation find the "global" id
     @*/

   unsigned vert_id;

   // sanity check input
   assert( ijk[0] < self->sizes[0] );
   assert( ijk[1] < self->sizes[1] );
   if ( self->dim == 3 ) {
       assert( ijk[2] < self->sizes[2] );
   } else {
       ijk[2] = 0;
   }

   // then project
   vert_id = ijk[0] +
             ijk[1] * self->sizes[0] +
             ijk[2] * (self->sizes[0] * self->sizes[1] );

   if( id_tag == GLOBAL )
      vert_id += self->v_offset;

   return vert_id;
}

void Sixth_Lift( Sixth *self, unsigned id, unsigned *ijk )
{
   /*@
     Given the id, find the ijk parametrisation for sixth
     @*/
   unsigned rem = id;
   unsigned d_i;
   int dim = self->dim;

   for( d_i = dim; d_i > 0; d_i-- ) {
      unsigned	dimInd = d_i - 1;
      div_t		divRes;

      divRes = div( rem, self->basis[dimInd] );
      ijk[dimInd] = divRes.quot;
      rem = divRes.rem;
   }
}

// nbr node indices -- ordering represents the nbr direction
enum NBR_TYPE {
   NBR_W = 0,
   NBR_N = 1,
   NBR_E = 2,
   NBR_S = 3,
   NBR_U = 4,
   NBR_D = 5
};

void build_sixth_e_n( Sixth *self )
{
   /*@
     function build a complete sixth element to node mapping
     @*/

   unsigned kk,jj,ii, n_depth_els, n_surface_els;
   int elId, lvid, gvid;
   int **e_n = self->my_e_n;
   int **nbrs = self->my_nadj;
   int dim = self->dim;

   // for each element there are 8 nodes

   /*
   vtu order for hexahedra i.e VTK_HEXAHEDRA = 12

     7 __________ 6
    / :         / \
   4 ________ 5    \
   |  :       |     \
   |  3 ------\-----2
   | /         \   /
    /           \ /
   0 ___________ 1


   */

   if( dim == 3 ) {
       n_surface_els=self->elSizes[0];
       n_depth_els=self->elSizes[2];

       elId = 0;
       gvid=self->v_offset;
       for( kk=0; kk<n_depth_els; kk++ ) {
          for( jj=0; jj<n_surface_els; jj++ ) {
             for( ii=0; ii<n_surface_els; ii++ ) {
                // get starting node for element
                e_n[elId][0] = gvid;
                // calc index for local vert array, i.e. lvid
                lvid = gvid - self->v_offset;

                // get neighbouring nodes in vtu order
                gvid = e_n[elId][1] = nbrs[lvid][NBR_E];
                lvid = gvid - self->v_offset;

                gvid = e_n[elId][2] = nbrs[lvid][NBR_N];
                lvid = gvid - self->v_offset;

                gvid = e_n[elId][3] = nbrs[lvid][NBR_W];
                lvid = gvid - self->v_offset;

                // going to node 7 now
                gvid = e_n[elId][7] = nbrs[lvid][NBR_U];
                lvid = gvid - self->v_offset;

                gvid = e_n[elId][4] = nbrs[lvid][NBR_S];
                lvid = gvid - self->v_offset;

                gvid = e_n[elId][5] = nbrs[lvid][NBR_E];
                lvid = gvid - self->v_offset;

                e_n[elId][6] = nbrs[lvid][NBR_N];

                gvid = e_n[elId][1]; //start on lvid 1 next
                elId++;
             }
             gvid = e_n[elId-n_surface_els][3]; //start in south-west corner again one node more north
          }
          gvid = e_n[elId-(n_surface_els*n_surface_els)][4]; //start on in south-west corner one node up
       }
   } else {
       n_surface_els=self->elSizes[0];
       n_depth_els=self->elSizes[1];

       elId = 0;
       gvid=self->v_offset;
          for( jj=0; jj<n_depth_els; jj++ ) {
             for( ii=0; ii<n_surface_els; ii++ ) {
                // get starting node for element
                e_n[elId][0] = gvid;
                // calc index for local vert array, i.e. lvid
                lvid = gvid - self->v_offset;

                // get neighbouring nodes in vtu order
                gvid = e_n[elId][1] = nbrs[lvid][NBR_E];
                lvid = gvid - self->v_offset;

                gvid = e_n[elId][2] = nbrs[lvid][NBR_N];
                lvid = gvid - self->v_offset;

                gvid = e_n[elId][3] = nbrs[lvid][NBR_W];
                lvid = gvid - self->v_offset;

                gvid = e_n[elId][1]; //start on lvid 1 next
                elId++;
             }
             gvid = e_n[elId-n_surface_els][3]; //start in south-west corner again one node more north
          }
   }

}

int Sixth_OwnsVert( Sixth *g, int vid )
{
   /*@
     Tests if a global vert id belongs to the given Sixth g
    @*/
   int val=1;
   if ( vid < g->v_offset ) val=0;
   if ( vid > (g->v_offset + g->nVerts) ) val=0;

   return val;
}

void addSixthToGiant( CS* giant, Sixth *g1 )
{
   /*@
     this function adds the adjacency list and element-node mapping into a global indices
     @*/

   int **e_n = g1->my_e_n;
   int **nbrs = g1->my_nadj;
   int **vattribs = g1->vattribs;
   double **coords = g1->my_coords;
   int v_i, e_i, index;

   int *v_count = &giant->nVerts;
   int *el_count = &giant->nEls;

   int **g_nadj = giant->my_nadj;
   int **g_e_n = giant->my_e_n;
   double **g_coords = giant->my_coords;
   int **g_vattribs = giant->vattribs;

   for( v_i=0; v_i<g1->nVerts; v_i++ ) {
      giant->activeList[(*v_count)] = g1->activeList[v_i];

      /* add nbrs list into g_nadj */
      memcpy(g_nadj[(*v_count)], nbrs[v_i], 6*sizeof(int) );

      /* add vert attributes into giant */
      memcpy(g_vattribs[(*v_count)], vattribs[v_i], giant->nvattribs*sizeof(int) );

      /* add node geometry to giant */
      memcpy(g_coords[(*v_count)], coords[v_i], 3*sizeof(double) );

      (*v_count)++;
   }
   for( e_i=0; e_i<g1->nEls ; e_i++) {
      memcpy( g_e_n[*el_count], g1->my_e_n[e_i], 8*sizeof(int) );
      (*el_count)++;
   }
}

void evaluate_sixth_geom( Sixth *self, double x_rot, double y_rot )
{
   /*@
     Builds the geometry of the verts for sixth self
     @*/

   unsigned ii,jj,kk, v_i;
   double d_xi, d_eta, d_r;
   double X,Y,r,xi,eta,d,R,a,b,c;
   double rot[9];
   double **coords=self->my_coords;
   int dim = self->dim;

   // build sixth's rotation matrix
   build_rotation_matrix( x_rot, y_rot, rot );
   v_i=0;

   // denominator is (size[]-1) because we are calculating the spacing
   d_xi = fabs( (self->xi_limits[1]-self->xi_limits[0])/(double)(self->sizes[0]-1) );
   d_eta = fabs( (self->eta_limits[1]-self->eta_limits[0])/(double)(self->sizes[1]-1) );
   d_r = fabs( (self->r_limits[1]-self->r_limits[0])/(double)(self->sizes[2]-1) );


   if( dim == 3 ) {
       for( kk=0; kk<self->sizes[2]; kk++ ) {
          r = self->r_limits[0] + kk*d_r;

          for( jj=0; jj<self->sizes[1]; jj++ ) {
             eta = self->eta_limits[0] + jj*d_eta;

             for( ii=0; ii<self->sizes[0]; ii++ ) {
                xi = self->xi_limits[0] + ii*d_xi;

                // calculate points in X-Y plane
                X = 1*tan( xi );
                Y = 1*tan( eta );

                // calc projection vars.
                d = sqrt( 1 + X*X + Y*Y );
                R = sqrt(3)*r;

                // project points onto spherical surface of radius R
                /*
                a = r/d * X;
                b = r/d * Y;
                c = r/d * 1;
                */
                a =  X;
                b =  Y;
                c = r * 1;

                coords[v_i][0] = rot[0]*a + rot[1]*b + rot[2]*c;
                coords[v_i][1] = rot[3]*a + rot[4]*b + rot[5]*c;
                coords[v_i][2] = rot[6]*a + rot[7]*b + rot[8]*c;

                v_i++;
             }
          }
       }
   } else {
    for( jj=0; jj<self->sizes[1]; jj++ ) {
         eta = self->eta_limits[0] + jj*d_eta;

         for( ii=0; ii<self->sizes[0]; ii++ ) {
            xi = self->xi_limits[0] + ii*d_xi;

            // calculate points in X-Y plane
            X = 1*tan( xi );
            Y = 1*tan( eta );

            // calc projection vars.
            d = sqrt( 1 + X*X + Y*Y );

            // project points onto spherical surface of radius R
            /*
            a = r/d * X;
            b = r/d * Y;
            c = r/d * 1;
            */
            a =  X;
            b =  Y;
            c = 0;

            coords[v_i][0] = rot[0]*a + rot[1]*b + rot[2]*c;
            coords[v_i][1] = rot[3]*a + rot[4]*b + rot[5]*c;
            coords[v_i][2] = rot[6]*a + rot[7]*b + rot[8]*c;

            v_i++;
         }
      }
   }
}

void build_sixth_nadj( Sixth *sixth )
{
   /*@

    builds the simple node connectivity for each sixth
    and saves information in nbrs

    nbrs order: nbrs[vert_id] = [ west, north, east, south, down, up ]

    NOTE: globalId represents a "global" nodeID
    @*/
   unsigned  s_verts, s_i, ijk[3], globalId, v_i;
   int **nbrs = sixth->my_nadj;
   int **vattribs = sixth->vattribs;
   int dim = sixth->dim;
   s_verts = sixth->nVerts;

   for( v_i=0; v_i<s_verts; v_i++ ) {
      // get local (sixth) ijk
      Sixth_Lift( sixth, v_i, ijk );

      // get west neighbour
      if( ijk[0]>0 ) {
         ijk[0]--;
         globalId = Sixth_Project( sixth, ijk, GLOBAL );
         nbrs[v_i][NBR_W] = globalId;
         ijk[0]++; //reset
      } else {
         nbrs[v_i][NBR_W] = -1;
      }

      // get north neighbour
      if( ijk[1]<sixth->sizes[1]-1 ) {
         ijk[1]++;
         globalId = Sixth_Project( sixth, ijk, GLOBAL );
         nbrs[v_i][NBR_N] = globalId;
         ijk[1]--; //reset
      } else {
         nbrs[v_i][NBR_N] = -1;
      }

      // get east neighbour
      if( ijk[0]<sixth->sizes[0]-1 ) {
         ijk[0]++;
         globalId = Sixth_Project( sixth, ijk, GLOBAL );
         nbrs[v_i][NBR_E] = globalId;
         ijk[0]--; //reset
      } else {
         nbrs[v_i][NBR_E] = -1;
      }

      // get south neighbour
      if( ijk[1]>0 ) {
         ijk[1]--;
         globalId = Sixth_Project( sixth, ijk, GLOBAL );
         nbrs[v_i][NBR_S] = globalId;
         ijk[1]++; //reset
      } else {
         nbrs[v_i][NBR_S] = -1;
      }

      if( dim == 3 ) {
          // get core neighbour
          if( ijk[2]>0 ) {
             ijk[2]--;
             globalId = Sixth_Project( sixth, ijk, GLOBAL );
             nbrs[v_i][NBR_D] = globalId;
             ijk[2]++; //reset
          } else {
             nbrs[v_i][NBR_D] = -1;
             vattribs[v_i][VA_BC] = 1; // 1 for inner radius node
          }

          // get surface neighbour
          if( ijk[2]<sixth->sizes[2]-1 ) {
             ijk[2]++;
             globalId = Sixth_Project( sixth, ijk, GLOBAL );
             nbrs[v_i][NBR_U] = globalId;
             ijk[2]--; //reset
          } else {
             nbrs[v_i][NBR_U] = -1;
             vattribs[v_i][VA_BC] = 2; // 1 for outer radius node
          }
      }
   }


}

int main(int argc, char** argv )
{
   unsigned n_surface_verts, n_depth_verts;
   unsigned n_surface_els, n_depth_els, dim, nParts;
   double outer_radius, inner_radius;
   Sixth sixths[6];
   CS giant;
   unsigned sixthSize[3], ijk[3], tmp;
   int nVerts, nEls, s_i, v_i, n_id, e_i;
   int **nbrs;
   double **verts;
   double r_limits[2], eta_limits[2], xi_limits[2];

   int **e_n, elsBuilt;

   double **global_coord;
   int nGlobalEls;
   int nGlobalVerts;
   int c;
   extern char *optarg;

   clock_t t2, t1 = clock();


   /***
     get input:
     geometry and number of elements
   ****/

   n_surface_els=-1;
   n_depth_els=-1;
   nParts=-1;
   while ( (c = getopt(argc, argv, "l:d:z:p:")) != -1) {
      switch( c ) {
      case 'l':
         n_surface_els=atoi(optarg);
         break;
      case 'p':
         nParts = atoi(optarg);
         break;
      case 'z':
         dim = atoi(optarg);
         break;
      case 'd':
         n_depth_els=atoi(optarg);
         break;
      case '?':
         printf("execute with\n ./go -l 11 -d 6 -z [2|3] -p 2\n");
         printf("where l is the length resolution and d is the depth resolution\n");
         exit(1);
      default:
         printf("execute with\n ./go -l 11 -d 6 -z [2|3] -p 2\n");
         printf("where l is the length resolution and d is the depth resolution\n");
         exit(1);
      }

   }
   if( n_surface_els==-1 || n_depth_els==-1 || dim ==-1 || nParts==-1) {
      printf("execute with\n ./go -l 11 -d 6 -z [2|3]\n");
      printf("where l is the number of elements in length (d depth) per sixth\n");
      exit(1);
   }

   giant.nParts=nParts;
   giant.dim=sixths[0].dim=dim;

   /***
     get e_n, epart, n_n, npart
   ***/
   /* set initial sixth definition */
   n_surface_verts = n_surface_els+1;
   n_depth_verts = n_depth_els+1;
   outer_radius=6;
   inner_radius=3;

   /* if 3D or 2D */
   if( dim == 3 ) {
       nGlobalEls = n_surface_els*n_surface_els*n_depth_els;
       nGlobalVerts = n_surface_verts*n_surface_verts*n_depth_verts;
   } else {
       nGlobalEls = n_surface_els*n_depth_els;
       nGlobalVerts = n_surface_verts*n_depth_verts;
   }

   // initialise counter as they will be used to count the overall number of verts and elements as the sixths get stitched together
   giant.nVerts=0;
   giant.nEls=0;
   giant.nvattribs=1; // BC attribute

   giant.my_e_n = malloc( nGlobalEls * sizeof(int*) );
   for( e_i=0; e_i<nGlobalEls; e_i++ ) {
      giant.my_e_n[e_i] = malloc( 8*sizeof(int) );
   }

   giant.activeList = malloc( nGlobalVerts * sizeof(int) );
   giant.my_nadj = malloc( nGlobalVerts * sizeof(int*) );
   giant.vattribs = malloc( nGlobalVerts * sizeof(int*) );
   giant.my_coords = malloc( nGlobalVerts * sizeof(double*) );
   for( v_i=0; v_i<nGlobalVerts; v_i++ ) {
      giant.my_nadj[v_i] = malloc( 6*sizeof(int) );
      giant.vattribs[v_i] = malloc( giant.nvattribs*sizeof(int) );
      giant.my_coords[v_i] = malloc( 6*sizeof(double) );
   }

   r_limits[0] = inner_radius;
   r_limits[1] = outer_radius;
   sixthSize[0] = n_surface_verts;
   sixthSize[1] = n_depth_verts;
   sixthSize[2] = n_surface_verts;
   xi_limits[0]  = -M_PI/4;
   xi_limits[1]  = M_PI/4;
   eta_limits[0] = -M_PI/4;
   eta_limits[1] = M_PI/4;

   Sixth_Setup( &(sixths[0]), 0, 0, sixthSize, xi_limits, eta_limits, r_limits );
   build_sixth_nadj( &sixths[0] );
   build_sixth_e_n( &sixths[0] );
   evaluate_sixth_geom( &sixths[0],0,0 );

   t2 = clock();
   printf("Time to initialise sixth %g sec\n", (double)(t2-t1)/CLOCKS_PER_SEC );

   // debug
   write_sixth_vtu( &sixths[0], "sixth-0.vtu" );
 
   /*
      run function to join walls
   */

   t1=clock();

   t2=clock();
   printf("Time to join sixths %g sec\n", (double)(t2-t1)/CLOCKS_PER_SEC );

   /* add all sixths together into a giant mesh */
   t1=clock();
   addSixthToGiant( &giant, &sixths[0] );
   t2=clock();
   printf("Time to add sixths to global  %g sec\n", (double)(t2-t1)/CLOCKS_PER_SEC );

   giant.epart = malloc( giant.nEls*sizeof(int) );
   giant.vpart = malloc( giant.nVerts*sizeof(int) );

   /* go METIS go ! */
   t1=clock();
   partitionCS( &giant );
   t2=clock();
   printf("Time to METIS global  %g sec\n", (double)(t2-t1)/CLOCKS_PER_SEC );

   write_vtu( &giant, "globe.vtu" );
   write_ascii( &giant );


   for( v_i=0; v_i<giant.nVerts; v_i++ ) {
      free(giant.my_nadj[v_i]);
      free(giant.vattribs[v_i]);
      free(giant.my_coords[v_i]);
   }
   free(giant.epart);
   free(giant.vpart);
   free(giant.my_nadj);
   free(giant.vattribs);
   free(giant.my_coords);
   for( e_i=0; e_i<nGlobalEls; e_i++ ) {
      free( giant.my_e_n[e_i]);
   }
   free(giant.my_e_n);
   free(giant.activeList);

   for( s_i=0; s_i<1; s_i++ ) {
      Sixth_FreeMem( &sixths[s_i] );
   }

   printf("Success\n");

}

void write_sixth_vtu( Sixth* self, char* name )
{
   /*
      function just writes a Sixth* to a given file called filename in vtu format
   */

   FILE *oFile=NULL;
   unsigned v_i, e_i, n_i;
   int num_conn=0; // the number of connections
   int nVerts=self->nVerts;
   int nEls=self->nEls;
   double **verts=self->my_coords;
   int **nbrs=self->my_nadj;
   int **e_n=self->my_e_n;
   int dim = self->dim;

   oFile=fopen(name, "w");
   assert(oFile);

   /* count the number of vtk connections */
   // 1st vert adjacency
   for( v_i=0; v_i<self->nVerts; v_i++ ) {
      if( nbrs[v_i] == NULL ) continue; //skip if non-existant: case for deleted verts
      for( n_i=0; n_i<6; n_i++ ) {
         if( nbrs[v_i][n_i] != -1 ) {
            num_conn++;
         }
      }
   }

   // write vtu header
   fprintf(oFile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
   fprintf(oFile, "<UnstructuredGrid>\n<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nVerts, num_conn+self->nEls);
   fprintf(oFile, "<Points>\n");
   fprintf(oFile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");

   for( v_i=0; v_i<nVerts; v_i++ ) {
      fprintf(oFile,  "%g %g %g\n", verts[v_i][0], verts[v_i][1], verts[v_i][2] );
   }

   fprintf(oFile, "</DataArray>\n");
   fprintf(oFile, "</Points>\n");

   fprintf(oFile, "<Cells>\n");
   fprintf(oFile, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");

   // do node nbr connectivity
   for( v_i=0; v_i<nVerts; v_i++ ) {
      if( nbrs[v_i] == NULL ) continue; //skip if non-existant
      for( n_i=0; n_i<6; n_i++ ) {
         if( nbrs[v_i][n_i] != -1 ) {
            fprintf(oFile, "%d %d\n", v_i, nbrs[v_i][n_i]-self->v_offset );
         }
      }
   }
   // do element node connectivity
   int nodesPerEl=8;
   if( dim==2 ) nodesPerEl=4;
   for( v_i=0; v_i<nEls; v_i++) {
      for(  n_i=0; n_i<nodesPerEl; n_i++ ) {
         fprintf( oFile, "%d ", e_n[v_i][n_i]-self->v_offset); // print node in element
      }
      fprintf(oFile, "\n");
   }

   fprintf(oFile,"</DataArray>\n");

   fprintf(oFile,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
   for( v_i=1; v_i<num_conn+1; v_i++) {
      fprintf(oFile, "%d\n", 2*v_i );
   }
   for( n_i=1; n_i<nEls+1; n_i++) {
      fprintf(oFile, "%d\n", (2*num_conn)+nodesPerEl*n_i );
   }
   fprintf(oFile,"</DataArray>\n");
   fprintf(oFile,"<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n");
   for( v_i=0; v_i<num_conn; v_i++) {
      fprintf(oFile, "%d\n", 3 ); // magic definition for VTK_LINE
   }
   if( dim == 3 ) {
       for( v_i=0; v_i<nEls; v_i++) {
          fprintf(oFile, "%d\n", 12 ); // magic definition for VTK_HEXAHEDRA
       }
   } else {
       for( v_i=0; v_i<nEls; v_i++) {
          fprintf(oFile, "%d\n", 9 ); // magic definition for VTK_HEXAHEDRA
       }
   }
   fprintf(oFile,"</DataArray>\n");

   fprintf(oFile,"</Cells>\n");
   fprintf(oFile,"</Piece>\n");
   fprintf(oFile,"</UnstructuredGrid>\n");
   fprintf(oFile,"</VTKFile>\n");

   fclose(oFile);

}

void partitionCS( CS *mesh )
{
   /*
      create compressed storage format (CSR) of mesh and METIS it
    */

   CS *self=mesh;

   int nVerts=mesh->nVerts;
   int nEls=mesh->nEls;
   int **e_n=mesh->my_e_n;
   int *epart, *vpart;
   int dim = mesh->dim;
   int nParts = mesh->nParts;
   int nodesPerEl;

   epart=mesh->epart;
   vpart=mesh->vpart;

   idx_t *eptr=malloc( (nEls+1)*sizeof(idx_t) );
   idx_t *eind=malloc( nEls*8*sizeof(idx_t) );
   idx_t *part=malloc( nVerts*sizeof(idx_t) ); // the vert
   idx_t ncon=1;
   real_t *tpwgts=NULL;
   idx_t v_i, nbr_i;
   idx_t num_conn=0;

   idx_t vtxdist=nVerts;
   idx_t elmdist=nEls;

   idx_t ncommonnodes = 4; 
   nodesPerEl = 8;
   if( dim == 2 ) {
    ncommonnodes = 2; // for quad
    nodesPerEl = 4;
   }

   // fill in element indices
   idx_t num_eind=0;
   idx_t e_i, t_i;
   for( e_i=0; e_i<nEls; e_i++ ) {
      eptr[e_i]=num_eind;
      for(nbr_i=0; nbr_i<nodesPerEl; nbr_i++ ) {
         eind[num_eind] = e_n[e_i][nbr_i];
         num_eind++;
      }
   }
   eptr[e_i]=num_eind;  // the final entry

   idx_t ncommon=4;
   idx_t nparts=nParts;

   if( dim == 2 )
       ncommon = 2;

   tpwgts = malloc( ncon*nparts*sizeof(real_t) );
   memset( tpwgts, 1, nparts*sizeof(real_t) );
   idx_t options[METIS_NOPTIONS];
   idx_t objval;


   METIS_SetDefaultOptions(options);
   options[METIS_OPTION_OBJTYPE] =METIS_OBJTYPE_VOL; // METIS_OBJTYPE_VOL; METIS_OBJTYPE_CUT;
   options[METIS_OPTION_NUMBERING] = 0; // c-style
   for( t_i=0; t_i<ncon*nparts; t_i++ ) {
      tpwgts[t_i]=1.0/(double)nparts;
   }

   /*
   int err = METIS_PartMeshNodal(
                (idx_t*)&nEls,    // number of elements
                (idx_t*)&nVerts,  // number of verts
                eptr,
                eind,
                NULL, // vwgt
                NULL, // vsize
                &nparts,
                NULL, //tpwgts
                options, //options....
                &objval,
                (idx_t*)epart,
                (idx_t*)vpart );
                */

   int err = METIS_PartMeshDual(
                (idx_t*)&nEls,    // number of elements
                (idx_t*)&nVerts,  // number of verts
                eptr,
                eind,
                NULL, // vwgt
                NULL, // vsize
                &ncommon, // 4
                &nparts,
                NULL, //tpwgts
                NULL, //options....
                &objval,
                (idx_t*)epart,
                (idx_t*)vpart );

   if( err==METIS_OK ) {
      printf("Wow, METIS_OK - the number of edges is %d\n", objval);
   } else {
      printf("METIS is not happy\n");
   }

   free( tpwgts );
   free( eptr );
   free( eind );
   free( part );
}

void write_vtu( CS *cs, char* filename )
{
   /*
      write_vtu() writes the CS data struct to a vtu file called, filename
   */

   FILE *file;
   double **verts;
   int **nbrs, nVerts, **e_n, nEls, *epart, *vpart;
   int dim = cs->dim;

   verts = cs->my_coords;
   nbrs = cs->my_nadj;
   e_n = cs->my_e_n;
   nVerts = cs->nVerts;
   nEls = cs->nEls;
   epart = cs->epart;
   vpart = cs->vpart;

   file = fopen(filename, "w");
   assert(file);

   unsigned e_i,v_i, n_i;
   int num_conn=0;
   int pervert=0;
   for( v_i=0; v_i<nVerts; v_i++ ) {
      pervert=0;
      if( nbrs[v_i] == NULL ) continue;

      for( n_i=0; n_i<6; n_i++ ) {
         if( nbrs[v_i][n_i] != -1 ) {
            pervert++;
            num_conn++;
         }
      }
   }

   fprintf(file, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
   fprintf(file, "<UnstructuredGrid>\n<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nVerts, (/*num_conn+*/nEls));

   /********* code to add scalars to points **************/
   fprintf(file, "<PointData Scalars=\"scalars\">\n");
   fprintf(file, "<DataArray type=\"Int32\" Name=\"bc\" format=\"ascii\">\n");
   for( v_i=0; v_i<nVerts; v_i++ ) {
      fprintf(file, "%d ", cs->vattribs[v_i][VA_BC] );
   }
   fprintf(file, "\n</DataArray>\n");
   fprintf(file, "<DataArray type=\"Int32\" Name=\"nParts\" format=\"ascii\">\n");
   for( v_i=0; v_i<nVerts; v_i++ ) {
      fprintf(file, "%d ", cs->vpart[v_i] );
   }
   fprintf(file, "\n</DataArray>\n");

   fprintf(file, "</PointData>\n");
   /**************************/

   /********* code to add scalars to cells **************/
   fprintf(file, "<CellData Scalars=\"scalars\">\n");
   fprintf(file, "<DataArray type=\"Int32\" Name=\"partition\" format=\"ascii\">\n");
   for( e_i=0; e_i<nEls; e_i++ ) {
      fprintf(file, "%d ", epart[e_i] );
   }
   fprintf(file, "\n</DataArray>\n");

   fprintf(file, "</CellData>\n");
   /**************************/


   /********* add points ****************/
   fprintf(file, "<Points>\n");
   fprintf(file, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");

   for( v_i=0; v_i<nVerts; v_i++ ) {
      fprintf(file,  "%g %g %g\n", verts[v_i][0], verts[v_i][1], verts[v_i][2] );
   }

   fprintf(file, "</DataArray>\n");
   fprintf(file, "</Points>\n");

   /******** cell data, is the element-node connectivity *****/
   fprintf(file, "<Cells>\n");
   fprintf(file, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");

   /*
   for( v_i=0; v_i<nVerts; v_i++ ) {
      if( nbrs[v_i] == NULL ) continue;

      for( n_i=0; n_i<6; n_i++ ) {
         if( nbrs[v_i][n_i] != -1 ) {
            fprintf(file, "%d %d\n", v_i, nbrs[v_i][n_i] );
         }
      }
   }
   */
   int nodesPerCell = 8;
   if( dim == 2 ) nodesPerCell = 4;
   for( v_i=0; v_i<nEls; v_i++) {
      for(  n_i=0; n_i<nodesPerCell; n_i++ ) {
         fprintf( file, "%d ", e_n[v_i][n_i]);
      }
      fprintf(file, "\n");
   }

   fprintf(file,"</DataArray>\n");

   fprintf(file,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
   /*
   for( v_i=1; v_i<num_conn+1; v_i++) {
     fprintf(file, "%d\n", 2*v_i );
   }
   */
   for( n_i=1; n_i<nEls+1; n_i++) {
      fprintf(file, "%d\n", /*(2*num_conn)+*/nodesPerCell*n_i );
   }
   fprintf(file,"</DataArray>\n");

   fprintf(file,"<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n");
   /*
   for( v_i=0; v_i<num_conn; v_i++) {
     fprintf(file, "%d\n", 3 ); // magic definition for VTK_LINE
   }
   */
   int magic_no = 12;
   if( dim == 2 ) magic_no = 9;
   for( v_i=0; v_i<nEls; v_i++) {
      fprintf(file, "%d\n", magic_no ); // magic definition for VTK_HEXAHEDRA
   }
   fprintf(file,"</DataArray>\n");

   fprintf(file,"</Cells>\n");
   fprintf(file,"</Piece>\n");
   fprintf(file,"</UnstructuredGrid>\n");
   fprintf(file,"</VTKFile>\n");

   fclose(file);
}

void write_ascii( CS *cs )
{
   /*
      write the mesh information
      vert coord - coords.dat
      element_vert - elements.dat
      vert_vert    - vert_conn.dat
   */

   FILE *file;
   double **verts;
   int **nbrs, nVerts, **vattribs, **e_n, nEls, *epart, *vpart;
   unsigned e_i,v_i, n_i;
   int num_conn=0;
   int pervert=0;

   verts = cs->my_coords;
   nbrs = cs->my_nadj;
   e_n = cs->my_e_n;
   nVerts = cs->nVerts;
   nEls = cs->nEls;
   epart = cs->epart;
   vpart = cs->vpart;
   vattribs = cs->vattribs;

   file = fopen("coords.dat", "w");
   assert(file);

   // output the total number of verts
   fprintf(file, "%d\n", nVerts );

   for(v_i=0; v_i<nVerts; v_i++ ) {
      fprintf(file, "%g %g %g %d\n", verts[v_i][0], verts[v_i][1], verts[v_i][2], vattribs[v_i][VA_BC] );
   }
   fclose(file);

   file = fopen("elements.dat", "w");
   assert(file);

   // output the total number of verts
   fprintf(file, "%d\n", nEls );

   for(e_i=0; e_i<nEls; e_i++ ) {
      for( v_i=0; v_i<8; v_i++ ) {
         fprintf(file, "%d ", e_n[e_i][v_i] );
      }
      fprintf(file, "\n");
   }
   fclose(file);

   file = fopen("vert_conn.dat", "w");
   assert(file);

   for(v_i=0; v_i<nVerts; v_i++ ) {
      for( n_i=0; n_i<6; n_i++ ) {
         fprintf(file, "%d ", nbrs[v_i][n_i]);
      }
      fprintf(file, "\n");
   }
   fclose(file);
}

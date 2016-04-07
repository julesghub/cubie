#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<unistd.h>
#include<time.h>

#define NOS 6

/* data structure for complete cubed sphere, CS */
typedef struct {
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

#if METIS_ENABLED
#include<metis.h>
void partitionCS( CS *mesh );
#endif

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


void Sixth_Setup( Sixth* self, unsigned e_offset, unsigned v_offset, unsigned* inputSizes, double *xi_limits, double *eta_limits, double *r_limits )
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

   memcpy( self->sizes, inputSizes, 3*sizeof(unsigned) );

   self->basis[0] = 1;
   self->basis[1] = self->sizes[0];
   self->basis[2] = self->sizes[0]*self->sizes[1];

   self->elSizes[0] = self->sizes[0]-1;
   self->elSizes[1] = self->sizes[1]-1;
   self->elSizes[2] = self->sizes[2]-1;

   // used to create globalIds
   self->v_offset = v_offset;
   self->e_offset = e_offset;

   self->nvattribs=1;

   self->nEls = self->elSizes[0] * self->elSizes[1] * self->elSizes[2];
   self->nVerts = self->sizes[0]*self->sizes[1]*self->sizes[2];
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
   assert( ijk[2] < self->sizes[2] );

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

   for( d_i = 3; d_i > 0; d_i-- ) {
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

void join_wall_verts( Sixth *g1, enum NBR_TYPE face1_def, Sixth *g2, enum NBR_TYPE face2_def, int flip )
/*@
  A general function to link to sixths together.
  The g1 sixth is a dominant sixth one, duplicate verts from g2 are deactivated as the join is made.

  1st) the interfacing verts of each sixth are built
  2nd) g2's verts are deactivated
  3rd) g2's elements swap deactived g2 verts for active g1 verts
  4th) update g1 and g2 nbr adjacency list to contain global active verts

  the order of how join_wall_verts is called for all joining walls is important!

  @*/

{
   unsigned fixed1_dof, link1_dof, fixed1_value;
   unsigned fixed2_dof, link2_dof, fixed2_value;
   int *inorder_wall_g1, *inorder_wall_g2;
   unsigned n_surface_verts, n_depth_verts;

   int newVertId, oldVertId, lvid, g1_lvid, nbrId, tmp_i;
   int tmplocal;
   unsigned ijk[3], abc[3];

   n_surface_verts = g1->sizes[0];
   n_depth_verts = g1->sizes[2];

   /* set up how to iterate on face1 */
   if( face1_def == NBR_N ) {
      fixed1_dof = 1;
      link1_dof = 0;
      fixed1_value = g1->sizes[fixed1_dof]-1;
   } else if( face1_def == NBR_S ) {
      fixed1_dof = 1;
      link1_dof = 0;
      fixed1_value = 0;
   } else if ( face1_def == NBR_E ) {
      fixed1_dof = 0;
      link1_dof = 1;
      fixed1_value = g1->sizes[fixed1_dof]-1;
   } else if( face1_def == NBR_W ) {
      fixed1_dof = 0;
      link1_dof = 1;
      fixed1_value = 0;
   }

   /* set up how to iterate on face2 */
   if( face2_def == NBR_N ) {
      fixed2_dof = 1;
      link2_dof = 0;
      fixed2_value = g2->sizes[fixed2_dof]-1;
   } else if( face2_def == NBR_S ) {
      fixed2_dof = 1;
      link2_dof = 0;
      fixed2_value = 0;
   } else if ( face2_def == NBR_E ) {
      fixed2_dof = 0;
      link2_dof = 1;
      fixed2_value = g2->sizes[fixed2_dof]-1;
   } else if( face2_def == NBR_W ) {
      fixed2_dof = 0;
      link2_dof = 1;
      fixed2_value = 0;
   }

   /* allocate list of wall verts */
   inorder_wall_g1 = malloc( n_surface_verts*n_depth_verts * sizeof(int) );
   inorder_wall_g2 = malloc( n_surface_verts*n_depth_verts * sizeof(int) );

   /* set fixed index */
   ijk[fixed1_dof] = fixed1_value;
   abc[fixed2_dof] = fixed2_value;


   /* create sets of joining wall nodes:
       inorder_wall_g1 represents the nodes involved from sixth g1, using global vert ids
       inorder_wall_g2 represents the nodes involved from sixth g2, using global vert ids
    */
   int kk,jj,node_i,v_i,count=0;
   for( kk=0; kk<n_depth_verts; kk++ ) {
      abc[2]=ijk[2]=kk; // depth is always the outer loop
      for( jj=0; jj<n_surface_verts; jj++ ) {
         ijk[link1_dof] = jj;
         if( flip ) abc[link2_dof] = (n_surface_verts-1)-jj; //if we are backwards on one wall
         else       abc[link2_dof] = jj;

         inorder_wall_g1[count] = Sixth_Project( g1, ijk, GLOBAL );
         inorder_wall_g2[count] = Sixth_Project( g2, abc, GLOBAL );
         count++;
      }
   }

   /*
      disactivate the vert adjacecy lists of the g2 verts no longer used
   */
   for( v_i=0; v_i < count; v_i++ ) {
      oldVertId = inorder_wall_g2[v_i]-g2->v_offset; // make local index
      g2->activeList[oldVertId] = 0;
   }

   /* fixed indices for elements, either 0 or fixed*_value-1 */
   ijk[fixed1_dof] = (fixed1_value != 0) ? (fixed1_value-1) : 0;
   abc[fixed2_dof] = (fixed2_value != 0) ? (fixed2_value-1) : 0;

   unsigned e_g2, nbr_i;
   int nid, active;

   /* update elements in g2 to contain verts in g1 */
   for( kk=0 ; kk< n_depth_verts-1; kk++) {
      abc[2]=ijk[2]=kk;
      for( jj=0 ; jj< n_surface_verts-1; jj++) {
         ijk[link1_dof] = jj;
         if( flip ) abc[link2_dof] = (n_surface_verts-2)-jj; //if we are backwards on one wall note -2
         else       abc[link2_dof] = jj;

         e_g2 = Sixth_eProject( g2, abc, LOCAL );

         /* find nodes in e_g2, that are in inorder_wall_g2 */
         for( node_i=0; node_i<8 ; node_i++ ) {
            nid = g2->my_e_n[e_g2][node_i];
            if( !Sixth_OwnsVert(g2, nid ) ) {
               //TODO    printf("won't update cause I don't own it\n");
               continue;
            }
            for( v_i=0; v_i<count; v_i++ ) { //search all - not efficient but gets the job done
               if( nid == inorder_wall_g2[v_i] ) {
                  // now swap the g2 vert with a g1 vert
                  g2->my_e_n[e_g2][node_i] = inorder_wall_g1[v_i];
                  break;
               }
            }
         }
      }
   }


   /*
      update g1 and g2 nbr adjacency list. so...
      1) no verts in g2 touches the verts in the inorder_wall_g2 set
         instead they touch the inorder_wall_g1 set
      2) verts in inorder_wall_g1 set touch with the appropriate g2 verts
    */
   int isIn1, isIn2;
   for( v_i=0; v_i < count ; v_i++ ) {
      newVertId = inorder_wall_g1[v_i];

      oldVertId = inorder_wall_g2[v_i];
      tmplocal = oldVertId - g2->v_offset; // convert to local id

      // find oldVertId's nbring verts to be updated to newVertId
      for( nbr_i=0; nbr_i<6; nbr_i++ ) {
         nbrId = g2->my_nadj[tmplocal][nbr_i]; // get global id of g2 nbr verts

         if( nbrId == -1 ) continue; // if empty connection do nothing

         isIn1 = Sixth_OwnsVert(g1, nbrId);
         isIn2 = Sixth_OwnsVert(g2, nbrId);

         if( !isIn1 && !isIn2 ) {
            //TODO     printf("Triple junction on global node %d - doing not updates\n", nbrId );
            continue;
         }

         lvid = nbrId - g2->v_offset; // convert to local id

         if( g2->activeList[lvid] == 0 ) {
            // if inactive nbr, skip change it's adjacency list
            //  g2->my_nadj[tmplocal][nbr_i]=-1;
            continue;

         }
         // find correct nbrs vert to update
         for( tmp_i=0; tmp_i<6; tmp_i++ ) {
            if( g2->my_nadj[lvid][tmp_i] == oldVertId ) {
               // update to vert from g1
               g2->my_nadj[lvid][tmp_i] = newVertId;

               //  g2->my_nadj[tmplocal][nbr_i] = -1; // disconnect
               break;
            }
         }

         g1_lvid = newVertId - g1->v_offset; // convert to local id
         // update the adjacency list of the g1 vert. Don't worry about ordering information in nadj list
         for( tmp_i=0; tmp_i<6; tmp_i++ ) {
            if (g1->my_nadj[g1_lvid][tmp_i] == -1 ) {
               g1->my_nadj[g1_lvid][tmp_i] = nbrId;
               break;
            }
         }
      }
   }

   free(inorder_wall_g1);
   free(inorder_wall_g2);
}


void reorderCS( CS* giant, int oopps )
{
   /*@
      used to reorder indices so the data structures can be contiguous
      and vtu output has no redundant information like nodes that don't
      belong to any element
    @*/

   int *activeList = giant->activeList;
   int *mapping = malloc( giant->nVerts*sizeof(int) );
   int **nbr=NULL; // the new adjacency table
   int **vattribs=NULL; // the new attributes table
   double **coords=NULL; // new coord array

   int activeCount;
   int v_i, nbr_i, oldId, newId;
   int mapid;

   // create mapping function
   activeCount=0;
   for( v_i=0; v_i< giant->nVerts; v_i++ ) {
      if( activeList[v_i] == 0 ) {
         mapping[v_i]=-1; // definte invalid as -1
         continue;
      }

      // mapping for oldId to newId
      mapping[v_i]=activeCount;
      activeCount++;

   }

   // allocate memory for replacement arrays
   nbr = malloc( activeCount*sizeof(int*) );
   vattribs = malloc( activeCount*sizeof(int*) );
   coords = malloc( activeCount* sizeof(double*) );
   for( v_i=0; v_i<activeCount; v_i++ ) {
      nbr[v_i] = malloc( 6*sizeof(int) );
      vattribs[v_i] = malloc( giant->nvattribs*sizeof(int) );
      coords[v_i] = malloc( 3*sizeof(double) );
   }

   int e_i;
   for( e_i=0; e_i < giant->nEls; e_i++ ) {
      for( v_i=0; v_i<8; v_i++ ) {
         oldId = giant->my_e_n[e_i][v_i];
         newId = mapping[oldId];

         if( newId == -1 ) {
            printf("Why the fuck\n");
         } else {
            giant->my_e_n[e_i][v_i] = newId;
         }
      }
   }

   /* now rebuilt the vert adjacency information with new indices */
   activeCount=0;
   for( v_i=0; v_i< giant->nVerts; v_i++ ) {
      if( activeList[v_i] == 0 ) continue; // skip vert if inactive

      // go through nbr list and redefine it using mapping
      for( nbr_i=0; nbr_i<6; nbr_i++ ) {
         oldId = giant->my_nadj[v_i][nbr_i];

         // no mapping required if oldId is -1
         if( oldId == -1 ) {
            nbr[activeCount][nbr_i] = -1;
            continue;
         }

         newId = mapping[oldId];
         nbr[activeCount][nbr_i] = newId;

      }
      activeCount++;
   }

   /* replace giant->my_coords, giant->my_nadj */
   for( v_i=0; v_i<giant->nVerts; v_i++ ) {
      // copy coordinate into new array
      mapid = mapping[v_i];

      if( mapid == -1 ) continue; // if inactive node don't copy it to giant
      memcpy(coords[mapid], giant->my_coords[v_i], 3*sizeof(double) );
      memcpy(vattribs[mapid], giant->vattribs[v_i], giant->nvattribs*sizeof(int) );
   }

   // free old giant definitions
   for(v_i=0; v_i<oopps; v_i++ ) {
      free(giant->my_nadj[v_i]);
      free(giant->vattribs[v_i]);
      free(giant->my_coords[v_i]);
   }
   free( giant->my_nadj );
   free( giant->vattribs );
   free( giant->my_coords );

   // replace ptrs
   giant->my_coords = coords;
   giant->my_nadj = nbr;
   giant->vattribs = vattribs;
   //replace count
   giant->nVerts=activeCount;

   free(mapping);
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

   // build sixth's rotation matrix
   build_rotation_matrix( x_rot, y_rot, rot );
   v_i=0;

   // denominator is (size[]-1) because we are calculating the spacing
   d_xi = fabs( (self->xi_limits[1]-self->xi_limits[0])/(double)(self->sizes[0]-1) );
   d_eta = fabs( (self->eta_limits[1]-self->eta_limits[0])/(double)(self->sizes[1]-1) );
   d_r = fabs( (self->r_limits[1]-self->r_limits[0])/(double)(self->sizes[2]-1) );


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
            a = r/d * X;
            b = r/d * Y;
            c = r/d * 1;

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

void txtOutput( int n_surface_elements, int n_depth_elements, double inner_radius, double outer_radius, CS* sphere ) {



}

int main(int argc, char** argv )
{
   unsigned n_surface_verts, n_depth_verts;
   unsigned n_surface_els, n_depth_els;
   double outer_radius, inner_radius;
   Sixth sixths[6];
   CS giant;
   unsigned sixthSize[3], ijk[3], tmp;
   int nVerts, nEls, s_i, v_i, n_id, e_i;
   int **nbrs;
   double **verts;
   double angleDelta;
   double r_limits[2], eta_limits[2], xi_limits[2];

   int **e_n, elsBuilt;

   double **global_coord;
   int nGlobalEls;
   int nGlobalVerts;
   int c;

   clock_t t2, t1 = clock();


   /***
     get input:
     geometry and number of elements
   ****/

   n_surface_els=-1;
   n_depth_els=-1;
   while ( (c = getopt(argc, argv, "l:d:")) != -1) {
      switch( c ) {
      case 'l':
         n_surface_els=atoi(optarg);
         break;
      case 'd':
         n_depth_els=atoi(optarg);
         break;
      case '?':
         printf("execute with\n ./go -l 11 -d 6\n");
         printf("where l is the length resolution and d is the depth resolution\n");
         exit(1);
      default:
         printf("execute with\n ./go -l 11 -d 6\n");
         printf("where l is the length resolution and d is the depth resolution\n");
         exit(1);
      }

   }
   if( n_surface_els==-1 || n_depth_els==-1 ) {
      printf("execute with\n ./go -l 11 -d 6\n");
      printf("where l is the number of elements in length (d depth) per sixth\n");
      exit(1);
   }



   /***
     get e_n, epart, n_n, npart
   ***/
   /* set initial sixth definition */
   n_surface_verts = n_surface_els+1;
   n_depth_verts = n_depth_els+1;
   outer_radius=6;
   inner_radius=3;

   angleDelta = (M_PI/2)/(double)(n_surface_verts-1);

   nGlobalEls = 6*n_surface_els*n_surface_els*n_depth_els;
   nGlobalVerts = 6*n_surface_verts*n_surface_verts*n_depth_verts;

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

   sixthSize[2] = n_depth_verts;
   r_limits[0] = inner_radius;
   r_limits[1] = outer_radius;
   sixthSize[0] = n_surface_verts;
   sixthSize[1] = n_surface_verts;
   xi_limits[0]  = -M_PI/4;
   xi_limits[1]  = M_PI/4;
   eta_limits[0] = -M_PI/4;
   eta_limits[1] = M_PI/4;

   Sixth_Setup(&(sixths[0]), 0, 0, sixthSize, xi_limits, eta_limits, r_limits );
   build_sixth_nadj( &sixths[0] );
   build_sixth_e_n( &sixths[0] );
   evaluate_sixth_geom( &sixths[0],0,0 );

   Sixth_Setup(&(sixths[1]), sixths[0].nEls, sixths[0].nVerts, sixthSize, xi_limits, eta_limits, r_limits );
   build_sixth_nadj( &sixths[1] );
   build_sixth_e_n( &sixths[1] );
   evaluate_sixth_geom( &sixths[1],0,-M_PI/2 );

   Sixth_Setup(&(sixths[2]), 2*sixths[0].nEls, 2*sixths[0].nVerts, sixthSize, xi_limits, eta_limits, r_limits );
   build_sixth_nadj( &sixths[2] );
   build_sixth_e_n( &sixths[2] );
   evaluate_sixth_geom( &sixths[2],0,-M_PI );

   Sixth_Setup(&(sixths[3]), 3*sixths[0].nEls, 3*sixths[0].nVerts, sixthSize, xi_limits, eta_limits, r_limits );
   build_sixth_nadj( &sixths[3] );
   build_sixth_e_n( &sixths[3] );
   evaluate_sixth_geom( &sixths[3],0,M_PI/2 );

   Sixth_Setup(&(sixths[4]), 4*sixths[0].nEls, 4*sixths[0].nVerts, sixthSize, xi_limits, eta_limits, r_limits );
   build_sixth_nadj( &sixths[4] );
   build_sixth_e_n( &sixths[4] );
   evaluate_sixth_geom( &sixths[4],-M_PI/2,0 );

   Sixth_Setup(&(sixths[5]), 5*sixths[0].nEls, 5*sixths[0].nVerts, sixthSize, xi_limits, eta_limits, r_limits );
   build_sixth_nadj( &sixths[5] );
   build_sixth_e_n( &sixths[5] );
   evaluate_sixth_geom( &sixths[5],M_PI/2,0 );

   t2 = clock();
   printf("Time to initialise sixth %g sec\n", (double)(t2-t1)/CLOCKS_PER_SEC );

   /*
   // debug
   write_sixth_vtu( &sixths[0], "sixth-0.vtu" );
   write_sixth_vtu( &sixths[1], "sixth-1.vtu" );
   write_sixth_vtu( &sixths[2], "sixth-2.vtu" );
   write_sixth_vtu( &sixths[3], "sixth-3.vtu" );
   write_sixth_vtu( &sixths[4], "sixth-4.vtu" );
   write_sixth_vtu( &sixths[5], "sixth-5.vtu" );
   */

   /*
      run function to join walls
   */

   t1=clock();

   join_wall_verts( &sixths[0], NBR_E, &sixths[1], NBR_W, 0 );
   join_wall_verts( &sixths[1], NBR_E, &sixths[2], NBR_W, 0 );
   join_wall_verts( &sixths[2], NBR_E, &sixths[3], NBR_W, 0 );
   join_wall_verts( &sixths[0], NBR_W, &sixths[3], NBR_E, 0 );

   join_wall_verts( &sixths[0], NBR_N, &sixths[4], NBR_S, 0 );
   join_wall_verts( &sixths[1], NBR_N, &sixths[4], NBR_E, 0 );
   join_wall_verts( &sixths[2], NBR_N, &sixths[4], NBR_N, 1 );
   join_wall_verts( &sixths[3], NBR_N, &sixths[4], NBR_W, 1 );

   join_wall_verts( &sixths[0], NBR_S, &sixths[5], NBR_N, 0 );
   join_wall_verts( &sixths[1], NBR_S, &sixths[5], NBR_E, 1 );
   join_wall_verts( &sixths[2], NBR_S, &sixths[5], NBR_S, 1 );
   join_wall_verts( &sixths[3], NBR_S, &sixths[5], NBR_W, 0 );

   t2=clock();
   printf("Time to join sixths %g sec\n", (double)(t2-t1)/CLOCKS_PER_SEC );

   /* add all sixths together into a giant mesh */
   t1=clock();
   addSixthToGiant( &giant, &sixths[0] );
   addSixthToGiant( &giant, &sixths[1] );
   addSixthToGiant( &giant, &sixths[2] );
   addSixthToGiant( &giant, &sixths[3] );
   addSixthToGiant( &giant, &sixths[4] );
   addSixthToGiant( &giant, &sixths[5] );
   t2=clock();
   printf("Time to add sixths to global  %g sec\n", (double)(t2-t1)/CLOCKS_PER_SEC );

   /* reorder the indices of the giant mesh as the removal of duplicate nodes
      - made during the join_wall_verts func - means the giant data structure is
      not contiguous.
      Reording makes is contiguous again
    */
   t1=clock();
   reorderCS(&giant, nGlobalVerts );
   t2=clock();
   printf("Time to add sixths to global  %g sec\n", (double)(t2-t1)/CLOCKS_PER_SEC );

   giant.epart = malloc( giant.nEls*sizeof(int) );
   giant.vpart = malloc( giant.nVerts*sizeof(int) );

#if METIS_ENABLED
   /* go METIS go ! */
   t1=clock();
   partitionCS( &giant );
   t2=clock();
   printf("Time to METIS global  %g sec\n", (double)(t2-t1)/CLOCKS_PER_SEC );
#else
   printf("no metis action because it wasn't found in compile\n" );
#endif

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

   for( s_i=0; s_i<6; s_i++ ) {
      Sixth_FreeMem( &sixths[s_i] );
   }

   printf("Success\n");
   exit(0);

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

   for( v_i=0; v_i<nVerts; v_i++ ) {
      if( nbrs[v_i] == NULL ) continue; //skip if non-existant
      for( n_i=0; n_i<6; n_i++ ) {
         if( nbrs[v_i][n_i] != -1 ) {
            fprintf(oFile, "%d %d\n", v_i, nbrs[v_i][n_i]-self->v_offset );
         }
      }
   }
   for( v_i=0; v_i<nEls; v_i++) {
      for(  n_i=0; n_i<8; n_i++ ) {
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
      fprintf(oFile, "%d\n", (2*num_conn)+8*n_i );
   }
   fprintf(oFile,"</DataArray>\n");

   fprintf(oFile,"<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n");
   for( v_i=0; v_i<num_conn; v_i++) {
      fprintf(oFile, "%d\n", 3 ); // magic definition for VTK_LINE
   }
   for( v_i=0; v_i<nEls; v_i++) {
      fprintf(oFile, "%d\n", 12 ); // magic definition for VTK_HEXAHEDRA
   }
   fprintf(oFile,"</DataArray>\n");

   fprintf(oFile,"</Cells>\n");
   fprintf(oFile,"</Piece>\n");
   fprintf(oFile,"</UnstructuredGrid>\n");
   fprintf(oFile,"</VTKFile>\n");

   fclose(oFile);

}

#if METIS_ENABLED
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

   idx_t ncommonnodes = 4; // for hexahedral

   // fill in element indices
   idx_t num_eind=0;
   idx_t e_i, t_i;
   for( e_i=0; e_i<nEls; e_i++ ) {
      eptr[e_i]=num_eind;
      for(nbr_i=0; nbr_i<8; nbr_i++ ) {
         eind[num_eind] = e_n[e_i][nbr_i];
         num_eind++;
      }
   }
   eptr[e_i]=num_eind;  // the final entry

   idx_t ncommon=4;
   idx_t nparts=4;

   tpwgts = malloc( ncon*nparts*sizeof(real_t) );
   memset( tpwgts, 1, nparts*sizeof(real_t) );
   idx_t options[METIS_NOPTIONS];
   idx_t objval;


   METIS_SetDefaultOptions(options);
   options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
   for( t_i=0; t_i<ncon*nparts; t_i++ ) {
      tpwgts[t_i]=1.0/(double)nparts;
   }

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
#endif

void write_vtu( CS *cs, char* filename )
{
   /*
      write_vtu() writes the CS data struct to a vtu file called, filename
   */

   FILE *file;
   double **verts;
   int **nbrs, nVerts, **e_n, nEls, *epart, *vpart;

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
   for( v_i=0; v_i<nEls; v_i++) {
      for(  n_i=0; n_i<8; n_i++ ) {
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
      fprintf(file, "%d\n", /*(2*num_conn)+*/8*n_i );
   }
   fprintf(file,"</DataArray>\n");

   fprintf(file,"<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n");
   /*
   for( v_i=0; v_i<num_conn; v_i++) {
     fprintf(file, "%d\n", 3 ); // magic definition for VTK_LINE
   }
   */
   for( v_i=0; v_i<nEls; v_i++) {
      fprintf(file, "%d\n", 12 ); // magic definition for VTK_HEXAHEDRA
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

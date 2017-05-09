/* lattice3-1-04
 * Evolution of a lattice with cell divisions
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>

#include "const.h" /* Main file where the parameters are defined */

#include "lattice.h" /* Declaration of the structures of the lattice and declaration of the functions defined in this file */
#include "locerror.h" /* Error function declaration */
#include "measurements.h" /* Declaration of the functions used for statistics */
#include "pngwrite.h"
#include "opti.h"
#include "compalone.h"
#include "kevin_tools.h"
#include "input.h"

/* Time is evaluated from effdivision() function */
long int itime=0;
double ttime;

/* # of t1 processes */
int it1process=0;

int last_first_lattice; /* variable used to build the regular lattice */
int nb_cell_tot; /* Total number of cells (including the exterior cell which is composed of the vertices on the boundary) */
int nb_clone_tot;
int nb_vertex_tot; /* Total number of vertices */
int nb_cell_killed; 
int nb_bond_tot;

/* Defines the division rate of clones and background :
 * 1 (or >1) is for a uniform division_rate
 * 0.5 makes it divide twice slower, etc
 * May be used in compalone.c
 */
double division_rates[2]={1, 1};
double cell_sizes[2]={1, 1};


unsigned int rand_disp=0;

unsigned int edges_dist[10]={0};

double add_factor = 0.0000;
double add_factor1 = 0.0000;
double add_factor2 = 0.0;

/* Pointer of the cell in the center of the initial lattice
*/
Dual *pcentercell;

/* Arrays of pointers to vertices or cells
*/
Node **web;
Dual **web_dual;
Bond **web_bond;

/* Structure used to decide whether a t1process is necessary
*/
struct t1struct {
   Node *pvertex;
   int ineighb;
} t1data;


/* Allocation functions
*/

/* Returns a pointer to a memory block of size dimension*(Bond*)*/
Bond **createbondvector(int dimension){
   Bond **a;
   a = (Bond**)malloc(dimension * sizeof(Bond*));
   if (!a) locerror("createbondvector", "Not enough memory");
   return a;
}

/* Reallocation (to different size)*/
Bond **reallocbondvector(Bond **vectortorealloc, int newdimension){
   Bond **a;
   a = (Bond**)realloc(vectortorealloc, newdimension*sizeof(Bond*));
   if (!a) locerror("reallocbondvector", "Not enough memory");
   return a;
}

/* Returns a pointer to a memory block of size dimension*(Node*)
*/
   Node
**createnodevector(int dimension)
{
   Node**a;
   a = (Node**)malloc(dimension * sizeof(Node*));
   if (!a) locerror("createnodevector", "Not enough memory");
   return a;
}

/* Reallocation (to a bigger size for example)
*/
   Node
**reallocnodevector(Node **vectortorealloc, int newdimension)
{
   Node**a;
   a = (Node**)realloc(vectortorealloc, newdimension * sizeof(Node*));
   if (!a) locerror("reallocnodevector", "Not enough memory");
   return a;
}

/* The same for cells
*/
   Dual
**createdualvector(int dimension)
{
   Dual **a;
   a = (Dual **)malloc(dimension * sizeof(Dual *));
   if (!a) locerror("createdualvector", "Not enough memory");
   return a;
}

   Dual
**reallocdualvector(Dual **vectortorealloc, int newdimension)
{
   Dual **a;
   a = (Dual **)realloc(vectortorealloc, newdimension * sizeof(Dual *));
   if (!a) locerror("reallocdualvector", "Not enough memory");
   return a;
}

/* Allocation of space for a vertex
*/
   Node
*allocnode()
{
   Node*a;
   a = (Node*)malloc(sizeof(Node));
   if (!a) locerror("allocnode", "Not enough memory");
   return a;
}

/* Allocation of space for a cell
*/
   Dual
*allocdual()
{
   Dual *a;
   a = (Dual *)malloc(sizeof(Dual));
   if (!a) locerror("allocdual", "Not enough memory");
   return a;
}

/* Allocation of space for a bond*/
Bond *allocbond(){
   Bond *a;
   a = (Bond *)malloc(sizeof(Bond));
   if (!a) locerror("allocbond", "Not enough memroy");
   return a;
}

/* Allacation of space for a list of vertices in a cell
*/
   Node
**allocvertexlist(int dimension)
{
   Node **a;
   a = (Node **)malloc(dimension * sizeof(Node *));
   if (!a) locerror("allocvertexlist", "Not enough memory");
   return a;
}

   Node
** reallocvertexlist(Node **listtorealloc, int newdimension)
{
   Node **a;
   a = (Node **)realloc(listtorealloc, newdimension * sizeof(Node *));
   if (!a) locerror("reallocvertexlist", "Not enough memory");
   return a;
}

/*Allocation of space for a list of neighboring cells*/
Dual **alloccelllist(int dimension){
   Dual **a;
   a = (Dual**)malloc(dimension*sizeof(Dual*));
   if (!a) locerror("alloccelllist", "Not enough memory");
   return a;
}

Dual **realloccelllist(Dual **listtorealloc, int newdimension){
   Dual **a;
   a = (Dual**)realloc(listtorealloc, newdimension*sizeof(Dual*));
   if (!a) locerror("realloccelllist", "Not enough memory");
   return a;
}

/*Allocation of space for list of boundary bonds*/
Bond **allocbondlist(int dimension){
   Bond **a;
   a = (Bond**)malloc(dimension * sizeof(Bond*));
   if (!a) locerror("allocbondlist","Not enough memory");
   return a;
}

Bond **reallocbondlist(Bond **listtorealloc, int newdimension){
   Bond **a;
   a = (Bond**)realloc(listtorealloc,newdimension*sizeof(Bond*));
   if (!a) locerror("reallocbondlist","Not enough memory");
   return a;
}

/* Free the memory used for a cell
*/
   void
freecell(Dual *pcell)
{
   free(pcell->vertexlist);
   free(pcell->celllist);
   free(pcell->bondlist);
   free(pcell);
}

/* Gives the indice of a vertix from the cartesian coordinates on the first lattice
*/
   int
vertex_from_coord(int row, int column)
{
   if (row>2*LENGTH_FIRST_LATTICE-1 || column>LENGTH_FIRST_LATTICE){
      locerror("vertex_from_coord", "coordinates went over the limits");
   }
   return LENGTH_FIRST_LATTICE*(row-1)+column;
}


/* Initializes the first regular lattice
*/
   void
init_lattice()
{
   int icell,j,k;
   Dual *pcell;
   char filename1[100],filename2[100];


   printf("before coord\n");
   coord_gen();//comment this out when using infoinput()
   printf("before shrink\n");
   sprintf(filename1,"info/0088.vertinfo");
   sprintf(filename2,"info/0088.cellinfo");
   //infoinput(filename1,filename2);
   //cut_cell_find();
   shrink_datafirst();//comment this out to use infoinput()
   //coord_mod();//to give the initial vertices an offset
   printf("after shrink\n");


   for (icell=0 ; icell<nb_cell_tot ; icell++){
      pcell=web_dual[icell];
      pcell->marker.size=0;
      pcell->marker.div_rate=0;
      pcell->marker.clone_index=0;
      pcell->area_soll=1;
      pcell->area_soll0=pcell->area_soll;
      pcell->marker.tension_index=1.0;
   }

   forcevector=(double**)malloc(3*sizeof(double*));
   for(k=0;k<3;k++){
      forcevector[k]=(double*)malloc(nb_vertex_tot*sizeof(double));
   }
   svglattice("svgs/0000.svg",1);
   latticeprint("dats/0000.dat",0);
   //cut_cell_multistep();
   for(icell=0;icell<nb_cell_tot;icell++){
      pcell=web_dual[icell];
      //if(icell==20 || icell==23 || icell==27 || icell==24 || icell==30){pcell->marker.tension_index+=add_factor;}
      //if(icell==0){pcell->area_soll=78;pcell->area_soll0=pcell->area_soll;}
      optimizet1();
   }
   itime = 0;
   it1process=0;
}

/* Free the memory used to store the lattice */
   void
kill_lattice()
{
   int i;

   for (i=0; i<nb_cell_tot; i++){
      freecell(web_dual[i]);
   }

   for (i=0; i<nb_vertex_tot; i++){
      free(web[i]);
   }

   free(web);
   free(web_dual);
}


/* Generation of the first lattice (some pointers to cells remains void and shrink_datafirst() must be applied)
*/
   void 
coord_gen()
{
   int row, column, column_lim, vertex, cell;
   Node *pvertex;
   Dual *pcell;

   last_first_lattice=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-1);

   web = createnodevector(last_first_lattice);

   web_dual=createdualvector(last_first_lattice);

   web_dual[0]=allocdual();
   web_dual[0] -> area =0;

   /* Allocation of memory for the vertices
    * The code is not optimized at all but it is only used for initialisation and it is more readable
    */

   for (row=2 ; row<2*LENGTH_FIRST_LATTICE-1; row++){

      if (row%2){
         column_lim=LENGTH_FIRST_LATTICE + 1;
         for (column=1 ; column<column_lim; column++){

            vertex=vertex_from_coord(row, column);
            web[vertex]=allocnode();
         }
      }

      if ((row+1)%2){
         column_lim=LENGTH_FIRST_LATTICE - 1;
         for (column=1 ; column<column_lim; column++){

            vertex=vertex_from_coord(row, column);
            web[vertex]=allocnode();
         }

         for (column=column_lim ; column < column_lim + 2; column++){
            vertex=vertex_from_coord(row, column);
            web[vertex]=NULL;
         }
      }
   }

   web[0]=NULL;
   web[1]=NULL;
   web[LENGTH_FIRST_LATTICE]=NULL;

   for (column=2 ; column < LENGTH_FIRST_LATTICE; column++){
      vertex=column;
      web[vertex]=allocnode();
   }

   web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+1]=NULL;

   for (column=2 ; column < LENGTH_FIRST_LATTICE; column++){
      vertex=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+column;
      web[vertex]=allocnode();
   }

   /* Allocation of memory for the cells */

   for (row=1 ; row<2*LENGTH_FIRST_LATTICE-1; row++){

      if (row%2){
         column=1;
         vertex=vertex_from_coord(row, column);
         web_dual[vertex]=NULL;

         column_lim=LENGTH_FIRST_LATTICE + 1;
         web_dual[LENGTH_FIRST_LATTICE*(row-1)+1]=NULL;
         for (column=2 ; column<column_lim; column++){
            vertex=vertex_from_coord(row, column);
            if (column%2){
               web_dual[vertex]=allocdual();
            }
            else {
               web_dual[vertex]=NULL;
            }
         }
      }

      if ((row+1)%2){
         column_lim=LENGTH_FIRST_LATTICE+1;
         for (column=1 ; column<column_lim; column++){
            vertex=vertex_from_coord(row, column);

            if (column%2){
               web_dual[vertex]=allocdual();
            }
            else {
               web_dual[vertex]=NULL;
            }
         }
      }
   }
   row=2*LENGTH_FIRST_LATTICE-1;
   column=1;
   vertex=vertex_from_coord(row, column);
   web_dual[vertex]=NULL;

   column_lim=LENGTH_FIRST_LATTICE;
   for (column=2 ; column<column_lim; column++){
      vertex=vertex_from_coord(row, column);
      if (column%2){
         web_dual[vertex]=allocdual();
      }
      else {
         web_dual[vertex]=NULL;
      }
   }


   web_dual[0]=allocdual();

   /* Effective construction of the lattice */

   for (row=2 ; row < 2*LENGTH_FIRST_LATTICE -1 ; row++){

      if (row%2){

         column_lim=LENGTH_FIRST_LATTICE;

         for (column=2 ; column<column_lim; column++){
            vertex=vertex_from_coord(row, column);
            pvertex=web[vertex];
            if (pvertex){

               pvertex -> y=(row-1)*sqrt(3)/2;

               if (column%2){
                  pvertex -> x=-1.5+1.5*column;

                  pvertex -> pneighb[0]=web[vertex+1];
                  pvertex -> pneighb[1]=web[vertex+LENGTH_FIRST_LATTICE-1];
                  pvertex -> pneighb[2]=web[vertex-LENGTH_FIRST_LATTICE-1];

                  pvertex -> pncell[0]=web_dual[vertex+LENGTH_FIRST_LATTICE];	
                  pvertex -> pncell[1]=web_dual[vertex];
                  pvertex -> pncell[2]=web_dual[vertex-LENGTH_FIRST_LATTICE];
               }

               if ((column+1)%2){
                  pvertex -> x=-2+1.5*column;

                  pvertex -> pneighb[0]=web[vertex-1];
                  pvertex -> pneighb[1]=web[vertex-LENGTH_FIRST_LATTICE-1];
                  pvertex -> pneighb[2]=web[vertex+LENGTH_FIRST_LATTICE-1];

                  pvertex -> pncell[0]=web_dual[vertex-LENGTH_FIRST_LATTICE-1];	
                  pvertex -> pncell[1]=web_dual[vertex+1];
                  pvertex -> pncell[2]=web_dual[vertex+LENGTH_FIRST_LATTICE-1];
               }

            }
            else {
               locerror("coord_gen", "The vertex does not exist");
            }
         }
      }


      if ((row+1)%2){

         column_lim=LENGTH_FIRST_LATTICE-1;

         for (column=1 ; column<column_lim; column++){
            vertex=vertex_from_coord(row, column);
            pvertex=web[vertex];
            if (pvertex){

               pvertex -> y=(row-1)*sqrt(3)/2;

               if (column%2){
                  pvertex -> x=1.5*column;

                  pvertex -> pneighb[0]=web[vertex+1];
                  pvertex -> pneighb[1]=web[vertex+LENGTH_FIRST_LATTICE+1];
                  pvertex -> pneighb[2]=web[vertex-LENGTH_FIRST_LATTICE+1];

                  pvertex -> pncell[0]=web_dual[vertex+LENGTH_FIRST_LATTICE+2];	
                  pvertex -> pncell[1]=web_dual[vertex];
                  pvertex -> pncell[2]=web_dual[vertex-LENGTH_FIRST_LATTICE+2];
               }

               if ((column+1)%2){
                  pvertex -> x=-0.5+1.5*column;

                  pvertex -> pneighb[0]=web[vertex-1];
                  pvertex -> pneighb[1]=web[vertex-LENGTH_FIRST_LATTICE+1];
                  pvertex -> pneighb[2]=web[vertex+LENGTH_FIRST_LATTICE+1];

                  pvertex -> pncell[0]=web_dual[vertex-LENGTH_FIRST_LATTICE+1];	
                  pvertex -> pncell[1]=web_dual[vertex+1];
                  pvertex -> pncell[2]=web_dual[vertex+LENGTH_FIRST_LATTICE+1];
               }

            }
            else {
               locerror("coord_gen", "The vertex does not exist");
            }

         }
      }
   }

   column_lim=LENGTH_FIRST_LATTICE-1;
   for (column=3 ; column<column_lim; column++){
      if (web[column] && web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+column]){
         if (column%2){
            row=1;
            vertex=vertex_from_coord(row, column);
            pvertex=web[vertex];
            pvertex -> x = -1.5 + 1.5*column;
            pvertex -> y = 0;

            pvertex -> pneighb[0]=web[vertex+1];
            pvertex -> pneighb[1]=web[vertex+LENGTH_FIRST_LATTICE-1];
            pvertex -> pneighb[2]=web[vertex-1];

            pvertex -> pncell[0]=web_dual[vertex+LENGTH_FIRST_LATTICE];	
            pvertex -> pncell[1]=web_dual[vertex];
            pvertex -> pncell[2]=web_dual[0];


            row=2*LENGTH_FIRST_LATTICE-1;
            vertex=vertex_from_coord(row, column);
            pvertex=web[vertex];

            pvertex -> x = -1.5 + 1.5*column;
            pvertex -> y = (row-1)*sqrt(3)/2;

            pvertex -> pneighb[0]=web[vertex+1];
            pvertex -> pneighb[1]=web[vertex-1];
            pvertex -> pneighb[2]=web[vertex-LENGTH_FIRST_LATTICE-1];

            pvertex -> pncell[0]=web_dual[0];	
            pvertex -> pncell[1]=web_dual[vertex];
            pvertex -> pncell[2]=web_dual[vertex-LENGTH_FIRST_LATTICE];
         }

         if ((column+1)%2){
            row=1;
            vertex=vertex_from_coord(row, column);
            pvertex=web[vertex];

            pvertex -> x = -2 + 1.5*column;
            pvertex -> y = 0;

            pvertex -> pneighb[0]=web[vertex-1];
            pvertex -> pneighb[1]=web[vertex+1];
            pvertex -> pneighb[2]=web[vertex+LENGTH_FIRST_LATTICE-1];

            pvertex -> pncell[0]=web_dual[0];	
            pvertex -> pncell[1]=web_dual[vertex+1];
            pvertex -> pncell[2]=web_dual[vertex + LENGTH_FIRST_LATTICE-1];


            row=2*LENGTH_FIRST_LATTICE-1;
            vertex=vertex_from_coord(row, column);
            pvertex=web[vertex];


            pvertex -> x = -2 + 1.5*column;
            pvertex -> y = (row-1)*sqrt(3)/2;

            pvertex -> pneighb[0]=web[vertex-1];
            pvertex -> pneighb[1]=web[vertex-LENGTH_FIRST_LATTICE-1];
            pvertex -> pneighb[2]=web[vertex+1];

            pvertex -> pncell[0]=web_dual[vertex-LENGTH_FIRST_LATTICE-1];
            pvertex -> pncell[1]=web_dual[vertex+1];
            pvertex -> pncell[2]=web_dual[0];
         }
      }
      else {
         locerror("coord_gen", "The vertex does not exist");
      }
   }

   for (row=5; row< 2 * LENGTH_FIRST_LATTICE - 4; row++){

      if (row%2){
         column=1;
         vertex=vertex_from_coord(row, column);
         pvertex=web[vertex];

         if (pvertex){

            pvertex -> x = 0;
            pvertex -> y = (row-1)*sqrt(3)/2;

            pvertex -> pneighb[0]=web[vertex+1];
            pvertex -> pneighb[1]=web[vertex+2*LENGTH_FIRST_LATTICE];
            pvertex -> pneighb[2]=web[vertex-2*LENGTH_FIRST_LATTICE];

            pvertex -> pncell[0]=web_dual[vertex+LENGTH_FIRST_LATTICE];
            pvertex -> pncell[1]=web_dual[0];
            pvertex -> pncell[2]=web_dual[vertex-LENGTH_FIRST_LATTICE];
         }
         else {
            locerror("coord_gen", "The vertex does not exist");
         }

         column=LENGTH_FIRST_LATTICE;
         vertex=vertex_from_coord(row, column);
         pvertex=web[vertex];

         if (pvertex){

            pvertex -> x = -2 + 1.5*LENGTH_FIRST_LATTICE;
            pvertex -> y = (row-1)*sqrt(3)/2;

            pvertex -> pneighb[0]=web[vertex-1];
            pvertex -> pneighb[1]=web[vertex-2*LENGTH_FIRST_LATTICE];
            pvertex -> pneighb[2]=web[vertex+2*LENGTH_FIRST_LATTICE];

            pvertex -> pncell[0]=web_dual[vertex-LENGTH_FIRST_LATTICE-1];
            pvertex -> pncell[1]=web_dual[0];
            pvertex -> pncell[2]=web_dual[vertex+LENGTH_FIRST_LATTICE-1];
         }
         else {
            locerror("coord_gen", "The vertex does not exist");
         }
      }
   }

   /* Corner generation */
   vertex=2;

   pvertex=web[vertex];

   pvertex -> x = 1;
   pvertex -> y = 0;

   pvertex -> pneighb[0]=web[2*LENGTH_FIRST_LATTICE+1];
   pvertex -> pneighb[1]=web[3];
   pvertex -> pneighb[2]=web[LENGTH_FIRST_LATTICE+1];

   pvertex -> pncell[0]=web_dual[0];
   pvertex -> pncell[1]=web_dual[3];
   pvertex -> pncell[2]=web_dual[LENGTH_FIRST_LATTICE + 1];


   vertex=2*LENGTH_FIRST_LATTICE+1;
   pvertex=web[vertex];

   pvertex -> x = 0;
   pvertex -> y = sqrt(3);

   pvertex -> pneighb[0]=web[2*LENGTH_FIRST_LATTICE+2];
   pvertex -> pneighb[1]=web[4*LENGTH_FIRST_LATTICE+1];
   pvertex -> pneighb[2]=web[2];

   pvertex -> pncell[0]=web_dual[3*LENGTH_FIRST_LATTICE+1];
   pvertex -> pncell[1]=web_dual[0];
   pvertex -> pncell[2]=web_dual[LENGTH_FIRST_LATTICE+1];


   vertex=LENGTH_FIRST_LATTICE-1;
   pvertex=web[vertex];

   pvertex -> x = -3 + 1.5*LENGTH_FIRST_LATTICE;
   pvertex -> y = 0;

   pvertex -> pneighb[0]=web[3*LENGTH_FIRST_LATTICE];
   pvertex -> pneighb[1]=web[2*LENGTH_FIRST_LATTICE-2];
   pvertex -> pneighb[2]=web[LENGTH_FIRST_LATTICE-2];

   pvertex -> pncell[0]=web_dual[2*LENGTH_FIRST_LATTICE-1];
   pvertex -> pncell[1]=web_dual[LENGTH_FIRST_LATTICE-1];
   pvertex -> pncell[2]=web_dual[0];


   vertex=3*LENGTH_FIRST_LATTICE;
   pvertex=web[vertex];

   pvertex -> x = -2 + 1.5*LENGTH_FIRST_LATTICE;
   pvertex -> y = sqrt(3);

   pvertex -> pneighb[0]=web[3*LENGTH_FIRST_LATTICE-1];
   pvertex -> pneighb[1]=web[LENGTH_FIRST_LATTICE-1];
   pvertex -> pneighb[2]=web[5*LENGTH_FIRST_LATTICE];

   pvertex -> pncell[0]=web_dual[2*LENGTH_FIRST_LATTICE-1];
   pvertex -> pncell[1]=web_dual[0];
   pvertex -> pncell[2]=web_dual[4*LENGTH_FIRST_LATTICE-1];


   vertex=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-4)+1;
   pvertex=web[vertex];

   pvertex -> x = 0;
   pvertex -> y = (2*LENGTH_FIRST_LATTICE-4)*sqrt(3)/2;

   pvertex -> pneighb[0]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-4)+2];
   pvertex -> pneighb[1]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+2];
   pvertex -> pneighb[2]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-6)+1];

   pvertex -> pncell[0]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)+1];
   pvertex -> pncell[1]=web_dual[0];
   pvertex -> pncell[2]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-5)+1];


   vertex=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+2;
   pvertex=web[vertex];

   pvertex -> x = 1;
   pvertex -> y = (2*LENGTH_FIRST_LATTICE-2)*sqrt(3)/2;

   pvertex -> pneighb[0]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-4)+1];
   pvertex -> pneighb[1]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)+1];
   pvertex -> pneighb[2]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+3];

   pvertex -> pncell[0]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)+1];
   pvertex -> pncell[1]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+3];
   pvertex -> pncell[2]=web_dual[0];


   vertex=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3);
   pvertex=web[vertex];

   pvertex -> x = -2 + 1.5*LENGTH_FIRST_LATTICE;
   pvertex -> y = (2*LENGTH_FIRST_LATTICE-4)*sqrt(3)/2;

   pvertex -> pneighb[0]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)-1];
   pvertex -> pneighb[1]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-5)];
   pvertex -> pneighb[2]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-1)-1];

   pvertex -> pncell[0]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-4)-1];
   pvertex -> pncell[1]=web_dual[0];
   pvertex -> pncell[2]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)-1];


   vertex=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-1)-1;
   pvertex=web[vertex];

   pvertex -> x = -3 + 1.5*LENGTH_FIRST_LATTICE;
   pvertex -> y = (2*LENGTH_FIRST_LATTICE-2)*sqrt(3)/2;

   pvertex -> pneighb[0]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)];
   pvertex -> pneighb[1]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-1)-2];
   pvertex -> pneighb[2]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)-2];

   pvertex -> pncell[0]=web_dual[0];
   pvertex -> pncell[1]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-1)-1];
   pvertex -> pncell[2]=web_dual[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)-1];

   /* List of vertices in cells */

   /* Hexagons */

   for (row=2 ; row<2*LENGTH_FIRST_LATTICE-1; row++){

      if (row%2){
         column_lim=LENGTH_FIRST_LATTICE;
         for (column=3 ; column<column_lim; column++){
            if (column%2){

               cell=vertex_from_coord(row, column);
               pcell=web_dual[cell];
               if (pcell){
                  pcell -> nb_vertices = 6;
                  pcell -> vertexlist = allocvertexlist(6);
                  pcell -> vertexlist[0]=web[cell];
                  pcell -> vertexlist[1]=web[cell+LENGTH_FIRST_LATTICE-1];
                  pcell -> vertexlist[2]=web[cell+LENGTH_FIRST_LATTICE-2];
                  pcell -> vertexlist[3]=web[cell-1];
                  pcell -> vertexlist[4]=web[cell-LENGTH_FIRST_LATTICE-2];
                  pcell -> vertexlist[5]=web[cell-LENGTH_FIRST_LATTICE-1];
               }
               else {
                  locerror("coord_gen", "The cell does not exist");
               }
            }

         }
      }

      if ((row+1)%2){
         column_lim=LENGTH_FIRST_LATTICE - 2;
         for (column=3 ; column<column_lim; column++){

            if (column%2){

               cell=vertex_from_coord(row, column);
               pcell=web_dual[cell];
               if (pcell){
                  pcell -> nb_vertices = 6;
                  pcell -> vertexlist = allocvertexlist(6);
                  pcell -> vertexlist[0]=web[cell];
                  pcell -> vertexlist[1]=web[cell+LENGTH_FIRST_LATTICE+1];
                  pcell -> vertexlist[2]=web[cell+LENGTH_FIRST_LATTICE];
                  pcell -> vertexlist[3]=web[cell-1];
                  pcell -> vertexlist[4]=web[cell-LENGTH_FIRST_LATTICE];
                  pcell -> vertexlist[5]=web[cell-LENGTH_FIRST_LATTICE+1];
               }
               else {
                  locerror("coord_gen", "The cell does not exist");
               }
            }
         }
      }
   }

   for (row=4 ; row<2*LENGTH_FIRST_LATTICE-3; row++){

      if ((row+1)%2){
         column=1;
         cell=vertex_from_coord(row, column);
         pcell=web_dual[cell];

         if (pcell){
            pcell -> nb_vertices = 5;
            pcell -> vertexlist = allocvertexlist(5);
            pcell -> vertexlist[0]=web[cell];
            pcell -> vertexlist[1]=web[cell+LENGTH_FIRST_LATTICE+1];
            pcell -> vertexlist[2]=web[cell+LENGTH_FIRST_LATTICE];
            pcell -> vertexlist[3]=web[cell-LENGTH_FIRST_LATTICE];
            pcell -> vertexlist[4]=web[cell-LENGTH_FIRST_LATTICE+1];
         }
         else {
            locerror("coord_gen", "The cell does not exist");
         }


         column=LENGTH_FIRST_LATTICE-2;
         cell=vertex_from_coord(row, column);
         pcell=web_dual[cell+1];

         if (pcell){
            pcell -> nb_vertices = 5;
            pcell -> vertexlist = allocvertexlist(5);
            pcell -> vertexlist[0]=web[cell];
            pcell -> vertexlist[1]=web[cell-LENGTH_FIRST_LATTICE+1];
            pcell -> vertexlist[2]=web[cell-LENGTH_FIRST_LATTICE+2];
            pcell -> vertexlist[3]=web[cell+LENGTH_FIRST_LATTICE+2];
            pcell -> vertexlist[4]=web[cell+LENGTH_FIRST_LATTICE+1];
         }
         else {
            locerror("coord_gen", "The cell does not exist");
         }

      }
   }

   column_lim=LENGTH_FIRST_LATTICE-1;
   for (column=2 ; column<column_lim; column++){
      if ((column+1)%2){
         row=1;
         cell=vertex_from_coord(row, column);
         pcell=web_dual[cell+1];
         if (pcell){
            pcell -> nb_vertices = 4;
            pcell -> vertexlist = allocvertexlist(4);
            pcell -> vertexlist[0]=web[cell];
            pcell -> vertexlist[1]=web[cell+1];
            pcell -> vertexlist[2]=web[cell+LENGTH_FIRST_LATTICE];
            pcell -> vertexlist[3]=web[cell+LENGTH_FIRST_LATTICE-1];
         }
         else {
            locerror("coord_gen", "The cell+1 does not exist");
         }


         row=2*LENGTH_FIRST_LATTICE-1;
         cell=vertex_from_coord(row, column);
         pcell=web_dual[cell+1];
         if (pcell){
            pcell -> nb_vertices = 4;
            pcell -> vertexlist = allocvertexlist(4);
            pcell -> vertexlist[0]=web[cell];
            pcell -> vertexlist[1]=web[cell-LENGTH_FIRST_LATTICE-1];
            pcell -> vertexlist[2]=web[cell-LENGTH_FIRST_LATTICE];
            pcell -> vertexlist[3]=web[cell+1];
         }
         else {
            locerror("coord_gen", "The cell+1 does not exist");
         }
      }
   }


   pcell=web_dual[LENGTH_FIRST_LATTICE+1];
   if (pcell){
      pcell -> nb_vertices = 4;
      pcell -> vertexlist = allocvertexlist(4);
      pcell -> vertexlist[0]=web[LENGTH_FIRST_LATTICE+1];
      pcell -> vertexlist[1]=web[2*LENGTH_FIRST_LATTICE+2];
      pcell -> vertexlist[2]=web[2*LENGTH_FIRST_LATTICE+1];
      pcell -> vertexlist[3]=web[2];
   }
   else {
      locerror("coord_gen", "The LENGTH_FIRST_LATTICE+1 does not exist");
   }

   pcell=web_dual[2*LENGTH_FIRST_LATTICE-1];
   if (pcell){
      pcell -> nb_vertices = 4;
      pcell -> vertexlist = allocvertexlist(4);
      pcell -> vertexlist[0]=web[3*LENGTH_FIRST_LATTICE];
      pcell -> vertexlist[1]=web[3*LENGTH_FIRST_LATTICE-1];
      pcell -> vertexlist[2]=web[2*LENGTH_FIRST_LATTICE-2];
      pcell -> vertexlist[3]=web[LENGTH_FIRST_LATTICE-1];
   }
   else {
      locerror("coord_gen", "The LENGTH_FIRST_LATTICE+1 does not exist");
   }

   cell=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)+1;
   pcell=web_dual[cell];
   if (pcell){
      pcell -> nb_vertices = 4;
      pcell -> vertexlist = allocvertexlist(4);
      pcell -> vertexlist[0]=web[cell];
      pcell -> vertexlist[1]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)+2];
      pcell -> vertexlist[2]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-4)+1];
      pcell -> vertexlist[3]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-4)+2];
   }
   else {
      locerror("coord_gen", "The LENGTH_FIRST_LATTICE+1 does not exist");
   }

   cell=LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)-1;
   pcell=web_dual[cell];
   if (pcell){
      pcell -> nb_vertices = 4;
      pcell -> vertexlist = allocvertexlist(4);
      pcell -> vertexlist[0]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-1)-1];
      pcell -> vertexlist[1]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-2)-2];
      pcell -> vertexlist[2]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)-1];
      pcell -> vertexlist[3]=web[LENGTH_FIRST_LATTICE*(2*LENGTH_FIRST_LATTICE-3)];
   }
   else {
      locerror("coord_gen", "The LENGTH_FIRST_LATTICE+1 does not exist");
   }


   /* Initialisation of the outside cell */

   pcell=web_dual[0];
   pcell -> nb_vertices = 4*(LENGTH_FIRST_LATTICE-2);
   pcell -> vertexlist = allocvertexlist(4*(LENGTH_FIRST_LATTICE-2));


   for (row=3 ; row<2*LENGTH_FIRST_LATTICE-2; row++){

      if (row%2){
         column=1;
         vertex=vertex_from_coord(row, column);
         pvertex=web[vertex];
         pcell -> vertexlist[(row-1)/2-1]=pvertex;

         column=LENGTH_FIRST_LATTICE;
         vertex=vertex_from_coord(row, column);
         pvertex=web[vertex];
         pcell -> vertexlist[3*LENGTH_FIRST_LATTICE-(row-1)/2-6]=pvertex;
      }
   }


   column_lim=LENGTH_FIRST_LATTICE;
   for (column=2 ; column<column_lim; column++){
      row=2*LENGTH_FIRST_LATTICE-1;
      vertex=vertex_from_coord(row, column);
      pvertex=web[vertex];
      pcell -> vertexlist[LENGTH_FIRST_LATTICE+column-4]=pvertex;

      row=1;
      vertex=vertex_from_coord(row, column);
      pvertex=web[vertex];
      pcell -> vertexlist[4*LENGTH_FIRST_LATTICE-column-7]=pvertex;
   }

   pcentercell=web_dual[vertex_from_coord(LENGTH_FIRST_LATTICE, LENGTH_FIRST_LATTICE/2)];
   if (!pcentercell){
      pcentercell=web_dual[vertex_from_coord(LENGTH_FIRST_LATTICE, LENGTH_FIRST_LATTICE/2) - 1];
      if (!pcentercell){
         locerror(" coord_gen ", " The central cell does not exist");
      }
   }
}

/* If a pointer to a cell or a vertex is void, the array is shrinked so that at the end the number of pointers is equal to the number of cells of vertices.
*/
   void
shrink_datafirst()
{
   int ivertex, idefvert, icell, idefcell;
   Node **webtemp;
   Dual **web_dualtemp;
   /* Node list shrinking */

   nb_vertex_tot=0;

   for (ivertex=0 ; ivertex < last_first_lattice ; ivertex++){
      if (web[ivertex]){
         nb_vertex_tot++;
      }
   }
   webtemp=createnodevector(nb_vertex_tot);

   //printf("Total number of vertices : %i\n", nb_vertex_tot);

   idefvert=0;

   for (ivertex=0 ; ivertex < last_first_lattice ; ivertex++){
      if (web[ivertex]){
         webtemp[idefvert]=web[ivertex];
         idefvert++;
      }
   }
   free(web);
   web=webtemp;

   /* Cell list shrinking */    
   nb_cell_tot=0;

   for (icell=0 ; icell<last_first_lattice ; icell++){
      if (web_dual[icell]){
         nb_cell_tot++;
      }
   }
   web_dualtemp=createdualvector(nb_cell_tot);

   //printf("Total number of cells : %i\n", nb_cell_tot);fflush(stdout);

   idefcell=0;

   for (icell=0 ; icell<last_first_lattice ; icell++){
      if (web_dual[icell]){
         web_dualtemp[idefcell]=web_dual[icell];
         idefcell++;
      }
   }
   free(web_dual);
   web_dual=web_dualtemp;
}

/* Shrink web[] and web_dual[] when a cell is being killed
*/
   void
shrink_data()
{
   int ivertex, icell;
   int nb_vertex_def=0;
   int nb_cell_def=0;
   int idefvert=0;
   int idefcell=0;
   Node **webtemp;
   Dual **web_dualtemp;

   /* Node list shrinking */
   for (ivertex=0 ; ivertex < nb_vertex_tot ; ivertex++){
      if (web[ivertex]){
         nb_vertex_def++;
      }
   }
   webtemp=createnodevector(nb_vertex_def);

   for (ivertex=0 ; ivertex < nb_vertex_tot ; ivertex++){
      if (web[ivertex]){
         webtemp[idefvert]=web[ivertex];
         idefvert++;
      }
   }
   nb_vertex_tot=nb_vertex_def;
   //printf("Total number of vertices : %d\n",nb_vertex_tot);

   free(web);
   web=webtemp;

   /* Cell list shrinking */    
   for (icell=0 ; icell<nb_cell_tot ; icell++){
      if (web_dual[icell]){
         nb_cell_def++;
      }
   }
   web_dualtemp=createdualvector(nb_cell_def);

   for (icell=0 ; icell<nb_cell_tot ; icell++){
      if (web_dual[icell]){
         web_dualtemp[idefcell]=web_dual[icell];
         idefcell++;
      }
   }
   nb_cell_tot=nb_cell_def;
   //printf("Total number of cells : %d\n",nb_cell_tot);fflush(stdout);

   free(web_dual);
   web_dual=web_dualtemp;
}

/* Computes the exterior cell composed of the vertices on the boundary
*/
   void
cell0gen(Node *pfirstvertex)
{
   int icell, vertcount;
   Node *pvertex;
   Dual *pcell0;

   pvertex=NULL;
   vertcount=1;
   pcell0=web_dual[0];

   icell=indexcell(pfirstvertex, pcell0);
   pvertex=pfirstvertex -> pneighb[icell];
   while(pvertex!=pfirstvertex){
      vertcount++;
      icell=indexcell(pvertex, pcell0);
      pvertex=pvertex -> pneighb[icell];
   }

   free(pcell0 -> vertexlist);
   pcell0 -> vertexlist = createnodevector(vertcount);
   pcell0 -> nb_vertices = vertcount;

   pvertex=NULL;
   vertcount=0;
   pcell0 -> vertexlist[0]=pfirstvertex;
   icell=indexcell(pfirstvertex, pcell0);
   pvertex=pfirstvertex -> pneighb[icell];
   while(pvertex!=pfirstvertex){
      vertcount++;
      icell=indexcell(pvertex, web_dual[0]);
      pvertex=pvertex -> pneighb[icell];
   }
}




/* indexcell return the index of cell for a given vertex and cell */ 
   int
indexcell(Node *pvertex, Dual *pcell)
{
   int i;

   for (i=0; i<3 ; i++ ){
      if ((pvertex -> pncell[i])==pcell) return i;
   }
   locerror("indexcell", "the cell does not belong to the vertex");
   return -1;
}


/* Area computation */

   void
cell_update()
{
   int i, icell, ivert;
   Dual *pcell;
   Node *pvert;

   for (ivert=0;ivert< nb_vertex_tot;ivert++){
      pvert=web[ivert];
      for (i=0 ; i<3 ; i++){
         pvert->dist_vert[i]=dist(pvert, pvert->pneighb[i]);
      }
   }
   for (icell=0;icell<nb_cell_tot;icell++){
      pcell=web_dual[icell];
      pcell -> area=area(pcell);

      pcell -> centroid=centroid2(pcell);
   }
}

/* Computes the centroid of a cell (considering vertices)
*/
   vector3
centroid(Dual *pcell)
{
   int i, nbvertices;
   Node *pvert;
   vector3 center;

   center.x=0;
   center.y=0;
   center.z=0;

   nbvertices=pcell->nb_vertices;

   for (i=0; i<nbvertices; i++){
      pvert=pcell->vertexlist[i];

      center.x += pvert->x;
      center.y += pvert->y;
      center.z += pvert->z;
   }

   center.x /= nbvertices;
   center.y /= nbvertices;
   center.z /= nbvertices;

   return center;
}

/* Computes the centroid of a cell (considering edges)
*/
   vector3
centroid2(Dual *pcell)
{
   register unsigned int i, nbvertices, index_cell;
   register double perim, distance;
   register Node *pvert, *pvertnext;
   register vector3 center;

   center.x=0;
   center.y=0;
   center.z=0;
   perim=0;

   nbvertices=pcell->nb_vertices;

   for (i=0; i<nbvertices-1; i++){
      pvert=pcell->vertexlist[i];
      pvertnext=pcell->vertexlist[i+1];

      index_cell=indexcell(pvert, pcell);

      distance=pvert->dist_vert[index_cell];

      center.x += distance*(pvert->x+pvertnext->x);
      center.y += distance*(pvert->y+pvertnext->y);
      center.z += distance*(pvert->z+pvertnext->z);
      perim += distance;
   }

   pvert=pcell->vertexlist[nbvertices-1];
   pvertnext=pcell->vertexlist[0];

   index_cell=indexcell(pvert, pcell);

   distance=pvert->dist_vert[index_cell];

   center.x += distance*(pvert->x+pvertnext->x);
   center.y += distance*(pvert->y+pvertnext->y);
   center.z += distance*(pvert->z+pvertnext->z);

   perim += distance;

   center.x /= perim*2;
   center.y /= perim*2;
   center.z /= perim*2;/*CHECK_HERE not sure about this... maybe *3?*/

   return center;
}

/* distance between 2 vertices */
   double
dist(Node *pvertex1, Node *pvertex2)
{
   register double distx, disty, distz;

   distx=pvertex1->x - pvertex2->x;
   disty=pvertex1->y - pvertex2->y;
   distz=pvertex1->z - pvertex2->z;
   return sqrt(distx*distx + disty*disty + distz*distz + LEN_CUTOFF);
}


   int
indexvertex(Node *pvertex, Dual *pcell)
{
   int nbvertices, n;

   nbvertices=pcell->nb_vertices;

   for (n=0; n < nbvertices ; n++){
      if ((pcell->vertexlist[n])==pvertex){
         return n;
      }
   }
   locerror("indexvertex", "The vertex does not belong to the cell");
   return -1;
}

   int
indexvertexinweb(Node *pvertex)
{
   int n;

   for (n=0; n < nb_vertex_tot ; n++){
      if (web[n]==pvertex){
         return n;
      }
   }
   locerror("indexvertexinweb", "The vertex does not belong to the web!!");
   return 0;
}

   int
indexcellinwebdual(Dual *pcell)
{
   int n;

   for (n=0; n < nb_cell_tot ; n++){
      if (web_dual[n]==pcell){
         return n;
      }
   }
   locerror("indexcellinwebdual", "The cell does not belong to the web_dual!!");
   return 0;
}

/* Insert a new vertex in vertexlist */
   void
insvertex(Dual *pcell, Node *pnewvert, Node *pnextvert) /* The new vertex is at the first position */
{
   int nbvertices, n, ivert;
   Node **pvertexlisttemp;

   if ((ivert=indexvertex(pnextvert, pcell)) < 0){
      printf("In insvertex ivert < 0\n");
      exit(1);
   };

   nbvertices=pcell->nb_vertices;
   pvertexlisttemp=createnodevector(nbvertices+1);
   pcell->nb_vertices =nbvertices+1;

   pvertexlisttemp[0]=pnewvert;
   for (n=0 ; n < nbvertices ; n++){
      pvertexlisttemp[n+1]=pcell->vertexlist[(ivert+n)%nbvertices];
   }

   free(pcell->vertexlist);
   pcell->vertexlist = pvertexlisttemp;
}

   void
delvertex(Dual *pcell, Node *pdelvert)
{
   int ivert, n, nbvertices;
   Node **pvertexlisttemp;

   if ((ivert=indexvertex(pdelvert, pcell)) < 0){
      printf("In delvertex ivert < 0\n");
      exit(1);
   };


   nbvertices=pcell->nb_vertices;
   pvertexlisttemp=createnodevector(nbvertices-1);
   pcell->nb_vertices =nbvertices-1;

   for (n=0 ; n < nbvertices-1 ; n++){
      pvertexlisttemp[n]=pcell->vertexlist[(ivert+n+1)%nbvertices];
   }

   free(pcell->vertexlist);
   pcell->vertexlist = pvertexlisttemp;
}


   int
kill_cell_p(Dual *pcell)
{
   int i, j, ix;
   Node *ptmpvert,*ptttvert;
   Node *cv[2];
   Dual *pcelltmp,*pcell0, *pcellp;
   double pres, presdif;


   pcell0 = pcellp = web_dual[0];

   // find neighboring cell with highest/lowest pressure
   pres = pressure(pcell);
   printf("%f\n",pres);
   presdif = 0.;
   for (i=0;i<pcell->nb_vertices;i++) {
      ptmpvert = pcell->vertexlist[i];
      for (j=0;j<3;j++){
         pcelltmp = ptmpvert->pncell[j];
         if (pcelltmp != pcell && pcelltmp != pcell0){
            printf("%i %i %f\n",i,j,pressure(pcelltmp));
            if (fabs(pressure(pcelltmp)-pres) > presdif) {
               pcellp  = pcelltmp;
               presdif = fabs(pressure(pcelltmp)-pres);
            }
         }
      }
   }
   printf("%f %f\n",pres, pressure(pcellp));

   if (pcell->nb_vertices == 3 && pcellp->nb_vertices == 3) {
      printf("pcell and pcellp have both 3 vertices\n");
      exit(1);
   }

   // find verticies commen two pcell and pcellp
   printf("pcell  has %i vertices\n",pcell->nb_vertices);
   printf("pcellp has %i vertices\n",pcellp->nb_vertices);
   ix = 0;
   for (i=0;i<pcell->nb_vertices;i++) {
      ptmpvert = pcell->vertexlist[i];
      for (j=0;j<pcellp->nb_vertices;j++) {
         ptttvert = pcellp->vertexlist[j];
         if (ptttvert == ptmpvert) {
            if (ix >1) {
               printf("pcell and pcellc hav more than 2 vertices in commen\n");
               exit(1);
            }
            cv[ix] = ptttvert;
            ix++;
         }
      }
   }

   // check neighboring cells
   for (j=0;j<3;j++){
      ptmpvert = cv[0];
      pcelltmp = ptmpvert->pncell[j];
      if (pcelltmp != pcell && pcelltmp != pcell0 && pcelltmp->nb_vertices == 3){
         printf("cv[1] neighb cell has both 3 vertices\n");
         exit(1);
      }
      ptmpvert = cv[1];
      pcelltmp = ptmpvert->pncell[j];
      if (pcelltmp != pcell && pcelltmp != pcell0 && pcelltmp->nb_vertices == 3){
         printf("cv[2] neighb cell has both 3 vertices\n");
         exit(1);
      }
   }

   for (j=0;j<3;j++) {
      ptmpvert = cv[0]->pneighb[j];
   }



   exit(1);
}


   int
find_smallc(Dual *pcell)
{
   int nvert, index=0;
   double areamin=100000;

   for (nvert=0;nvert<pcell->nb_vertices;nvert++){
      if (pcell->area<areamin){
         areamin=pcell->area;
         index=nvert;
      }
   }
   return index;
}

   int
kill_cell_area(Dual *pcell)
{
   int index;
   vector3 center;

   index=find_smallc(pcell);
   center=centroid(pcell);
   kill_cell(pcell,index);
   return 0;
}

/* Kill a cell */
   int
kill_cell(Dual *pcell, int ivertex)
{
   int i, icell, icell1, icell2, icella, icellb, icellc, icelld, ivertexbig;
   int nbvertices1, nbvertices2, clone;
   //char filename[256];
   Node *pvert;
   Node *pvertex1, *pvertex2, *pvertexa, *pvertexb, *pvertexc, *pvertexd;
   Node **vertexlisttemp;
   Dual *pcellbig,*pcelltmp,*pcell0;
   double pres, presdif;


   pcell0 = web_dual[0];
   clone = 0;
   if (pcell->marker.div_rate == 1) {
      clone = 1;
   }

   // find neighboring cell with highest/lowest pressure
   pres = pressure(pcell);
   //printf("%f\n",pres);
   presdif = 0.;
   /*for (i=0;i<pcell->nb_vertices;i++) {
     ptmpvert = pcell->vertexlist[i];
     j=indexcell(ptmpvert, pcell);
     pcelltmp=ptmpvert->pncell[(j+2)%3];

     if (pcelltmp != pcell && pcelltmp != pcell0){
   //printf("%i %i %f\n",i,j,pressure(pcelltmp));
   if (fabs(pressure(pcelltmp)-pres) > presdif) {
   ivertex = j;
   }
   }
   }*/

   pvertex1=pcell->vertexlist[ivertex];
   icell1=indexcell(pvertex1, pcell);
   pcellbig=pvertex1->pncell[(icell1+2)%3];


   if ((ivertexbig=indexvertex(pvertex1, pcellbig)) < 0){
      printf("In kill_cell ivertexbig < 0\n");
      exit(1);
   };

   pvertexa=pvertex1->pneighb[(icell1+1)%3];
   icella=indexcell(pvertexa, pcell);
   pvertexb=pvertex1->pneighb[(icell1+2)%3];
   icellb=indexcell(pvertexb, pcellbig);
   pvertex2=pvertex1->pneighb[icell1];
   icell2=indexcell(pvertex2, pcell);

   //check if either vertex 1 or 2 consists of only three edges
   /*
      if (pcell->nb_vertices == 3) {
      printf("- cell has 3 edges -> not killed\n");
      for (i=0;i<3;i++){
      ptmpvert=pcell->vertexlist[i];
      for (j=0;j<3;j++){
      pcelltmp = ptmpvert->pncell[j];
      if (pcelltmp == pcell){
      }
      else if (pcelltmp->nb_vertices == 3){
      printf("- cell has %i edges -> not killed\n",pcelltmp->nb_vertices);
      return 0;
      }
      }

      }
      }
      */
   for (i=0;i<3;i++){
      pcelltmp=pvertex1->pncell[(icell1+i)%3];
      if (pcelltmp->nb_vertices <=3 && pcelltmp != pcell && pcelltmp != pcellbig) {
         printf("- cell has %i edges -> not killed\n",pcelltmp->nb_vertices);
         return 0;
      }
   }
   for (i=0;i<3;i++){
      pcelltmp=pvertex2->pncell[(icell1+i)%3];
      if (pcelltmp->nb_vertices <=3 && pcelltmp != pcell && pcelltmp != pcellbig) {
         printf("- cell has %i edges -> not killed\n",pcelltmp->nb_vertices);
         return 0;
      }
   }

   pvertexc=pvertex2->pneighb[(icell2+2)%3];
   icellc=indexcell(pvertexc, pcellbig);
   pvertexd=pvertex2->pneighb[icell2];
   icelld=indexcell(pvertexd, pcell);

   pvertexa->pneighb[icella]=pvertexb;
   pvertexb->pneighb[(icellb+1)%3]=pvertexa;

   pvertexc->pneighb[icellc]=pvertexd;
   pvertexd->pneighb[(icelld+1)%3]=pvertexc;

   nbvertices1=pcell->nb_vertices;
   nbvertices2=pcellbig->nb_vertices;

   for(i=0; i < nbvertices1-2; i++){
      pvert=pcell->vertexlist[(ivertex+2+i)%nbvertices1];
      icell=indexcell(pvert, pcell);
      pvert->pncell[icell]=pcellbig;
   }

   vertexlisttemp=createnodevector(nbvertices1+nbvertices2-4);

   for(i=0; i < nbvertices1-2; i++){
      vertexlisttemp[i]=pcell->vertexlist[(ivertex+2+i)%nbvertices1];
   }
   for(i=0; i < nbvertices2-2; i++){
      vertexlisttemp[nbvertices1-2+i]=pcellbig->vertexlist[(ivertexbig+1+i)%nbvertices2];
   }


   web_dual[indexcellinwebdual(pcell)]=NULL;
   freecell(pcell);

   web[indexvertexinweb(pvertex1)]=NULL;
   web[indexvertexinweb(pvertex2)]=NULL;

   free(pcellbig->vertexlist);

   pcellbig->vertexlist=vertexlisttemp;
   pcellbig->nb_vertices=nbvertices1+nbvertices2-4;

   pcellbig->area_soll = 1.; //(1+cell_area(pcellbig))*0.5;
   //printf("kill area %f\n",pcellbig->area_soll);

   delvertex(pvertex1->pncell[(icell1+1)%3], pvertex1);
   delvertex(pvertex2->pncell[(icell2+2)%3], pvertex2);
   free(pvertex1);
   free(pvertex2);

   shrink_data();

   //snprintf(filename, 256, "./lattice.eps");
   //epsfile(filename,0);

   if (clone == 1) {
      nb_clone_tot--;
   }

   return 1;
}



/* Division functions */

   void
division_cell(Dual *pcell, int ivertex)
{
   int icell;
   Node *pvertex;

   pvertex=pcell->vertexlist[ivertex];
   icell=indexcell(pvertex, pcell);
   division_eq(pvertex, icell);
}


   void
division_eq(Node *pvertex, int icell)
{
   int n, icelltemp, ivert0cell0, nbvertices, last_vert, icell1, icell2, icell3;
   int ivertex, j,count,k;
   Node *pvert, *pvertprev, *pvertex1, *pvertex2, *pvertex3, *pnewvertex1, *pnewvertex2, **pvertexlisttemp;
   Dual *pcell, *pnewcell, *pncell0m1, *pncell2m1;
   double a, a0, x, y, z;

   nb_vertex_tot += 2;
   nb_cell_tot++;

   web = reallocnodevector(web, nb_vertex_tot);
   web[nb_vertex_tot - 2]=allocnode();
   web[nb_vertex_tot - 1]=allocnode();
   pnewvertex1=web[nb_vertex_tot - 2];
   pnewvertex2=web[nb_vertex_tot - 1];

   web_dual= reallocdualvector(web_dual, nb_cell_tot);
   web_dual[nb_cell_tot - 1]=allocdual();
   pnewcell=web_dual[nb_cell_tot - 1];

   pcell=pvertex->pncell[icell];


   if ((ivert0cell0=indexvertex(pvertex, pcell)) < 0){
      printf("In division ivert0cell0 < 0\n");
      exit(1);
   };

   // choose cleavage plane such that areas of daughter cells are close to 0.5
   a0 = area(pcell);
   nbvertices=pcell->nb_vertices;
   a = 0.;
   for (j=1;j <= nbvertices && a<a0/2; j++) {

      pvert    =pcell -> vertexlist[(ivert0cell0+1)%nbvertices];
      pvertprev=pcell -> vertexlist[ivert0cell0%nbvertices];
      x = 0.5*(pvertprev -> x + pvert -> x);
      y = 0.5*(pvertprev -> y + pvert -> y); 
      z = 0.5*(pvertprev -> z + pvert -> z);

      //currently calculated like a cross product, with fixed orientation.
      /*    a=(y + pvert -> y) * (x - pvert -> x); //this needs to be changed.
            for (ivertex=1; ivertex < j ; ivertex++){
            pvertprev=pvert;
            pvert=pcell -> vertexlist[(ivert0cell0+ivertex+1)%nbvertices];
            a += (pvert -> y + pvertprev -> y) * (pvertprev -> x - pvert -> x);
            }
            a +=(y + pvert -> y) * (pvert -> x - x);
            a *= 0.5;*/
      a = part_area(pcell,j,ivert0cell0,x,y,z);//kevin_tools.c partial area calc
   }
   last_vert = j-2;
   pvertex =pcell->vertexlist[(ivert0cell0 )%nbvertices];
   pvertex1=pcell->vertexlist[(ivert0cell0+1)%nbvertices];
   pvertex2=pcell->vertexlist[(ivert0cell0+last_vert)%nbvertices];
   pvertex3=pcell->vertexlist[(ivert0cell0+last_vert+1)%nbvertices];

   icell1=indexcell(pvertex1, pcell);
   icell2=indexcell(pvertex2, pcell);
   icell3=indexcell(pvertex3, pcell);

   pncell0m1=pvertex->pncell[(icell+2)%3];

   pncell2m1=pvertex2->pncell[(icell2+2)%3];

   pnewvertex1->x=(pvertex->x + pvertex1->x)/2;
   pnewvertex1->y=(pvertex->y + pvertex1->y)/2;
   pnewvertex1->z=(pvertex->z + pvertex1->z)/2;
   pnewvertex1->pneighb[0]=pvertex;
   pnewvertex1->pneighb[1]=pvertex1;
   pnewvertex1->pneighb[2]=pnewvertex2;
   pnewvertex1->pncell[0]=pncell0m1;
   pnewvertex1->pncell[1]=pnewcell;
   pnewvertex1->pncell[2]=pcell;

   pnewvertex2->x=(pvertex2->x + pvertex3->x)/2;
   pnewvertex2->y=(pvertex2->y + pvertex3->y)/2;
   pnewvertex2->z=(pvertex2->z + pvertex3->z)/2;
   pnewvertex2->pneighb[0]=pvertex2;
   pnewvertex2->pneighb[1]=pvertex3;
   pnewvertex2->pneighb[2]=pnewvertex1;
   pnewvertex2->pncell[0]=pncell2m1;
   pnewvertex2->pncell[1]=pcell;
   pnewvertex2->pncell[2]=pnewcell;


   pvertex->pneighb[icell]=pnewvertex1;

   pvertex1->pneighb[(icell1+1)%3]=pnewvertex1;

   pvertex2->pneighb[icell2]=pnewvertex2;

   pvertex3->pneighb[(icell3+1)%3]=pnewvertex2;

   for (n=1 ; n < last_vert +1 ; n++){
      pvert=pcell->vertexlist[(ivert0cell0+n)%nbvertices];
      icelltemp=indexcell(pvert, pcell);
      pvert->pncell[icelltemp]=pnewcell;
   }

   pnewcell->vertexlist=createnodevector(last_vert+2);
   pnewcell->nb_vertices=last_vert+2;
   pnewcell->marker=pcell->marker;


   pnewcell->vertexlist[0]=pnewvertex1;
   for (n=1 ; n < last_vert +1 ; n++){
      pnewcell->vertexlist[n]=pcell->vertexlist[(ivert0cell0+n)%nbvertices];
   }
   pnewcell->vertexlist[last_vert +1]=pnewvertex2;


   pvertexlisttemp=createnodevector(nbvertices-last_vert+2);
   pcell->nb_vertices=nbvertices-last_vert+2;

   pvertexlisttemp[0]=pnewvertex2;
   for (n=1 ; n < nbvertices-last_vert +1 ; n++){
      pvertexlisttemp[n]=pcell->vertexlist[(ivert0cell0+last_vert+n)%nbvertices];
   }
   pvertexlisttemp[nbvertices-last_vert+1]=pnewvertex1;

   free(pcell->vertexlist);
   pcell->vertexlist=pvertexlisttemp;

   insvertex(pncell0m1, pnewvertex1 ,pvertex);

   insvertex(pncell2m1, pnewvertex2 ,pvertex2);

   pnewcell->marker.div_rate = pcell->marker.div_rate;
   pnewcell->marker.tension_index = pcell->marker.tension_index;
   pnewcell->marker.border = pcell->marker.border;
   pnewcell->area_soll0 = pcell->area_soll0;
   pnewcell->area_soll = pnewcell->area_soll0;
   pcell->area_soll = pcell->area_soll0;

   //    pnewcell->area_soll=area(pnewcell);
   //pcell->area_soll    = 1.; //area(pcell);
   //pnewcell->area_soll = 1.;
   //printf("areas %f %f\n",area(pcell),area(pcell)/(area(pcell)+area(pnewcell)));
   count = 0;
   for(k=0;k<nb_vertex_tot;k++){
      if(web[k]->border==1){count++;}
   }
   printf("areas %f %f %d\n",area(pcell),area(pcell)/(area(pcell)+area(pnewcell)),count);


}

/* substitute the area() function from kevin_tools.h.  more general for 3d case
   double
   area(Dual *pcell)
   {
   int ivert, nbvert;
   double a;
   Node *pv,*pvprev;

   nbvert=pcell -> nb_vertices;
   pvprev=pcell -> vertexlist[nbvert-1];
   pv=pcell -> vertexlist[0];
   a=sqrt(pow(pvprev->x*pv->y-pvprev->y*pv->x,2)+pow(pvprev->z*pv->x-pvprev->x*pv->z,2)+pow(pvprev->y*pv->z-pvprev->z*pv->y,2));

   for (ivert=1; ivert < nbvert; ivert++){
   pvprev=pv;
   pv=pcell -> vertexlist[ivert];
   a += sqrt(pow(pvprev->x*pv->y-pvprev->y*pv->x,2)+pow(pvprev->z*pv->x-pvprev->x*pv->z,2)+pow(pvprev->y*pv->z-pvprev->z*pv->y,2));
   }
   return 0.5*a;
   }*/




/* Computes the initialization vertex to divide a cell along its main axis
*/
   int
vertex_axis(Dual *pcell)
{
   int i, ivertex;
   int nbvertices;
   double det, nextdet;
   vector3 eigenvect;
   vector3 vertexpos, nextvertexpos;
   vector3 center;
   Node *pnextvert;

   center=centroid(pcell);
   eigenvect=eigenvecmin(stress(pcell));

   nbvertices=pcell->nb_vertices;
   ivertex = (int)(gsl_rng_uniform(rng)*nbvertices);	/* To be sure ivertex is properly initialized */

   pnextvert=pcell -> vertexlist[0];
   nextvertexpos.x = pnextvert->x - center.x;
   nextvertexpos.y = pnextvert->y - center.y;
   nextvertexpos.z = pnextvert->z - center.z;
   nextdet = nextvertexpos.x*eigenvect.y - nextvertexpos.y*eigenvect.x;//CHECK_HERE this needs to be changed

   for (i=0; i<nbvertices; i++){
      vertexpos.x=nextvertexpos.x;
      vertexpos.y=nextvertexpos.y;
      vertexpos.z=nextvertexpos.z;

      pnextvert=pcell -> vertexlist[(i+1)%nbvertices];
      nextvertexpos.x=pnextvert->x-center.x;
      nextvertexpos.y=pnextvert->y-center.y;
      nextvertexpos.z=pnextvert->z-center.z;

      det=nextdet;
      nextdet=nextvertexpos.x*eigenvect.y - nextvertexpos.y*eigenvect.x;//CHECK_HERE this needs to be changed too
      if (det*nextdet<0){
         ivertex=i;
         break;
      }
   }

   return ivertex;
}


/* Effective division : it takes into account the probability of division and call division() */
   void
effdivision()
{
   int celltodivide; 	/* index of the cell to divide in web_dual */
   int nbvertices; 	/* # of vertices in the cell, used to determine the stating point of the division */
   int ivertex; 	/* index of the starting point of the division */
   static int itime_pressure;	/* Time counting for pressure killing process */

   Dual *pcell;


   if (gsl_rng_uniform(rng) < DIVISION_RATE){
      celltodivide = (int)(gsl_rng_uniform(rng) * (nb_cell_tot-1))+1;

      pcell=web_dual[celltodivide];
      nbvertices=pcell->nb_vertices;


      /* Checking point activation if CHECKED_DIVISION != 0 */
#if CHECKED_DIVISION
      if (!(pcell->marker.checked_division)){
#if DIVISION_RATE_CHANGEABLE
         if (pcell -> marker.div_rate==1 && gsl_rng_uniform(rng) < division_rates[1]){
            pcell->marker.checked_division=1;
         }
         else if (gsl_rng_uniform(rng) < division_rates[0]){
            pcell->marker.checked_division=1;
         }
         return;
#else

         pcell->marker.checked_division=1;
         return;
#endif /* DIVISION_RATE_CHANGEABLE */
      }
#endif /* CHECKED_DIVISION */



      /* Determination of the vertex near axis division :
       * mechanically driven if AXIS_DIVISION != 0
       * stochastically driven otherwise */
#if AXIS_DIVISION
      ivertex=vertex_axis(pcell);
#else
      ivertex = (int)(gsl_rng_uniform(rng)*nbvertices);
#endif /* AXIS_DIVISION */



#if DIVISION_RATE_CHANGEABLE
      if (pcell -> marker.div_rate==1  && gsl_rng_uniform(rng) < division_rates[1]){
         division_cell(pcell, ivertex);
         itime_pressure=1;
      }
      else if (gsl_rng_uniform(rng) < division_rates[0]){
         division_cell(pcell, ivertex);
         itime_pressure=1;
      }
#else
      division_cell(pcell, ivertex);
      itime_pressure=1;
#endif
   }

   /* Time based on stochastic behavior */
   itime++;

   itime_pressure++;

   /* Kill cells which are too much stressed */
#if PRESSURE_KILL
   if (!(itime_pressure%PRESSURE_KILL_INTER)){
      pressure_kill();
   }
#endif /* PRESSURE_KILL */
}

/* pressure_kill() :
 * kill cells whose pressure is too high
 */
   void
pressure_kill()
{
   int icell;	/* Loop variable */

   Dual *pcell;

   double pressure_cell;

   for (icell=1 ; icell<nb_cell_tot ; icell++){
      pcell=web_dual[icell];

      pressure_cell=pressure(pcell);

      if (pressure_cell > PRESSURE_THRESHOLD){
         kill_cell(pcell, 0);
         nb_cell_killed++;
         printf("KILLED cell %i pressure %f total killed %i\n",icell, pressure_cell, nb_cell_killed);
      }
   }
}



/* Implementation of the t1 process
 * It has to be noticed that the choice of the vertex c is not arbitrary 
 * and can induce a rotation on the lattice. This effect should be negligible though (it is in practice)
 */
   void
t1transform(Node *pvertex, int neighbor)
{
   int indexbneighbc, indexcneighbb, indexdneighba;
   Node *pvertexb, *pvertexc, *pvertexd, *pvertexe, *pvertexf; /* The explanation of the variables is on the figure in the documentation */
   Dual *pcellalpha, *pcellbeta, *pcellgamma, *pcelldelta;

   pvertexb=pvertex->pneighb[neighbor];

   pcellalpha=pvertex->pncell[neighbor];
   pcellbeta=pvertex->pncell[(neighbor+2)%3];
   indexbneighbc=indexcell(pvertexb, pcellalpha);
   pvertexc=pvertexb->pneighb[indexbneighbc];
   pvertexd=pvertex->pneighb[(neighbor+2)%3];

   pcellgamma=pvertex->pncell[(neighbor+1)%3];
   pcelldelta=pvertexb->pncell[(indexbneighbc+2)%3];

   indexcneighbb=indexcell(pvertexc, pcelldelta);
   indexdneighba=indexcell(pvertexd, pcellgamma);

   pvertexe=pvertex->pneighb[(neighbor+1)%3];
   pvertexf=pvertexb->pneighb[(indexbneighbc+2)%3];

   //if(pcellalpha->nb_vertices<=4 || pcellbeta->nb_vertices<=4){return;}

   /* Changes */

   pvertexd->pneighb[indexdneighba]=pvertexb;
   pvertexc->pneighb[indexcneighbb]=pvertex;

   insvertex(pcellgamma, pvertexb, pvertex);
   insvertex(pcelldelta, pvertex, pvertexb);

   delvertex(pcellalpha, pvertexb);
   delvertex(pcellbeta, pvertex);

   pvertex->pneighb[0]=pvertexc;
   pvertex->pneighb[1]=pvertexe;
   pvertex->pneighb[2]=pvertexb;
   pvertex->pncell[0]=pcellalpha;
   pvertex->pncell[1]=pcellgamma;
   pvertex->pncell[2]=pcelldelta;

   pvertexb->pneighb[0]=pvertexd;
   pvertexb->pneighb[1]=pvertexf;
   pvertexb->pneighb[2]=pvertex;
   pvertexb->pncell[0]=pcellbeta;
   pvertexb->pncell[1]=pcelldelta;
   pvertexb->pncell[2]=pcellgamma;

   it1process++;
}


/* Statistics on the edges of cells
*/
   void
stat_edges()
{
   int icell;
   Dual *pcell;

   for (icell=1;icell<nb_cell_tot;icell++){
      pcell=web_dual[icell];
      edges_dist[pcell->nb_vertices]++;
   }
}


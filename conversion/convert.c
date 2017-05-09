///////////////////////////////////////////////////////////////
//
//  Reads in a file with bond, vertex and cell structure in a precise format and 
//  writes out the configuration file ready to be read by GMO
//  
// File vb (v=contains vertex and bond structure)
//
// Total # of vertices
// Total # of bonds
// Total # of cells
// x_vertex y_vertex  (the coordinates are referred to the matlab image, i.e. origin is top left)
//    .        .
//    .        .
//    .        .
// ind_neighb_cell1 ind_neighb_cell2 x_vertex_1 y_vertex_1 x_vertex_2 y_vertex_2 
//    .        .        .                .          .          .          .
//    .        .        .                .          .          .          .
//    .        .        .                .          .          .          .
//
//
//  Alberto, 10/31/2006
//
//
/////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include "convert.h"
#include "locerror.h"
#include "save.h"
#include "kevin_tools.h"
#include "read.h"
#include "correction.h"


//Variables for Node and Dual structures

Dual **web_dual;
Node **web;
Bond **web_bond;

int nb_bond_tot,nb_cell_tot,nb_vertex_tot;


int converta(){
  char filename[255],file1[255],file2[255];
  int i,j,k,l,count,flag,chk;
  Dual *pncell[2],*pcell0,*pcell;
  Node *pvert,*pnvert[2],**newvertexlist;

  scanf("%s",filename);
  readweb(filename);

  //this next section is to correct the vertex order of cell0.
  //for whatever reason, that part is not correct.  however
  //note that this assumes all other cells have correct order
  orientation_correction();

  /*
  sprintf(file1,"%s.vertinfo",filename);
  sprintf(file2,"%s.cellinfo",filename);
  infoinput(file1,file2);
  */

  for(k=0;k<nb_cell_tot;k++){
    pcell=web_dual[k];
    pcell->marker.size=0;
    pcell->marker.div_rate=0;
    pcell->area_soll0=pcell->area_soll;
    pcell->marker.tension_index=1.0;
  }

  forcevector=(double**)malloc(3*sizeof(double*));
  for(k=0;k<3;k++){
    forcevector[k]=(double*)malloc(nb_vertex_tot*sizeof(double));
  }

  border_check();

  updatedata();
  printf("Printing data to info files...");fflush(stdout);
  infoprint(filename,0);
  printf("done!\n");fflush(stdout);

  cleanup();
  for(k=0;k<3;k++){free(forcevector[k]);}
  free(forcevector);
}

int
main()
{
  converta();  
  return 0;
}
///////////////////////////////////////////////////////////////
//
//   Utilities
//
///////////////////////////////////////////////////////////////
Bond **createbondvector(int dimension){
  Bond **a;
  a = (Bond**)malloc(dimension * sizeof(Bond*));
  if (!a) locerror("createbondvector","Not enough memory");
  return a;
}
Node 
**createnodevector(int dimension)
{
  Node**a;
  a = (Node**)malloc(dimension * sizeof(Node*));
  if (!a) locerror("createnodevector", "Not enough memory");
  return a;
}
Dual
**createdualvector(int dimension)
{
  Dual **a;
  a = (Dual **)malloc(dimension * sizeof(Dual *));
  if (!a) locerror("createdualvector", "Not enough memory");
  return a;
}
Bond *allocbond(){
  Bond *a;
  a=(Bond*)malloc(sizeof(Bond));
  if (!a) locerror("allocbond","Not enough memory");
  return a;
}
Node
*allocnode()
{
  Node*a;
  a = (Node*)malloc(sizeof(Node));
  if (!a) locerror("allocnode", "Not enough memory");
  return a;
}

Node
**allocvertexlist(int dimension)
{
  Node **a;
  a = (Node **)malloc(dimension * sizeof(Node *));
  if (!a) locerror("allocvertexlist", "Not enough memory");
  return a;
}
Dual
*allocdual()
{
  Dual *a;
  a = (Dual *)malloc(sizeof(Dual));
  if (!a) locerror("allocdual", "Not enough memory");
  return a;
}
Dual **alloccelllist(int dimension){
  Dual **a;
  a = (Dual **)malloc(dimension*sizeof(Dual*));
  if (!a) locerror("alloccelllist","Not enough memory");
  return a;
}
Bond **allocbondlist(int dimension){
  Bond **a;
  a=(Bond **)malloc(dimension*sizeof(Bond*));
  if (!a) locerror("allocbondlist","Not enough memory");
  return a;
}
int
indexvertexinweb(Node *pvertex) 
{
  int n;
  for (n=0; n < nb_vertex_tot ; n++){
    if (web[n]==pvertex) return n;
  }
  locerror("indexvertexinweb", "The vertex does not belong to the web!!");
  return 0;
}
int
indexcellinwebdual(Dual *pcell)
{
    int n;
    
    for (n=0; n < nb_cell_tot ; n++){
      if (web_dual[n]==pcell) return n;
    }
    locerror("indexcellinwebdual", "The cell does not belong to the web_dual!!");
    return 0;
}

void
locerror(char *function, char *message)
{
  fprintf(stderr, "Fatal error : function %s, %s\n", function, message);
  exit(1);
}
/* Computes the centroid of one cell (considering vertices)
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
FILE
*openfile(char *filename)
{
    FILE *pfile;

    errno=0;

    if (!(pfile=fopen(filename, "w"))){
        fprintf(stderr,"function openfile : Cannot open file\n");
        perror(NULL);
        return NULL;
    }

    return pfile;
}
double dist(Node *pvertex1, Node *pvertex2){
  register double distx, disty, distz;

  distx=pvertex1->x - pvertex2->x;
  disty=pvertex1->y - pvertex2->y;
  distz=pvertex1->z - pvertex2->z;

  return sqrt(distx*distx + disty*disty + distz*distz+0.00001);
}
double perimeter(Dual *pcell){
  int i,nbvertices;
  double perim=0;

  nbvertices=pcell->nb_vertices;

  for(i=0;i<nbvertices-1;i++){
    perim+=dist(pcell->vertexlist[i],pcell->vertexlist[i+1]);
  }

  perim+=dist(pcell->vertexlist[0],pcell->vertexlist[nbvertices-1]);
  return perim;
}
void freecell(Dual *pcell){
  free(pcell->vertexlist);
  free(pcell->celllist);
  free(pcell->bondlist);
  free(pcell);
}
void kill_lattice(){
  int i;

  for(i=0;i<nb_cell_tot;i++){
    freecell(web_dual[i]);
  }

  for(i=0;i<nb_vertex_tot;i++){
    free(web[i]);
  }

  free(web);
  free(web_dual);
}
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

int indexvertex(Node *pvertex, Dual *pcell)
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

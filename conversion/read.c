#include <stdio.h>
#include <errno.h>
#include "save.h"
#include "read.h"
#include "convert.h"

int
readweb(char *filename)
{
  int n, m;
  int ivert,icell,ivert0,ivert1,ivert2,icell0,icell1,icell2;
  FILE *pfile;
  Node *pvert;
  Dual *pcell;

  
  if ( !(pfile=fopen(filename, "ra")) ){
    fprintf(stderr,"function save : Cannot open file\n lattice not saved : ");
    perror(NULL);
    return 1;
  }
  
  /* Number of cells and vertices */
  fscanf(pfile,"%i\n",&nb_vertex_tot);
  fscanf(pfile,"%i\n",&nb_cell_tot);
  
  // Allocate memory
  web      = createnodevector(nb_vertex_tot);
  web_dual = createdualvector(nb_cell_tot);
  for (n=0;n<nb_vertex_tot;n++){ 
    web[n]=allocnode(); 
  } 
  for (n=0;n<nb_cell_tot;n++){
    web_dual[n]=allocdual();
  }
  
  //Node
  for (n=0 ; n < nb_vertex_tot ; n++){
    fscanf(pfile,"%lf %lf %i %lf %lf %lf \n",&web[n]->x,&web[n]->y,&web[n]->move,&web[n]->dist_vert[0],&web[n]->dist_vert[1],&web[n]->dist_vert[2]);
    fscanf(pfile,"%i %i %i\n",&ivert0,&ivert1,&ivert2);
    web[n]->pneighb[0]=web[ivert0];
    web[n]->pneighb[1]=web[ivert1];
    web[n]->pneighb[2]=web[ivert2];
    
    fscanf(pfile,"%i %i %i\n",&icell0,&icell1,&icell2);
    web[n]->pncell[0]=web_dual[icell0];
    web[n]->pncell[1]=web_dual[icell1];
    web[n]->pncell[2]=web_dual[icell2];

    web[n]->z=0.0;
  }
  
  
  //Dual
  for (n=0 ; n < nb_cell_tot ; n++){ 
    pcell=web_dual[n];
    fscanf(pfile,"%i\n",&pcell->nb_vertices);
    
    web_dual[n]->vertexlist=allocvertexlist(pcell->nb_vertices);
    for (m=0 ; m < pcell->nb_vertices ; m++){
      fscanf(pfile,"%i\n",&ivert0);
      pcell->vertexlist[m]=web[ivert0];
    }
  }
  for(n=0;n<nb_cell_tot;n++){
    pcell=web_dual[n];
    pcell->sqrtarea=sqrt(area(pcell));
    pcell->perimeter = perimeter(pcell);
    pcell->centroid = centroid(pcell);
    pcell->marker.clone_index=0;
  }

  int numclone;
  int tempvar;
  fscanf(pfile,"%i\n",&numclone);
  for(n=0;n<numclone;n++){
      fscanf(pfile,"%i\n",&tempvar);
      pcell=web_dual[tempvar];
      pcell->marker.clone_index=1;
  }
  
  
  fclose(pfile);

  shrink_data();  
  return 0;
}

#include <stdio.h>
#include <errno.h>
#include "save.h"
#include "convert.h"

/* Record the lattice information in a binary file
 */
int
saveweb(char *filename)
{
  int n, m;
  int ivert,icell;
  FILE *pfile;
  Dual *pcell;
  Node *pvert;
  
  
  if ( !(pfile=fopen(filename, "wa")) ){
      fprintf(stderr,"function save : Cannot open file\n lattice not saved : ");
      perror(NULL);
      return 1;
    }
  
  /* Number of cells and vertices */
  
  fprintf(pfile,"%i\n",nb_vertex_tot);
  fprintf(pfile,"%i\n",nb_cell_tot);
  
  //Node
  for (n=0 ; n < nb_vertex_tot ; n++){
    fprintf(pfile,"%lf %lf %i %lf %lf %lf \n",web[n]->x,web[n]->y,web[n]->move,web[n]->dist_vert[0],web[n]->dist_vert[1],web[n]->dist_vert[2]);
    fprintf(pfile,"%i %i %i\n",indexvertexinweb(web[n]->pneighb[0]),indexvertexinweb(web[n]->pneighb[1]),indexvertexinweb(web[n]->pneighb[2]));
    fprintf(pfile,"%i %i %i\n",indexcellinwebdual(web[n]->pncell[0]),indexcellinwebdual(web[n]->pncell[1]),indexcellinwebdual(web[n]->pncell[2]));
  }
  
  //Dual
  for (n=0 ; n < nb_cell_tot ; n++){ 
    pcell=web_dual[n];
    fprintf(pfile,"%i\n",pcell->nb_vertices);
    
    for (m=0 ; m < pcell->nb_vertices ; m++){
      fprintf(pfile,"%i\n",indexvertexinweb(pcell->vertexlist[m]));
    }
  }
  
  fclose(pfile);
  
  return 0;
}

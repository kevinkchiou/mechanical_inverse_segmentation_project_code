#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "convert.h"
#include "locerror.h"
#include "save.h"
#include "correction.h"


void orientation_correction(){
  int i,j,k,l,count,flag,chk;
  Dual *pncell[2],*pcell0;
  Node *pnvert[2],*pvert,**newvertexlist;

  pcell0=web_dual[0];
  newvertexlist=allocvertexlist(pcell0->nb_vertices); //new vertex list we want for pcell0
  int smallestvert=nb_vertex_tot;
  int smallestvertindex=nb_vertex_tot;
  for(k=0;k<pcell0->nb_vertices;k++){
    //find the smallest global indexed vertex and start there
    if(indexvertexinweb(pcell0->vertexlist[k])<smallestvert){
      smallestvert = indexvertexinweb(pcell0->vertexlist[k]);
      smallestvertindex = indexvertex(pcell0->vertexlist[k],pcell0);
    }
  }
  if(smallestvertindex >= pcell0->nb_vertices){locerror("Bad counting!","correction");}
  newvertexlist[0] = pcell0->vertexlist[smallestvertindex];

  //now we've found the smallest vertex. we next iterate a process around pcell0 to
  //determine the following vertices in the newvertexlist[] array
  for(l=0;l<pcell0->nb_vertices-1;l++){
    pvert = newvertexlist[l];count=0;flag=0;chk=0;
    //find neighboring cells (that are not pcell0!)
    for(i=0;i<3;i++){
      if(pvert->pncell[i]!=pcell0){pncell[count] = pvert->pncell[i];count++;}
    }
    //finds the next vertex in each cell assuming counterclockwise orientation and progression.
    //printf("we are at vertex = %d, pncell0 = %d, pncell1 = %d\n",indexvertexinweb(pvert),indexcellinwebdual(pncell[0]),indexcellinwebdual(pncell[1]));fflush(stdout);
    pnvert[0] = pncell[0]->vertexlist[(indexvertex(pvert,pncell[0])+1) % pncell[0]->nb_vertices];
    pnvert[1] = pncell[1]->vertexlist[(indexvertex(pvert,pncell[1])+1) % pncell[1]->nb_vertices];

    //if that vertex borders pcell0, then it's the cell we want and the vertex we want
    for(i=0;i<3;i++){
      if(pnvert[0]->pncell[i] == pcell0){flag=1;chk=2;}
      if(pnvert[1]->pncell[i] == pcell0){flag=2;}
    }
    if(flag==0){
	printf("We have pnvert1 = %d, pnvert2 = %d\n",indexvertexinweb(pnvert[0]),indexvertexinweb(pnvert[1]));fflush(stdout);
	locerror("We have a problem","correction");}
    if(flag==chk){locerror("We have another problem","correction");}
    if(flag==1){newvertexlist[l+1] = pnvert[0];}
    if(flag==2){newvertexlist[l+1] = pnvert[1];}
    fflush(stdout);
  }
  pcell0->vertexlist = newvertexlist; //reset the list to the new one
}

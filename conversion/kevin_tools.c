#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>

#include "const.h"
#include "convert.h"
#include "save.h"
#include "locerror.h"
#include "kevin_tools.h"
#include "input.h"

int flagt1;
int flagt1vertex;
int flagt1neighb;

double **forcevector;

/*calculates the scalar area of a cell specified */
double area(Dual *pcell){
	int ivert,nbvert;
	vector3 av,avtemp,center;
	Node *pv,*pvprev;
	
	nbvert = pcell->nb_vertices;
	pvprev = pcell->vertexlist[nbvert-1];
	pv = pcell->vertexlist[0];
	center = centroid(pcell);
	av = area_calc(pvprev,pv,center);
	
	for(ivert=1;ivert<nbvert;ivert++){
		pvprev=pv;
		pv=pcell->vertexlist[ivert];
		avtemp=area_calc(pvprev,pv,center);
		av.x+=avtemp.x;
		av.y+=avtemp.y;
		av.z+=avtemp.z;
	}
	return 0.5*sqrt(pow(av.x,2)+pow(av.y,2)+pow(av.z,2));
}
	
//calculates the partial area of a cell specified starting at ivert0cell0
//and coords (x,y,z) up to a maximum number of vertices max
double part_area(Dual *pcell, int max, int ivert0cell0, double x, double y, double z){
	int ivert,nbvert;
	vector3 av,avtemp,center;
	Node *pv,*pvprev;
   double x1,x2,y1,y2,z1,z2;
	
   nbvert = pcell->nb_vertices;
	pvprev = pcell->vertexlist[(ivert0cell0)%nbvert];
	pv = pcell->vertexlist[(ivert0cell0+1)%nbvert];
	center = centroid(pcell);

   x1=x-center.x;x2=(pv->x)-center.x;
   y1=y-center.y;y2=(pv->y)-center.y;
   z1=z-center.z;z2=(pv->z)-center.z;
   av.x = y1*z2-z1*y2;
   av.y = z1*x2-x1*z2;
   av.z = x1*y2-y1*x2;
	
	for(ivert=1;ivert<max;ivert++){
		pvprev=pv;
		pv=pcell->vertexlist[(ivert0cell0+ivert+1)%nbvert];
		avtemp=area_calc(pvprev,pv,center);
		av.x+=avtemp.x;
		av.y+=avtemp.y;
		av.z+=avtemp.z;
	}
   pvprev=pv;
   pv=pcell->vertexlist[(ivert0cell0+max)%nbvert];
   x=0.5*(pvprev->x + pv->x);
   y=0.5*(pvprev->y + pv->y);
   z=0.5*(pvprev->z + pv->z);
   //note here we switched the order to maintain orientation
   x2=x-center.x;x1=(pv->x)-center.x;
   y2=y-center.y;y1=(pv->y)-center.y;
   z2=z-center.z;z1=(pv->z)-center.z;
   av.x += y1*z2-z1*y2;
   av.y += z1*x2-x1*z2;
   av.z += x1*y2-y1*x2;
	return 0.5*sqrt(pow(av.x,2)+pow(av.y,2)+pow(av.z,2));
}

/*calculates area vector of triangle spanned by points represented
  by the triangle calculated by the cross product pv1 x pv2 about
  the point centroid */
vector3 area_calc(Node *pv1, Node *pv2, vector3 centroid){
	double x1,x2,y1,y2,z1,z2;
	vector3 av; //area vector
	
	x1=(pv1->x)-(centroid.x);
	x2=(pv2->x)-(centroid.x);
	y1=(pv1->y)-(centroid.y);
	y2=(pv2->y)-(centroid.y);
	z1=(pv1->z)-(centroid.z);
	z2=(pv2->z)-(centroid.z); 
	
	av.x = y1*z2-z1*y2;
	av.y = z1*x2-x1*z2;
	av.z = x1*y2-y1*x2;
	
	return av;
}

void latticeprint(char *filename, unsigned int type){
	int i,ivert;
	FILE *pfile;
	Dual *pcell;
	Node *pvert;
	vector3 center;
	double xmin,xmax,ymin,ymax,zmin,zmax,x,y,z;
	
	pfile = openfile(filename);
	fprintf(pfile,"\n");
	
	for(i=1;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		fprintf(pfile,"%lf\n",0); //normally the pressure that is printed here
		pvert=pcell->vertexlist[0];
		for(ivert=0;ivert<pcell->nb_vertices;ivert++){
			pvert=pcell->vertexlist[ivert];
			if(ivert==0){
				xmin = pvert->x;
				xmax = pvert->x;
				ymin = pvert->y;
				ymax = pvert->y;
				zmin = pvert->z;
				zmax = pvert->z;
			}
			else{
				if (xmin > pvert->x) xmin = pvert->x;
				if (xmax < pvert->x) xmax = pvert->x;
				if (ymin > pvert->y) ymin = pvert->y;
				if (ymax < pvert->y) ymax = pvert->y;
				if (zmin > pvert->z) zmin = pvert->z;
				if (zmax < pvert->z) zmax = pvert->z;
			}
			x=pvert->x;y=pvert->y;z=pvert->z;
			fprintf(pfile,"%lf %lf %lf\n",x,y,z);
		}
		fprintf(pfile,"\n");
		fprintf(pfile,"\n");
	}
	fprintf(pfile,"END\n");
	fclose(pfile);
}

int indexbondinwebbond(Bond *pbond){
   int i;

   for(i=0;i<nb_bond_tot;i++){
      if(web_bond[i]==pbond){return i;}
   }
   locerror("indexbondinwebbond","The bond does not exist");
   return 0;
}

void infoprint(char *filename, unsigned int type){
   int i,j,k,num;
   char filename1[256],filename2[256],filename3[256];
   FILE *pfile1,*pfile2,*pfile3;
   Dual *pcell,*pcell0,*pcell2;
   Node *pvert;
   Bond *pbond;
   vector3 center;
   double xmin,xmax,ymin,ymax,zmin,zmax,x,y,z;
   int debug=0;//change this to 1 to output some debugger stuff

   snprintf(filename1,255,"info/%s.cellinfo",filename);
   snprintf(filename2,255,"info/%s.bondinfo",filename);
   snprintf(filename3,255,"info/%s.vertinfo",filename);

   pfile1=openfile(filename1);pfile2=openfile(filename2);pfile3=openfile(filename3);
   fprintf(pfile1,"%d\n",nb_cell_tot);
   fprintf(pfile2,"%d\n",nb_bond_tot);
   fprintf(pfile3,"%d\n",nb_vertex_tot);
   for(i=0;i<nb_cell_tot;i++){
      pcell=web_dual[i];
      num = pcell->nb_vertices;
      fprintf(pfile1,"%d\t%d\t",i,num);
      if(debug==1){printf("cell i = %d, neighbor ",i);fflush(stdout);}
      for(j=0;j<num;j++){
	 if(debug==1){printf("%d...",j);fflush(stdout);}
         if(j!=num-1){fprintf(pfile1,"%d ",indexcellinwebdual(pcell->celllist[j]));}
         if(j==num-1){fprintf(pfile1,"%d\t",indexcellinwebdual(pcell->celllist[j]));}
      }
      for(j=0;j<num;j++){
         if(j!=num-1){fprintf(pfile1,"%d ",indexbondinwebbond(pcell->bondlist[j]));}
         if(j==num-1){fprintf(pfile1,"%d\t",indexbondinwebbond(pcell->bondlist[j]));}
      }
      for(j=0;j<num;j++){
         if(j!=num-1){fprintf(pfile1,"%d ",indexvertexinweb(pcell->vertexlist[j]));}
         if(j==num-1){fprintf(pfile1,"%d\t",indexvertexinweb(pcell->vertexlist[j]));}
      }
      if(debug==1){printf("\n");fflush(stdout);}
      center=centroid(pcell);
      fprintf(pfile1,"%lf %lf %lf %d ",0.0,area(pcell),pcell->area_soll,pcell->marker.clone_index); //normally the pressure is printed here
      fprintf(pfile1,"%lf %lf %lf\n",center.x,center.y,center.z);
   }
   fclose(pfile1);
   for(i=0;i<nb_bond_tot;i++){
      pbond = web_bond[i];
      fprintf(pfile2,"%d\t",i);
      for(j=0;j<4;j++){
         if(j!=3){fprintf(pfile2,"%d ",indexcellinwebdual(pbond->pncell[j]));}
         if(j==3){fprintf(pfile2,"%d\t",indexcellinwebdual(pbond->pncell[j]));}
      }
      for(j=0;j<4;j++){
         if(j!=3){fprintf(pfile2,"%d ",indexbondinwebbond(pbond->pnbond[j]));}
         if(j==3){fprintf(pfile2,"%d\t",indexbondinwebbond(pbond->pnbond[j]));}
      }
      for(j=0;j<2;j++){
         if(j!=1){fprintf(pfile2,"%d ",indexvertexinweb(pbond->pnvert[j]));}
         if(j==1){fprintf(pfile2,"%d\t",indexvertexinweb(pbond->pnvert[j]));}
      }
      fprintf(pfile2,"%lf\t",dist(pbond->pnvert[0],pbond->pnvert[1]));
      pcell0=pbond->pncell[0];pcell2=pbond->pncell[2];
      fprintf(pfile2,"%lf\n",pcell0->marker.tension_index+pcell2->marker.tension_index);
   }
   fclose(pfile2);
   for(i=0;i<nb_vertex_tot;i++){
      pvert = web[i];
      fprintf(pfile3,"%d\t",i);
      for(j=0;j<3;j++){
         if(j!=2){fprintf(pfile3,"%d ",indexcellinwebdual(pvert->pncell[j]));}
         if(j==2){fprintf(pfile3,"%d\t",indexcellinwebdual(pvert->pncell[j]));}
      }
      for(j=0;j<3;j++){
         if(j!=2){fprintf(pfile3,"%d ",indexbondinwebbond(pvert->pnbond[j]));}
         if(j==2){fprintf(pfile3,"%d\t",indexbondinwebbond(pvert->pnbond[j]));}
      }
      for(j=0;j<3;j++){
         if(j!=2){fprintf(pfile3,"%d ",indexvertexinweb(pvert->pneighb[j]));}
         if(j==2){fprintf(pfile3,"%d\t",indexvertexinweb(pvert->pneighb[j]));}
      }
      fprintf(pfile3,"%lf %lf %lf\t",pvert->x,pvert->y,pvert->z);
      fprintf(pfile3,"%lf %lf %lf\n",forcevector[0][i],forcevector[1][i],forcevector[2][i]);
   }
   fclose(pfile3);
}

//counts the total number of bonds and sets all the
//pointers to bonds for each vertex to NULL to make
//the next part (where i call initbonds()) simpler
void countbonds(){
   int i,j,nb_bonds=0;
   Node *pvert;

   for(i=0;i<nb_vertex_tot;i++){
      for(j=0;j<3;j++){
         pvert = web[i]->pneighb[j];
         if(indexvertexinweb(pvert) > i){nb_bonds++;}
      }
   }
   nb_bond_tot = nb_bonds;
}

//initializes the memory for each bond as well as
//the pointers for the endpoint vertices as well
//as bond pointers for the vertices
void initbonds(){
   int i,j,k=0,l;
   Node *pvert,*pnvert;
   Bond *pbond;

   for(i=0;i<nb_vertex_tot;i++){
      pvert = web[i];
      for(j=0;j<3;j++){
         pnvert = pvert->pneighb[j];
         //sees if the neighboring vertex has greater overall index
         //if so, we create a new bond allocation and initialize
         //all the parameters as well as if the neighboring vertex
         //has index j, then the corresponding bond between the vertices
         //also has index j as a neighbor
         if(indexvertexinweb(pnvert) > i){
            web_bond[k] = allocbond();
            pbond = web_bond[k];
            pbond->pnvert[0] = pvert;
            pbond->pnvert[1] = pnvert;
            pvert->pnbond[j] = pbond;
            for(l=0;l<3;l++){
               //keeps a consistent index between neighboring vertices
               //and neighboring bonds, but in reverse as a neighbor
               if(pnvert->pneighb[l]==pvert){pnvert->pnbond[l] = pbond;}
            }
            k++;//increment the web_bond number
         }
      }
   }
   if(k!=nb_bond_tot){locerror("initbonds","Error in bond counting");}
}

//initializes the celllist and bondlist
//pointer arrays in the node structure
//also contains the algorithm for finding
//the bounding bonds of a cell
void dualneighb(){
   int i,j,k,nb_vert;
   Dual *pcell,*pcell1;
   Node *pvert,*pnvert;
   Bond *pbond;

   for(i=0;i<nb_cell_tot;i++){
      pcell = web_dual[i];
      nb_vert = pcell->nb_vertices;
      pcell->celllist = alloccelllist(nb_vert);
      pcell->bondlist = allocbondlist(nb_vert);
      for(j=0;j<nb_vert;j++){
         //looks at the vertices of a cell
         pvert = pcell->vertexlist[j];
         pnvert = pcell->vertexlist[(j+1)%nb_vert];
         for(k=0;k<3;k++){
            //sees which neighbor of the vertex qualifies as the corresponding
            //neighbor in the cell.  then calls the bond between vertices j and j+1
            //the bordering bond j of the cell. can be done since above the corresp
            //bond and vertex neighbors have the same index.
            if(pnvert == pvert->pneighb[k]){
               pcell->bondlist[j] = pvert->pnbond[k];
            }
         }
      }
   }
}

//algorithm for determining the neighboring bonds and cells to a bond
void bondneighb(){
   int i,j,k,l=0,m=2,n;
   Node *pvert1, *pvert2;
   Dual *pcell1, *pcell2, *pcell3, *pcell4;
   Bond *pbond;

   for(i=0;i<nb_bond_tot;i++){
      pbond = web_bond[i];
      pvert1 = pbond->pnvert[0];pvert2 = pbond->pnvert[1];
      l=0;m=2;//initialize the array counting
      for(k=0;k<3;k++){ //this does the neighboring bonds
         if(pvert1->pnbond[k]!=pbond){pbond->pnbond[l] = pvert1->pnbond[k];l++;}
         if(pvert2->pnbond[k]!=pbond){pbond->pnbond[m] = pvert2->pnbond[k];m++;}
         if(l>2){
	 printf("l = %d, pvert1 = %d, pvert2 = %d, pbond = %d\n",l,indexvertexinweb(pvert1),indexvertexinweb(pvert2),indexbondinwebbond(pbond));fflush(stdout);
	 locerror("bondneighb(l)","Error in neighboring bond assignment");}
         if(m>4){
	 printf("m = %d, pvert1 = %d, pvert2 = %d, pbond = %d\n",m,indexvertexinweb(pvert1),indexvertexinweb(pvert2),indexbondinwebbond(pbond));fflush(stdout);
	 locerror("bondneighb(m)","Error in neighboring bond assignment");}
      }
      l=0; //initialize the array index
      for(k=0;k<3;k++){ //this does the neighboring cells
         n=0;
         while(n<3){ //if the neighbors match, we assign them into pncell[0] or pncell[2]
            if(pvert1->pncell[k]==pvert2->pncell[n]){pbond->pncell[l]=pvert1->pncell[k];l+=2;break;}
            n++;
         }
      }
      for(k=0;k<3;k++){ //if the neighbors aren't either of the two above, they are the ends
         if((pvert1->pncell[k]!=pbond->pncell[0]) && (pvert1->pncell[k]!=pbond->pncell[2])){
            pbond->pncell[1] = pvert1->pncell[k];
         }
         if((pvert2->pncell[k]!=pbond->pncell[0]) && (pvert2->pncell[k]!=pbond->pncell[2])){
            pbond->pncell[3] = pvert2->pncell[k];
         }
      }
   }
}

//this is here since this is much easier to do after the bond
//neighbors have been established due to our configuration
void neighbcell(){
   int i,j,k,nb_bonds;
   Bond *pbond;
   Dual *pcell,*pncell;

   for(i=0;i<nb_cell_tot;i++){
      pcell = web_dual[i];
      nb_bonds = pcell->nb_vertices;
      for(j=0;j<nb_bonds;j++){
         pbond = pcell->bondlist[j];
	 if (pbond==NULL || pbond->pncell==NULL){
	    printf("Null pointer!\n");fflush(stdout);
	    printf("It occured for pcell = %d and at ",indexcellinwebdual(pcell));
	    printf("index j = %d.\n",j);
	    fflush(stdout);}
         if(pbond->pncell[0]!=pcell){pcell->celllist[j] = pbond->pncell[0];}
         if(pbond->pncell[2]!=pcell){pcell->celllist[j] = pbond->pncell[2];}
      }
   }
}

void updatedata(){
   countbonds();
   printf("There are %d bonds, %d cells, and %d vertices\n",nb_bond_tot,nb_cell_tot,nb_vertex_tot);fflush(stdout);
   web_bond = createbondvector(nb_bond_tot);
   initbonds(); //initialize bonds
   dualneighb(); //find neighbor info for cells
   bondneighb(); //find neighboring info for bonds
   neighbcell(); //find neighboring cells to cells
}

//cleans up the stuff that's created by kevin_tools.c
void cleanup(){
   int i;
   Dual *pcell;

   for(i=0;i<nb_cell_tot;i++){
      pcell = web_dual[i];
      free(pcell->bondlist);
      free(pcell->celllist);
   }
   for(i=0;i<nb_bond_tot;i++){
      free(web_bond[i]);
   }
   free(web_bond);
}

//gives each vertex an offset
void coord_mod(){
    int i;
    Node *pvert;

    for(i=0;i<nb_vertex_tot;i++){
       pvert = web[i];
       pvert->z = 12.0;
    }
}

void border_check(){
   int i,j,k;
   Dual *pcell,*pcellnb;
   Node *pvert;

   for(i=0;i<nb_vertex_tot;i++){
      pvert = web[i];
      pvert->border = 0;
   }
   for(i=0;i<nb_cell_tot;i++){
      pcell = web_dual[i];
      pcell->marker.border=0;
   }

   /*
   for(i=0;i<nb_cell_tot;i++){
      pcell = web_dual[i];
      if((pcell->marker.tension_index)>1.001){
         pcell->area_soll=0.8;
         pcell->area_soll0=pcell->area_soll;
         for(j=0;j<(pcell->nb_vertices);j++){
            pvert=pcell->vertexlist[j];
            for(k=0;k<3;k++){
               pcellnb = pvert->pncell[k];
               if((pcellnb->marker.tension_index)<1.001){
                  pvert->border=1;
                  pcell->marker.border=1;
                  pcell->area_soll=0.8;
                  pcell->area_soll0=pcell->area_soll;
               }
            }
         }
      }
   }
   */
}


/* The comments directly come from the libpng web site...
 *
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "convert.h"
#include "pngwrite.h"
#include "kevin_tools.h"

FILE *fp;


/* Svg file generation : takes a snapshop of the lattice
 * It is very usefull to test the program and/or to know what happen when there is a bug
 */
void
svglattice(char *filename, unsigned int type)
{
   int i;
    FILE *pfile;
    Dual *pcell;

    /* File initialisation */
    pfile = openfile(filename);

    fprintf(pfile, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n\
<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n\
   \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\
<svg\n   xmlns=\"http://www.w3.org/2000/svg\"\n\
   version=\"1.1\"\n\
   x=\"0mm\"\n\
   y=\"0mm\"\n\
   width=\"250mm\"\n\
   height=\"250mm\">\n\
   <desc>Growing lattice\n\
   </desc>\n\
   <defs>\n\
   </defs>\n\
   <g transform=\"matrix(1,0,0,-1,0,0)\">\n\
   <g transform=\"translate(100,-500)\">\n\
   <g transform=\"scale(15)\"\n\
      style=\"fill:white;stroke:black;stroke-width:0.05;stroke-linecap:round;stroke-linejoin:round\">\n");

    for (i=1;i<nb_cell_tot;i++){
       pcell=web_dual[i];
       svgprcell(pcell,pfile,type);
    }
    printgrad(pfile);
    fprintf(pfile, "      </g>\n </g>\n  </g>\n</svg>");
    fclose(pfile);
}

/* Eps file generation : takes a snapshop of the lattice
 * It is very usefull to test the program and/or to know what happen when there is a bug
 */
void
epsfile(char *filename, unsigned int type)
{
    FILE *pfile;
    int icell;	/* Loop variables */
    time_t epstime;

    /* File initialisation */
    pfile = openfile(filename);

    epstime=time(NULL);
    /* Header of the postscript file */
    /* Comments */
    fprintf(pfile, "%%!PS-Adobe-2.0 EPSF-2.0\n");
    fprintf(pfile, "%%%%Title: %s\n", filename);
    fprintf(pfile, "%%%%Creator: Morphogenesis\n");
    fprintf(pfile, "%%%%CreationDate: %s", asctime(gmtime(&epstime)));
    fprintf(pfile, "%%%%BoundingBox: -200 -200 300 300\n");
    fprintf(pfile, "%%%%EndComments\n");

    fprintf(pfile, "/unit {10 mul} def\n");
/*    fprintf(pfile, "100 100 translate\n");*/
    fprintf(pfile, "/clonedot {\ngsave\nnewpath\n0 setgray\nexch unit exch unit translate\n0 0 .08 unit 0 360 arc closepath fill\ngrestore\n}def\n");
    fprintf(pfile, "/Times-Roman findfont\n");
    fprintf(pfile, "/ccol{0 0 0 setrgbcolor}def\n");
    fprintf(pfile, "/lcol{0 0 0 setrgbcolor}def\n");
    fprintf(pfile, "/cslw{0.5 setlinewidth}def\n");
    fprintf(pfile, "/lslw{0.5 setlinewidth}def\n");
    fprintf(pfile, "/m{exch unit exch unit moveto}def\n");
    fprintf(pfile, "/l{exch unit exch unit lineto}def\n");
    fprintf(pfile, "/n{newpath}def\n");
    fprintf(pfile, "/c{closepath}def\n");
    fprintf(pfile, "/s{stroke}def\n");
    fprintf(pfile, "/f{fill}def\n");
    fprintf(pfile, "/pf{2.5}def\n");
    fprintf(pfile, "/scor{neg 2 mul 3 add exch neg 2 mul 3 add setrgbcolor}def\n");
    fprintf(pfile, "/scob{neg 1 add exch neg 1 add  1. setrgbcolor}def\n");

    fprintf(pfile, "6 scalefont\n");
    fprintf(pfile, "setfont\n");
    fprintf(pfile, "0.2 setlinewidth\n");
    

    /* Draw all cells */
    for (icell=1 ; icell<nb_cell_tot; icell++){
	  //print_cell(web_dual[icell], pfile, type);
    }

    /* File closure */
    fprintf(pfile, "showpage\n");
    fprintf(pfile, "%%%%Trailer\n");
    fprintf(pfile, "%%EOF\n");
    fclose(pfile);
}

void
svgprcell(Dual *pcell, FILE *pfile, unsigned int type)
{
  int ivert;
  Node *pvert;

  fprintf(pfile, "         <path\n");
  fprintf(pfile, "            id=\"cell%i\"\n", indexcellinwebdual(pcell));
  switch (type){
     case 0:
	printdefects(pcell, pfile);
	break;
     case 1:
	printpressure(pcell, pfile);
	break;
     case 2:
	printthick(pcell, pfile);
	break;
     case 3:
	printgradcell(pcell, pfile);
        break;
     case 4:
	printampli(pcell, pfile);
	break;
     case 5:
        break;
  }
  pvert=pcell->vertexlist[0];
  fprintf(pfile, "            d=\"M %.2f,%.2f", pvert -> x, pvert -> y);
  for (ivert=1; ivert< pcell->nb_vertices; ivert++){
    pvert=pcell->vertexlist[ivert];
    fprintf(pfile, " L %f,%f", pvert->x, pvert->y);
  }
  fprintf(pfile, " z\" />\n");
  if (pcell->marker.clone_index){
     pcell->centroid=centroid(pcell);
     fprintf(pfile, "            <circle cx=\"%.2f\" cy=\"%.2f\" r=\"0.1\" fill=\"black\" />",\
	   pcell->centroid.x, pcell->centroid.y);
  }
}

void
printpressure(Dual *pcell, FILE *pfile)
{
   double pres;
   
   //pres=pressure(pcell);
   pres=1.;
   
   if (pres<0) fprintf(pfile, "            fill=\"rgb(%1$i,%1$i,255)\"\n", (int)(255*(1+0.5*pres)));
   else fprintf(pfile, "            fill=\"rgb(255,%1$i,%1$i)\"\n", (int)(255*(1-0.5*pres)));
}

void
printgradcell(Dual *pcell, FILE *pfile)
{
   double grad;
   
   //grad=diffgrad(pcell);
   grad=1.;
   
   fprintf(pfile, "            fill=\"rgb(%1$i,%1$i,255)\"\n", (int)(255.0*(1.0-0.08*grad)));
}

void
printthick(Dual *pcell, FILE *pfile)
{
   if (pcell->thick>1.0) fprintf(pfile, "            fill=\"rgb(255,%1$i,%1$i)\"\n", (int)(255*(1-2*(pcell->thick-1))));
   else fprintf(pfile, "            fill=\"rgb(%1$i,%1$i,255)\"\n", (int)(255*(1+6*(pcell->thick-1))));
}

void
printampli(Dual *pcell, FILE *pfile)
{
   fprintf(pfile, "            fill=\"rgb(255,%1$i,%1$i)\"\n", (int)(255*(1-pcell->ampli)));
}

void
printgrad(FILE *pfile)
{
   int i, j;
   Node *pvert, *pvertn;
   Dual *pcell0, *pcell1, *pcell2;
   pcell0=web_dual[0];
   fprintf(pfile, "         <g stroke=\"moccasin\" stroke-width=\"0.15\">\n");
   for (i=0;i<nb_vertex_tot;i++){//vert!!!
      pvert=web[i];
      for (j=0;j<3;j++){
	 pcell1=pvert->pncell[j];
	 pcell2=pvert->pncell[(j+2)%3];
	 pvertn=pvert->pneighb[j];
	 if (pcell1!=pcell0 && pcell2!=pcell0 && pcell1->marker.difftype!=pcell2->marker.difftype){
	    fprintf(pfile, "            <path d=\"M %f,%f L %f,%f\"/>\n", pvert->x, pvert->y, pvertn->x, pvertn->y);
	 }
      }
   }
   fprintf(pfile, "         </g>\n");
}

void
printdefects(Dual *pcell, FILE *pfile)
{
   switch (pcell->nb_vertices){
      case 5:
	 fprintf(pfile, "            fill=\"red\"\n");
	 break;
      case 7:
	 fprintf(pfile, "            fill=\"blue\"\n");
	 break;
   }
}


/* Eps file generation : takes a snapshop of the lattice
 * It is very usefull to test the program and/or to know what happen when there is a bug
 */
void
epsfile_pressure_dpp(char *filename)
{

    FILE *pfile;
    int icell;	/* Loop variables */
    double avg_press;
    time_t epstime;
    Dual *pcell;


    /* File initialisation */
    pfile = openfile(filename);

    epstime=time(NULL);
    /* Header of the postscript file */
    /* Comments */
    fprintf(pfile, "%%!PS-Adobe-2.0 EPSF-2.0\n");
    fprintf(pfile, "%%%%Title: %s\n", filename);
    fprintf(pfile, "%%%%Creator: Morphogenesis\n");
    fprintf(pfile, "%%%%CreationDate: %s", asctime(gmtime(&epstime)));
    fprintf(pfile, "%%%%BoundingBox: -200 -200 300 300\n");
    fprintf(pfile, "%%%%EndComments\n");

    fprintf(pfile, "/unit {10 mul} def\n");
/*    fprintf(pfile, "100 100 translate\n");*/
    fprintf(pfile, "/clonedot {\ngsave\nnewpath\n0 setgray\nexch unit exch unit translate\n0 0 .08 unit 0 360 arc closepath fill\ngrestore\n}def\n");
    fprintf(pfile, "/Times-Roman findfont\n");
    fprintf(pfile, "/ccol{0 0 0 setrgbcolor}def\n");
    fprintf(pfile, "/lcol{0 exch 0 setrgbcolor}def\n");
    fprintf(pfile, "/cslw{0.5 setlinewidth}def\n");
    fprintf(pfile, "/lslw{0.5 setlinewidth}def\n");
    fprintf(pfile, "/m{exch unit exch unit moveto}def\n");
    fprintf(pfile, "/l{exch unit exch unit lineto}def\n");
    fprintf(pfile, "/n{newpath}def\n");
    fprintf(pfile, "/c{closepath}def\n");
    fprintf(pfile, "/s{stroke}def\n");
    fprintf(pfile, "/f{fill}def\n");
    fprintf(pfile, "/pf{2.5}def\n");
    fprintf(pfile, "/scor{pf mul neg 1 add exch pf mul neg 1 add  1. 3 1 roll setrgbcolor}def\n");
    fprintf(pfile, "/scob{pf mul 1 add exch pf mul 1 add  1. setrgbcolor}def\n");

    fprintf(pfile, "6 scalefont\n");
    fprintf(pfile, "setfont\n");
    fprintf(pfile, "0.2 setlinewidth\n");
    
    avg_press = 0.;

    /* Draw all background cells */
    for (icell=1 ; icell<nb_cell_tot; icell++){
	pcell=web_dual[icell];
	//fprintf(pfile, "%% vertex %i %f\n",icell,pcell->label);
	//print_cell_dpp(pcell,avg_press, pfile);
    }


    /* File closure */
    fprintf(pfile, "showpage\n");
    fprintf(pfile, "%%%%Trailer\n");
    fprintf(pfile, "%%EOF\n");
    fclose(pfile);
}



void
print_cell_dpp(Dual *pcell, double avg_press, FILE *pfile)
{
  int ivertex;
  int nbvertices;
  double press;
  double dist, lcol;
  vector2 centroid_cell;
  Node *pvertex;

  nbvertices=pcell->nb_vertices;
  //press=pressure(pcell);
  press=1.;
  //press=pcell->pressure;
  //press -= avg_press;

  centroid_cell=centroid(pcell);
  dist = centroid_cell.x*centroid_cell.x + centroid_cell.y*centroid_cell.y;


  pvertex=pcell->vertexlist[0];
  fprintf(pfile, "n\n");
  fprintf(pfile, "%f %f m\n", pvertex -> x, pvertex -> y);

  for (ivertex=1; ivertex < nbvertices; ivertex++){
    pvertex=pcell->vertexlist[ivertex];
    fprintf(pfile, "%f %f l\n", pvertex -> x, pvertex -> y);
  }

  fprintf(pfile, "c\n");
  fprintf(pfile, "gsave\n");

  if (press>0){
    //    fprintf(pfile, "1.0 %f %f setrgbcolor\n", 1.0-(prefac*press), 1.0-(prefac*press));
    fprintf(pfile, "%f %f scor\n", press, press);
  }
  else {
    //    fprintf(pfile, "%f %f 1.0 setrgbcolor\n", 1.0+(prefac*press), 1.0+(prefac*press));
    fprintf(pfile, "%f %f scob\n", press, press);
  }
  fprintf(pfile, "f\n");

  fprintf(pfile, "grestore\n");
  if (dist > DISC_RADIUS){
    lcol = 0;
  }
  else {
    lcol = 1;
  }
  lcol = 1.-dist/200.;
  /*if (ttime-pcell->birth < 1000) {
    fprintf(pfile, "%f %f clonedot\n", pcell->centroid_cell.x, pcell->centroid_cell.y);
  }*/
  fprintf(pfile, "lslw\n");
  fprintf(pfile, "%f lcol\n",lcol);
  fprintf(pfile, "stroke\n");
}

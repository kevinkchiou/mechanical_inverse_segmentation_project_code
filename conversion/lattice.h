#ifndef SRC_LATTICE_LATTICE_H
#define SRC_LATTICE_LATTICE_H 1

#include "const.h" /* Many preprocessor macros are used in this file */

/*********************************************************************/
/* We expect this new code to have 3d coordinates, but it is written */
/* with the expectation that the lattice structure lies primarily in */
/* 2d. So coords are 3d, but cells are coded considering 2d structure*/
/*********************************************************************/

/* These two lines define the # of vertices and cells on the lattice :
 * it changes at each division
 */
extern int nb_vertex_tot;
extern int nb_cell_tot;
extern int nb_bond_tot;
extern int nb_clone_tot;
extern int nb_cell_killed;
extern int flagt1;
extern int flagt1vertex;
extern int flagt1neighb;

extern double add_factor1;
extern double add_factor2;
extern double add_factor;

/* Definition of a 3x3 matrix
 */
typedef struct {double c[3][3];} matrix33;

/* Definition of a vector with 3 components (x, y, and z)
 */
typedef struct
{
    double x;
    double y;
    double z;
} vector3;

//definition of a bond!
typedef struct bond
{
   struct node *pnvert[2]; //pointers to vertex endpoints
   struct dual *pncell[4]; //pointers to neighboring cells
   struct bond *pnbond[4]; //pointers to neighboring bonds
} Bond;

/* Definition of a node !! */
typedef struct node
{
    double x;
    double y;
    double z;  /* position of the vertex */

    double d_x;
    double d_y;
    double d_z;   /* Recorded displacements */

    double dist_vert[3];	/* distance between the vertex and neighb[i] */

    int border; /*integer: 1 for on the border, 0 for negatory */

    struct node *pneighb[3];		/* pointer to the the neighbor i */
    struct dual *pncell[3];		/* pointer to the the neighboring cell i */
    struct bond *pnbond[3];      /* pointer to the neighboring bond i */
} Node;

/* Structure used to define the state of a cell */
struct strmarker {
    unsigned int size;
    unsigned int div_rate; /* type of division rate */
    unsigned int clone_index; /* index of a clone when studying several clones on the same lattice */
    double tension_index;
    unsigned int border;
};

typedef struct dual
{
   double area; /* Area of the cell */
   double area_soll; /* Area of the cell */
   double area_soll0; /* Natural area of the cell */
   double sqrtarea;
   double perimeter; /* perimeter of the cell */
   vector3 centroid; /* coordinates of the center of the cell */

   int nb_vertices; /* # of vertices in the cell */
   struct node **vertexlist; /* List of vertices in the cell */
   struct dual **celllist; /*list of neighboring cells*/
   struct bond **bondlist; /*list of bonds on the boundary*/
   struct strmarker marker;
} Dual;

extern unsigned int edges_dist[10];

extern unsigned int rand_disp;

/* pointer to the array of pointers to vertices */
extern Node **web;
/* pointer to the array of pointers to cells */
extern Dual **web_dual;
/* pointer to the array of pointers to bonds */
extern Bond **web_bond;

/* Array of division rate values (associated to marker.div_rate)
 * It is first initialised in lattice.c then reinitialised when cells
 * differentiate */
extern double division_rates[2];

/* Array of cell sizes values (associated to marker.size)
 * It is first initialised in lattice.c then reinitialised when cells
 * differentiate */
extern double cell_sizes[2];


/* Time on the lattice (number of displacement steps) */
extern long int itime;
extern double ttime;

/* Number of T1 processes */
extern int it1process;

/* Allocation of memory for the structures on the lattice */

/* Allocation of an array of pointers to vertices */
Node **createnodevector(int dimension);
/* Reallocation of an array of pointers to vertices */
Node **reallocnodevector(Node**vectortorealloc, int newdimension);

/* Allocation of an array of pointers to cells */
Dual **createdualvector(int dimension);
/* Reallocation of an array of pointers to cells */
Dual **reallocdualvector(Dual **vectortorealloc, int newdimension);

Bond **createbondvector(int dimension);
Bond **reallocbondvector(Bond **vectortorealloc, int newdimension);

/* Allocation of a vertex */
Node *allocnode();
/* Allocation of a cell */
Dual *allocdual();

Bond *allocbond();

/* Allocation of a list of vertices */
Node **allocvertexlist(int dimension);
/* Rellocation of a list of vertices */
Node ** reallocvertexlist(Node **listtorealloc, int newdimension);

Bond **allocbondlist(int dimension);
Bond **reallocbondlist(Bond **listtorealloc, int newdimension);

Dual **alloccelllist(int dimension);
Dual **realloccelllist(Dual **listtorealloc, int newdimension);

/* See document file */
int vertex_from_coord(int row, int column);

void init_lattice();

void kill_lattice();

void coord_gen();

void shrink_datafirst();

void shrink_data();

void cell0gen(Node *pfirstvertex);

void circlize(double radius, Node *pcenter);
    
void random_disp(double max_disp);

int indexcell(Node *pvertex, Dual *pcell);

void cell_update();

vector3 centroid(Dual *pcell);

vector3 centroid2(Dual *pcell);

double dist(Node *pvertex1, Node *pvertex2);

double sumtheta(Dual *pcell);

void disp(Node *pvertex);

int indexvertex(Node *pvertex, Dual *pcell);

int indexvertexinweb(Node *pvertex);

int indexcellinwebdual(Dual *pcell);

void insvertex(Dual *pcell, Node *pnewvert, Node *pnextvert);

void delvertex(Dual *pcell, Node *pdelvert);

int kill_cell(Dual *pcell, int ivertex);

int kill_cell_p(Dual *pcell);

void division_cell(Dual *pcell, int ivertex);

void division(Node *pvertex, int icell);
void division_eq(Node *pvertex, int icell);

void t1transform(Node *pvertex, int neighbor);

double area(Dual *pcell);

void stat_edges();

int kill_cell_area(Dual *pcell);

#endif /* SRC_LATTICE_LATTICE_H */

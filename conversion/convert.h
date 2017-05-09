typedef struct {double c[3][3];} matrix33;

typedef struct node
{
  double x;
  double y;             /* position of the vertex */
  double z;

  int move;             /* movement flag */

  double dist_vert[3];  /* distance between the vertex and neighb[i] */

  int border; /*integer: 1 for on the border, 0 for not*/
  
  struct node *pneighb[3];            /* pointer to the the neighbor i */
  struct dual *pncell[3];             /* pointer to the the neighborhood cell i */
  struct bond *pnbond[3];	      /* pointer to the neighboring bond i */
} Node;

typedef struct strmarker 
{
  unsigned int size;
  unsigned int div_rate; /* type of division rate */
  unsigned int clone_index; /* index of a clone when studying several clones on the same lattice */
  unsigned int daughter;
  unsigned int difftype;
  unsigned int tokill; /* set to one when killing the cell */
  double tension_index;
  unsigned int border;
} strmarker;

typedef struct
{
    double x;
    double y;
} vector2;

typedef struct
{
    double x;
    double y;
    double z;
} vector3;

typedef struct dual
{   
  double thick;
  double thick_0;
  double area; /* Area of the cell */
  double area_soll; /* Area of the cell */
  double area_soll0;
  double vol_soll;
  double sqrtarea;
  vector3 centroid; /* coordinates of the center of the cell */
  double perimeter; /* perimeter of the cell */
  int nb_vertices; /* # of vertices in the cell */
  double randfx, randfy;
  double ampli;
  struct node **vertexlist; /* List of vertices in the cell */
  struct dual **celllist;
  struct bond **bondlist;
  struct strmarker marker;
} Dual;

typedef struct bond
{
  struct node *pnvert[2]; //pointers to vertex endpoints
  struct dual *pncell[4]; //pointers to neighboring cells
  struct bond *pnbond[4]; //pointers to neighboring bonds
} Bond;

extern int nb_bond_tot,nb_cell_tot,nb_vertex_tot;
extern Node **web;
extern Dual **web_dual;
extern Bond **web_bond;


vector3 centroid(Dual *pcell);
Dual **createdualvector(int dimension);
Node *allocnode();
Node **creatnodevector(int dimension);
int converta();

/* Allocation of memory for the structures on the lattice */
/* Allocation of an array of pointers to vertices */
Node **createnodevector(int dimension);
/* Allocation of an array of pointers to cells */
Dual **createdualvector(int dimension);
/* Allocation of an array of pointers to bonds */
Bond **createbondvector(int dimension);
/* Allocation of a vertex */
Node *allocnode();
/* Allocation of a cell */
Dual *allocdual();
/* Allocation of a bond */
Bond *allocbond();
/* Allocation of a list of vertices */
Node **allocvertexlist(int dimension);
/* Allocation of a list of cells */
Dual **alloccelllist(int dimension);
/* Allocation of a list of bonds */
Bond **allocbondlist(int dimension);

double dist(Node *pvertex1, Node *pvertex2);
double perimeter(Dual *pcell);

int indexvertexinweb(Node *pvertex);
int indexcellinwebdual(Dual *pcell);

void freecell(Dual *pcell);
void kill_lattice();
void shrink_data();

int indexvertex(Node *pvertex, Dual *pcell);

//to help keep track of force balance
extern double **forcevector;
//different types of overall area calculations.  both call area_calc() below
double area(Dual * pcell);
double part_area(Dual *pcell,int max,int ivert0cell0,double x,double y,double z);

//calculates the area of the triangle created by pv1 and pv2 about centroid
vector3 area_calc(Node *pv1, Node *pv2, vector3 centroid);
//finds the bond index in web_bond
int indexbondinwebbond(Bond *pbond);
//prints files that are useful for kevin's opengl vis program
void latticeprint(char * filename, unsigned int type);
//prints kevin-patented *.*info files
void infoprint(char * filename, unsigned int type);

void countbonds();//counts the # of bonds
void initbonds();//initializes the bond structure
void nodeneighb();//finds neighbors to cells
void bondneighb();//finds neighbors to bonds
void neighbcell();//finds neighboring cells to cells
void updatedata();//calls all the previous functions
void cleanup();//cleans up the new data structures

void coord_mod();//modifies coords by an offset in z

void border_check();//checks to see which vertices are on the border for opti

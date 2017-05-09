void epsfile_pressure(char *filename, unsigned int type);
void epsfile_pressure_dpp(char *filename);

void epsfile(char *filename, unsigned int type);
void print_cell(Dual *pcell, FILE *pfile, unsigned int type);
void print_cell_dpp(Dual *pcell, double avg_press, FILE *pfile);


void printdefects(Dual *pcell, FILE *pfile);
void printpressure(Dual *pcell, FILE *pfile);
void printthick(Dual *pcell, FILE *pfile);
void printbound(Dual *pcell, FILE *pfile);
void printgrad(FILE *pfile);
void svglattice(char *filename, unsigned int type);
void svgprcell(Dual *pcell, FILE *pfile, unsigned int type);
void printgradcell(Dual *pcell, FILE *pfile);
void printampli(Dual *pcell, FILE *pfile);

#define DISC_RADIUS 60

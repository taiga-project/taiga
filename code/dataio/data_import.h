int read_vector(double **name, char *folder, char *shotname, char *filename, bool warning_on);
int read_vector(double **name, char *folder, char *shotname, char *filename);
int read_vector(double **name, char *folder, char *shotname, char *filename, int *successful);
int matrixColoumnReader(double **name, char *folder, char *shotname, char *filename, int coloumn_id, int total_coloumn);
int matrixColoumnReader(double **name, char *path, int coloumn_id, int total_coloumn);

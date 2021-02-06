void export_table(char *folder, char *runnumber, char *filename, int Ndat, 
                  double *dat1, char *head1, double *dat2, char *head2, double *dat3, char *head3, double *dat4, char *head4, double *dat5, char *head5, double *dat6, char *head6 );

void export_data(double *dat, int Ndat, char *folder, char *runnumber, char *filename);
void export_data(double *dat, int Ndat, char *folder, char *runnumber, char *subdir, char *filename);
void export_data(double *dat, int Ndat, char *folder, char *runnumber, char *filename, int dat_per_line);
void export_data(double *dat, int Ndat, char *folder, char *runnumber, char *subdir, char *filename, int dat_per_line);

void export_data(int *dat, int Ndat, char *folder, char *runnumber, char *filename);
void export_data(int *dat, int Ndat, char *folder, char *runnumber, char *subdir, char *filename);
void export_data(int *dat, int Ndat, char *folder, char *runnumber, char *filename, int dat_per_line);
void export_data(int *dat, int Ndat, char *folder, char *runnumber, char *subdir, char *filename, int dat_per_line);

void export_data(long *dat, int Ndat, char *folder, char *runnumber, char *filename);
void export_data(long *dat, int Ndat, char *folder, char *runnumber, char *subdir, char *filename);
void export_data(long *dat, int Ndat, char *folder, char *runnumber, char *filename, int dat_per_line);
void export_data(long *dat, int Ndat, char *folder, char *runnumber, char *subdir, char *filename, int dat_per_line);

void export_header(char *text, char *folder, char *runnumber);
void export_header(char *dataname,char *unitname,double dat, char *folder, char *runnumber);
void export_header(char *dataname, char *unitname, double dat, double dat2, char *folder, char *runnumber);
void export_header_addline(char *folder, char *runnumber);

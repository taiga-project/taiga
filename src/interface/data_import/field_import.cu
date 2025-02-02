int spline_read_and_init(ShotProp shot, RunProp run, char* field_name, double ***return_s_ptr, int dimRZ){
    char *spline_folder = "input/fieldSpl";
    int suc[1] = {1};
    
    double *S0,  *s0;  read_vector(&S0, spline_folder, shot.name, concat(field_name ,".spl11", NULL), suc);    cudaMalloc((void **) &s0,  dimRZ);
    double *S1,  *s1;  read_vector(&S1, spline_folder, shot.name, concat(field_name ,".spl12", NULL), suc);    cudaMalloc((void **) &s1,  dimRZ);
    double *S2,  *s2;  read_vector(&S2, spline_folder, shot.name, concat(field_name ,".spl13", NULL), suc);    cudaMalloc((void **) &s2,  dimRZ);
    double *S3,  *s3;  read_vector(&S3, spline_folder, shot.name, concat(field_name ,".spl14", NULL), suc);    cudaMalloc((void **) &s3,  dimRZ);
    double *S4,  *s4;  read_vector(&S4, spline_folder, shot.name, concat(field_name ,".spl21", NULL), suc);    cudaMalloc((void **) &s4,  dimRZ);
    double *S5,  *s5;  read_vector(&S5, spline_folder, shot.name, concat(field_name ,".spl22", NULL), suc);    cudaMalloc((void **) &s5,  dimRZ);
    double *S6,  *s6;  read_vector(&S6, spline_folder, shot.name, concat(field_name ,".spl23", NULL), suc);    cudaMalloc((void **) &s6,  dimRZ);
    double *S7,  *s7;  read_vector(&S7, spline_folder, shot.name, concat(field_name ,".spl24", NULL), suc);    cudaMalloc((void **) &s7,  dimRZ);
    double *S8,  *s8;  read_vector(&S8, spline_folder, shot.name, concat(field_name ,".spl31", NULL), suc);    cudaMalloc((void **) &s8,  dimRZ);
    double *S9,  *s9;  read_vector(&S9, spline_folder, shot.name, concat(field_name ,".spl32", NULL), suc);    cudaMalloc((void **) &s9,  dimRZ);
    double *S10, *s10; read_vector(&S10,spline_folder, shot.name, concat(field_name ,".spl33", NULL), suc);    cudaMalloc((void **) &s10,  dimRZ);
    double *S11, *s11; read_vector(&S11,spline_folder, shot.name, concat(field_name ,".spl34", NULL), suc);    cudaMalloc((void **) &s11,  dimRZ);
    double *S12, *s12; read_vector(&S12,spline_folder, shot.name, concat(field_name ,".spl41", NULL), suc);    cudaMalloc((void **) &s12,  dimRZ);
    double *S13, *s13; read_vector(&S13,spline_folder, shot.name, concat(field_name ,".spl42", NULL), suc);    cudaMalloc((void **) &s13,  dimRZ);
    double *S14, *s14; read_vector(&S14,spline_folder, shot.name, concat(field_name ,".spl43", NULL), suc);    cudaMalloc((void **) &s14,  dimRZ);
    double *S15, *s15; read_vector(&S15,spline_folder, shot.name, concat(field_name ,".spl44", NULL), suc);    cudaMalloc((void **) &s15,  dimRZ);
    
    size_t dimB = 16*sizeof(double*);
    double *S_PTR[16];  double **s_ptr; cudaMalloc((void **) &s_ptr,  dimB);
    
    S_PTR[0]  = s0;     S_PTR[1]  = s1 ;    S_PTR[2]  = s2;     S_PTR[3]  = s3;
    S_PTR[4]  = s4;     S_PTR[5]  = s5;     S_PTR[6]  = s6;     S_PTR[7]  = s7;
    S_PTR[8]  = s8;     S_PTR[9]  = s9;     S_PTR[10] = s10;    S_PTR[11] = s11;
    S_PTR[12] = s12;    S_PTR[13] = s13;    S_PTR[14] = s14;    S_PTR[15] = s15;
    
    if (suc[0] == 1){
        cudaMemcpy(s0, S0, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s1, S1, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s2, S2, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s3, S3, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s4, S4, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s5, S5, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s6, S6, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s7, S7, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s8, S8, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s9, S9, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s10, S10, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s11, S11, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s12, S12, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s13, S13, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s14, S14, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s15, S15, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s_ptr, S_PTR, dimB, cudaMemcpyHostToDevice);
        free(S0);   free(S1);   free(S2);   free(S3);
        free(S4);   free(S5);   free(S6);   free(S7);
        free(S8);   free(S9);   free(S10);  free(S11);
        free(S12);  free(S13);  free(S14);  free(S15);
        
        if (run.debug == 1){
            for (int i=0; i<10; ++i){
                printf("%s spline spl11 %d %lf\n", field_name, i, S0[i]);
            }
        }
        *return_s_ptr = s_ptr;
    }
    return suc[0];
}

int magnetic_field_read_and_init(ShotProp shot, RunProp run, TaigaCommons *s_host, TaigaCommons *s_shared){
    
    size_t dimB = 16*sizeof(double*);
    
    size_t dimRZ = s_host->grid_size[0]*s_host->grid_size[1]*sizeof(double);
    double **br_ptr;    cudaMalloc((void **) &br_ptr,  dimB);
    double **bt_ptr;    cudaMalloc((void **) &bt_ptr,  dimB);
    double **bz_ptr;    cudaMalloc((void **) &bz_ptr,  dimB);
    
    int s = 0;

    s = spline_read_and_init(shot, run, "brad", &br_ptr, dimRZ);
    spline_read_and_init(shot, run, "bz",   &bz_ptr, dimRZ);
    spline_read_and_init(shot, run, "btor", &bt_ptr, dimRZ);

    s_shared->brad = br_ptr;
    s_shared->bz   = bz_ptr;
    s_shared->btor = bt_ptr;
    
    if (s == 0){
        printf("Fatal error in the memory allocation of the magnetic field\n");
        exit(1);
    }
    
    return s;
}

int electric_field_read_and_init(ShotProp shot, RunProp run, TaigaCommons *s_host, TaigaCommons *s_shared){
    
    size_t dimB = 16*sizeof(double*);
    
    size_t dimRZ = s_host->grid_size[0]*s_host->grid_size[1]*sizeof(double);
    double **er_ptr;    cudaMalloc((void **) &er_ptr,  dimB);
    double **et_ptr;    cudaMalloc((void **) &et_ptr,  dimB);
    double **ez_ptr;    cudaMalloc((void **) &ez_ptr,  dimB);
    
    int is_electric_field_on = spline_read_and_init(shot, run, "erad", &er_ptr, dimRZ);
    spline_read_and_init(shot, run, "ez",   &ez_ptr, dimRZ);
    spline_read_and_init(shot, run, "etor", &et_ptr, dimRZ);
    
    s_shared->erad = er_ptr;
    s_shared->ez   = ez_ptr;
    s_shared->etor = et_ptr;
    s_shared->is_electric_field_on = is_electric_field_on;
    
    return is_electric_field_on;
}

int poloidal_flux_read_and_init(ShotProp shot, RunProp run, TaigaCommons *s_host, TaigaCommons *s_shared){

    size_t dimB = 16*sizeof(double*);

    size_t dimRZ = s_host->grid_size[0]*s_host->grid_size[1]*sizeof(double);
    double **psi_n;    cudaMalloc((void **) &psi_n, dimB);

    int is_magnetic_field_perturbation = spline_read_and_init(shot, run, "psi_n", &psi_n, dimRZ);

    s_shared->psi_n = psi_n;
    s_shared->is_magnetic_field_perturbation = is_magnetic_field_perturbation;

    return is_magnetic_field_perturbation;
}

//B-spline
int bspline_read_and_init(ShotProp shot, RunProp run, char* field_name, double ***return_s_ptr, int dimRZ){
    char *spline_folder = "input/fieldSpl";
    int suc[1] = {1};

    double *S0,  *s0;  read_vector(&S0, spline_folder, shot.name, concat(field_name ,".bspl", NULL), suc);
    cudaMalloc((void **) &s0,  dimRZ);

    size_t dimB = sizeof(double*);
    double *S_PTR[1];  double **s_ptr; cudaMalloc((void **) &s_ptr,  dimB);
    S_PTR[0]  = s0;

    if (suc[0] == 1){
        cudaMemcpy(s0, S0, dimRZ, cudaMemcpyHostToDevice);
        cudaMemcpy(s_ptr, S_PTR, dimB, cudaMemcpyHostToDevice);
        free(S0);
        if (run.debug == 1){
            for (int i=0; i<10; ++i){
                printf("%s bspline %d %lf\n", field_name, i, S0[i]);
            }
        }
        *return_s_ptr = s_ptr;
    }
    return suc[0];
}

int magnetic_field_read_and_init_with_bsplines(ShotProp shot, RunProp run, TaigaCommons *s_host, TaigaCommons *s_shared){

    size_t dimB = sizeof(double*);
    size_t dimRZ = s_host->grid_size[0]*s_host->grid_size[1]*sizeof(double);
    double **br_ptr;    cudaMalloc((void **) &br_ptr,  dimB);
    double **bt_ptr;    cudaMalloc((void **) &bt_ptr,  dimB);
    double **bz_ptr;    cudaMalloc((void **) &bz_ptr,  dimB);

    int s = 0;
    s = bspline_read_and_init(shot, run, "brad", &br_ptr, dimRZ);
    bspline_read_and_init(shot, run, "bz",   &bz_ptr, dimRZ);
    bspline_read_and_init(shot, run, "btor", &bt_ptr, dimRZ);

    s_shared->brad = br_ptr;
    s_shared->bz   = bz_ptr;
    s_shared->btor = bt_ptr;

    if (s == 0){
        printf("Fatal error in the memory allocation of the magnetic field\n");
        exit(1);
    }

    return s;
}

int poloidal_flux_read_and_init_with_bsplines(ShotProp shot, RunProp run, TaigaCommons *s_host, TaigaCommons *s_shared){

    size_t dimB = 16*sizeof(double*);

    size_t dimRZ = s_host->grid_size[0]*s_host->grid_size[1]*sizeof(double);
    double **psi_n;    cudaMalloc((void **) &psi_n, dimB);

    int is_magnetic_field_perturbation = bspline_read_and_init(shot, run, "psi_n", &psi_n, dimRZ);

    s_shared->psi_n = psi_n;
    s_shared->is_magnetic_field_perturbation = is_magnetic_field_perturbation;

    return is_magnetic_field_perturbation;
}
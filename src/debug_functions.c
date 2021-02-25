void debug_message_init(TaigaGlobals *g){
    printf("----------------------------------------------------------\n");
    printf("Starting coordinates\n");
    printf("X[0]:\t %lf\t %lf\t %lf\n",g->rad[0], g->z[0], g->tor[0]);
    printf("X[1]:\t %lf\t %lf\t %lf\n",g->rad[1], g->z[1], g->tor[1]);
    printf("V[0]:\t %lf\t %lf\t %lf\n",g->vrad[0], g->vz[0], g->vtor[0]);
    printf("----------------------------------------------------------\n");
    for(int i=1; i<20; ++i){
        printf("X[%2d]:\t %le\t %le\t %le\n", i, g->rad[i], g->z[i], g->tor[i]);
    }
    printf("----------------------------------------------------------\n");
}

void debug_message_run(TaigaGlobals *g){
    printf("----------------------------------------------------------\n");
    printf("Ending coordinates\n");
    printf("X[0]\t %lf\t %lf\t %lf\n", g->rad[0], g->z[0], g->tor[0]);
    printf("X[1]\t %lf\t %lf\t %lf\n", g->rad[1], g->z[1], g->tor[1]);
    printf("V[0]\t %lf\t %lf\t %lf\n", g->vrad[0], g->vz[0], g->vtor[0]);
    printf("----------------------------------------------------------\n");
}

void debug_service_vars(double *SERVICE_VAR){
    for(int i=0; i<SERVICE_VAR_LENGTH; ++i){
        printf("SERVICE VAR %d\t%lf\n", i, SERVICE_VAR[i]);
    }
}

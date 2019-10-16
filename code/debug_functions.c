void debug_message_init(double* XR, double* XZ, double* XT, double* VR, double* VZ, double* VT){
    printf("ionV:  0.\t %lf\t %lf\t %lf\n",VR[0], VZ[0], VT[0]);
    printf("ionX:  0.\t %lf\t %lf\t %lf\n",XR[0], XZ[0], XT[0]);
    printf("ionX:  1.\t %lf\t %lf\t %lf\n",XR[1], XZ[1], XT[1]);

    printf("----------------------------------------------------------\n");
    printf("ion:  0.\t %lf\t %lf\t %lf\n",XR[0], XZ[0], XT[0]);
    printf("----------------------------------------------------------\n");
    for(int i=1; i<20; i++){
        printf("ion: %2d.\t %le\t %le\t %le\n", i, XR[i], XZ[i], XT[i]);
    }
    printf("----------------------------------------------------------\n");  
}

void debug_message_run(double* XR, double* XZ, double* XT, double* VR, double* VZ, double* VT){
    printf("Xion:  0.\t %lf\t %lf\t %lf\n", XR[0], XZ[0], XT[0]);
    printf("Xion:  1.\t %lf\t %lf\t %lf\n", XR[1], XZ[1], XT[1]);
    printf("Vion:  0.\t %lf\t %lf\t %lf\n", VR[0], VZ[0], VT[0]);
}

void debug_service_vars(double *SERVICE_VAR){
    for(int i=0; i<SERVICE_VAR_LENGTH; i++){
        printf("SERVICE VAR %d\t%lf\n", i, SERVICE_VAR[i]);
    }
}

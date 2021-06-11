#include "lorentz.cu"
#include "localise_field.cuh"

//header
__device__ double get_dr_with_polynomials(TaigaCommons *c, int *local_spline_indices, double R);
__device__ double get_dz_with_polynomials(TaigaCommons *c, int *local_spline_indices, double Z);

__device__ void (*solve_diffeq)(double *X, double *a, double *B, double *E, double eperm, double timestep);

__device__ void (*get_coefficients)(TaigaCommons *c,
                                    int *local_spline_indices,
                                    double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                    double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                    double *local_polflux, int zgrid_length);

__device__ double (*calculate_local_field)(TaigaCommons *c, int *local_spline_indices,
                                           double *local_spline, double dr, double dz);

__device__ double (*get_dr)(TaigaCommons *c, int *local_spline_indices, double R);
__device__ double (*get_dz)(TaigaCommons *c, int *local_spline_indices, double Z);

__device__ double get_dr_with_polynomials(TaigaCommons *c, int *local_spline_indices, double R){
    return R - c->spline_rgrid[local_spline_indices[0]];
}
__device__ double get_dz_with_polynomials(TaigaCommons *c, int *local_spline_indices, double Z){
    return Z - c->spline_zgrid[local_spline_indices[1]];
}

__device__ int calculate_trajectory(TaigaCommons *c, double X[6], int detcellid){
    // next grid
    int local_spline_indices[2];
    local_spline_indices[0] = SPLINE_INDEX_ERROR;
    local_spline_indices[1] = SPLINE_INDEX_ERROR;

    double local_spline_brad[16];
    double local_spline_bz[16];
    double local_spline_btor[16];

    double local_spline_erad[16];
    double local_spline_ez[16];
    double local_spline_etor[16];

    double local_polflux[16];

    double local_bfield[3], local_efield[3];
    double dr, dz;
    double R;

    double eperm = c->eperm;
    double timestep = c->timestep;
    bool is_electric_field_on = c->is_electric_field_on;
    int magnetic_field_mode = c->magnetic_field_mode;

    double X_prev[6];
    double a[3] = {0, 0, 0};

    switch(c->solver){
        case SOLVER_RK45:
            solve_diffeq = &solve_diffeq_by_rk4;
            break;
        case SOLVER_VERLET:
            solve_diffeq = &solve_diffeq_by_verlet;
            break;
        case SOLVER_YOSHIDA:
            solve_diffeq = &solve_diffeq_by_yoshida;
            break;
    }
    int spline_mode = 0;
    switch(spline_mode){//(c->spline_mode){
        case 0:
            get_coefficients = &get_coefficients_with_polynomials;
            calculate_local_field = &calculate_local_field_with_polynomials;
            get_dr = &get_dr_with_polynomials;
            get_dz = &get_dz_with_polynomials;
            break;
    }

    if (is_electric_field_on){
        get_acceleration_from_lorentz_force = &get_acceleration_from_lorentz_force_with_electric_field;
    }else{
        get_acceleration_from_lorentz_force = &get_acceleration_from_lorentz_force_without_electric_field;
    }

    for (int loopi=0; (loopi < c->max_step_number && (detcellid == CALCULATION_NOT_FINISHED)); ++loopi){
        R = get_major_radius(X[0], X[2]);
        copy_local_field(c, R, X[1], local_spline_indices,
                            local_spline_brad, local_spline_bz, local_spline_btor,
                            local_spline_erad, local_spline_ez, local_spline_etor,
                            local_polflux);

        dr = (*get_dr)(c, local_spline_indices, R);
        dz = (*get_dz)(c, local_spline_indices, X[1]);

        local_bfield[0] = (*calculate_local_field)(c, local_spline_indices, local_spline_brad, dr, dz);
        local_bfield[1] = (*calculate_local_field)(c, local_spline_indices, local_spline_bz,   dr, dz);
        local_bfield[2] = (*calculate_local_field)(c, local_spline_indices, local_spline_btor, dr, dz);

        if (magnetic_field_mode == MAGNETIC_FIELD_FROM_FLUX){
            local_bfield[0] /= -X[0];    //Brad = -dPsi_dZ / R
            local_bfield[1] /=  X[0];    //Bz   =  dPsi_dR / R
            local_bfield[2] /=  X[0];    //Btor = (R*Btor) / R
        }

        local_bfield[0] = get_rad_from_poloidal(R, local_bfield[0], local_bfield[2], X[0], X[2]);
        local_bfield[2] = get_tor_from_poloidal(R, local_bfield[0], local_bfield[2], X[0], X[2]);

        if (is_electric_field_on){
            local_efield[0] = (*calculate_local_field)(c, local_spline_indices, local_spline_erad, dr, dz);
            local_efield[1] = (*calculate_local_field)(c, local_spline_indices, local_spline_ez,   dr, dz);
            local_efield[2] = (*calculate_local_field)(c, local_spline_indices, local_spline_etor, dr, dz);
            local_efield[0] = get_rad_from_poloidal(R, local_efield[0], local_efield[2], X[0], X[2]);
            local_efield[2] = get_tor_from_poloidal(R, local_efield[0], local_efield[2], X[0], X[2]);
        }

        // archive coordinates
        for(int i=0; i<6; ++i)  X_prev[i] = X[i];

        (*solve_diffeq)(X, a, local_bfield, local_efield, eperm, timestep);

        detcellid = calculate_detection_position(X, X_prev, c->detector_geometry);
    }

    return detcellid;
}

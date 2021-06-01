#define SPLINE_INDEX_ERROR -1

#include "lorentz.cu"

__device__ void copy_local_field(TaigaCommons *c,
                                 double position_rad, double position_z,
                                 int *local_spline_indices,
                                 double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                 double *local_spline_erad, double *local_spline_ez, double *local_spline_etor){
    int rci, zci;
    int i, i2;
    int rgrid_length = c->grid_size[0];
    int zgrid_length = c->grid_size[1];

    for(rci=0; (c->spline_rgrid[rci+1]<position_rad)&&(rci<rgrid_length-1); ++rci){;}

    for(zci=0; (c->spline_zgrid[zci+1]<position_z)&&(zci<zgrid_length-1); ++zci){;}

    // Particle leave out the cell
    if ((local_spline_indices[0] != rci) || (local_spline_indices[1] != zci)){
        local_spline_indices[0] = rci;
        local_spline_indices[1] = zci;

        for(i=0; i<16; ++i){
            i2 = (local_spline_indices[0])*(zgrid_length-1)+local_spline_indices[1];
            local_spline_brad[i] = c->brad[i][i2];
            local_spline_bz[i]   = c->bz[i][i2];
            local_spline_btor[i] = c->btor[i][i2];
        }
        if (c->is_electric_field_on){
            for(i=0; i<16; ++i){
                local_spline_erad[i] = c->erad[i][i2];
                local_spline_ez[i]   = c->ez[i][i2];
                local_spline_etor[i] = c->etor[i][i2];
            }
        }
    }
}

__device__ double calculate_local_field(double *local_spline, double dr, double dz){
    /* MATLAB CODE:
    sample2(3) =c11(bs1,bs2)*dsx^3*dsy^3 + c12(bs1,bs2)*dsx^3*dsy^2 + c13(bs1,bs2)*dsx^3*dsy + c14(bs1,bs2)*dsx^3 + ...
                c21(bs1,bs2)*dsx^2*dsy^3 + c22(bs1,bs2)*dsx^3*dsy^2 + c23(bs1,bs2)*dsx^2*dsy + c24(bs1,bs2)*dsx^2 + ...
                c31(bs1,bs2)*dsx  *dsy^3 + c32(bs1,bs2)*dsx  *dsy^2 + c33(bs1,bs2)*dsx  *dsy + c34(bs1,bs2)*dsx    + ...
                c41(bs1,bs2)      *dsy^3 + c42(bs1,bs2)      *dsy^2 + c43(bs1,bs2)      *dsy + c44(bs1,bs2);*/
    double local_field = 0.0, local_field_comp[16] ;
    for(int i=0; i<4; ++i){
        for(int j=0; j<4; ++j){
            local_field_comp[i*4+j] = local_spline[i*4+j]*pow(dr,3-i)*pow(dz,3-j);
        }
    }

    for(int i=0; i<4; ++i){
        for(int j=0; j<4; ++j){
            local_field += local_field_comp[i*4+j];
        }
    }
    return local_field;
}

__device__ void (*solve_diffeq)(double *X, double *a, double *B, double *E, double eperm, double timestep);

__device__ int traj(TaigaCommons *c, double X[6], int detcellid){
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

    if (is_electric_field_on){
        get_acceleration_from_lorentz_force = &get_acceleration_from_lorentz_force_with_electric_field;
    }else{
        get_acceleration_from_lorentz_force = &get_acceleration_from_lorentz_force_without_electric_field;
    }

    for (int loopi=0; (loopi < c->max_step_number && (detcellid == CALCULATION_NOT_FINISHED)); ++loopi){
        R = get_major_radius(X[0], X[2]);
        copy_local_field(c, R, X[1], local_spline_indices,
                         local_spline_brad, local_spline_bz, local_spline_btor,
                         local_spline_erad, local_spline_ez, local_spline_etor);

        dr = R-c->spline_rgrid[local_spline_indices[0]];
        dz = X[1]-c->spline_zgrid[local_spline_indices[1]];

        local_bfield[0] = calculate_local_field(local_spline_brad, dr, dz);
        local_bfield[1] = calculate_local_field(local_spline_bz,   dr, dz);
        local_bfield[2] = calculate_local_field(local_spline_btor, dr, dz);

        if (magnetic_field_mode == MAGNETIC_FIELD_FROM_FLUX){
            local_bfield[0] /= -X[0];    //Brad = -dPsi_dZ / R
            local_bfield[1] /=  X[0];    //Bz   =  dPsi_dR / R
            local_bfield[2] /=  X[0];    //Btor = (R*Btor) / R
        }

        local_bfield[0] = get_rad_from_poloidal(R, local_bfield[0], local_bfield[2], X[0], X[2]);
        local_bfield[2] = get_tor_from_poloidal(R, local_bfield[0], local_bfield[2], X[0], X[2]);

        if (is_electric_field_on){
            local_efield[0] = calculate_local_field(local_spline_erad, dr, dz);
            local_efield[1] = calculate_local_field(local_spline_ez,   dr, dz);
            local_efield[2] = calculate_local_field(local_spline_etor, dr, dz);
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

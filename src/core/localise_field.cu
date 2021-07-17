#include "localise_field.cuh"
#include "bspline.cuh"

__device__ void (*get_coefficients)(TaigaCommons *c,
                                    int *local_spline_indices,
                                    double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                    double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                    double *local_psi_n, int rgrid_length, int zgrid_length);

__device__ void get_coefficients_with_bsplines(
                                TaigaCommons *c,
                                int *local_spline_indices,
                                double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                double *local_psi_n, int rgrid_length, int zgrid_length){
    const int k = 3;
    const int k_plus_1 = 4;
    int i, j, index;

    if (local_spline_indices[0] < k_plus_1){
        local_spline_indices[0] = k;
    }
    if (local_spline_indices[1] < k_plus_1){
        local_spline_indices[1] = k;
    }
    for(i=0; i<k_plus_1; ++i){
        for(j=0; j<k_plus_1; ++j){
            index = (local_spline_indices[1]-k+j)*(rgrid_length-k_plus_1)
                    +local_spline_indices[0]-k+i;
            local_spline_brad[i*k_plus_1+j] = c->brad[0][index];
            local_spline_bz[i*k_plus_1+j]   = c->bz[0][index];
            local_spline_btor[i*k_plus_1+j] = c->btor[0][index];
            if (c->is_electric_field_on) {
                local_spline_erad[i*k_plus_1+j] = c->erad[0][index];
                local_spline_ez[i*k_plus_1+j]   = c->ez[0][index];
                local_spline_etor[i*k_plus_1+j] = c->etor[0][index];
            }
            if (c->is_magnetic_field_perturbation) {
                local_psi_n[i * k_plus_1 + j] = c->psi_n[0][index];
            }
        }
    }
}

__device__ void get_coefficients_with_splines(
        TaigaCommons *c,
        int *local_spline_indices,
        double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
        double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
        double *local_psi_n, int rgrid_length, int zgrid_length){
    int i, i2;
    for(i=0; i<16; ++i){
        i2 = (local_spline_indices[0])*(zgrid_length-1)+local_spline_indices[1];
        local_spline_brad[i] = c->brad[i][i2];
        local_spline_bz[i]   = c->bz[i][i2];
        local_spline_btor[i] = c->btor[i][i2];
        if (c->is_electric_field_on){
            local_spline_erad[i] = c->erad[i][i2];
            local_spline_ez[i]   = c->ez[i][i2];
            local_spline_etor[i] = c->etor[i][i2];
        }
        if (c->is_magnetic_field_perturbation){
            local_psi_n[i] = c->psi_n[i][i2];
        }
    }
}

__device__ void copy_local_field_coefficients(TaigaCommons *c,
                                              double position_rad, double position_z,
                                              int *local_spline_indices,
                                              double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                              double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                              double *local_psi_n){
    int rci, zci;
    int rgrid_length = c->grid_size[0];
    int zgrid_length = c->grid_size[1];

    for(rci=0; (c->spline_rgrid[rci+1]<position_rad)&&(rci<rgrid_length-1); ++rci){;}
    for(zci=0; (c->spline_zgrid[zci+1]<position_z)&&(zci<zgrid_length-1); ++zci){;}

    // Particle leave out the cell
    if ((local_spline_indices[0] != rci) || (local_spline_indices[1] != zci)){
        local_spline_indices[0] = rci;
        local_spline_indices[1] = zci;

        (*get_coefficients)(c, local_spline_indices,
                            local_spline_brad, local_spline_bz, local_spline_btor,
                            local_spline_erad, local_spline_ez, local_spline_etor,
                            local_psi_n, rgrid_length, zgrid_length);
    }
}

__device__ void get_local_field(double *local_bfield, double *local_efield, bool is_electric_field_on,
                                TaigaCommons *c,
                                double R, double *X,
                                int *local_spline_indices,
                                double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                double *local_psi_n){

    double dr, dz;
    dr = (*get_dr)(c, local_spline_indices, R);
    dz = (*get_dz)(c, local_spline_indices, X[1]);

    local_bfield[0] = (*calculate_local_field)(c, local_spline_indices, local_spline_brad, dr, dz);
    local_bfield[1] = (*calculate_local_field)(c, local_spline_indices, local_spline_bz,   dr, dz);
    local_bfield[2] = (*calculate_local_field)(c, local_spline_indices, local_spline_btor, dr, dz);

    local_bfield[0] = get_rad_from_poloidal(R, local_bfield[0], local_bfield[2], X[0], X[2]);
    local_bfield[2] = get_tor_from_poloidal(R, local_bfield[0], local_bfield[2], X[0], X[2]);

    if (is_electric_field_on){
        local_efield[0] = (*calculate_local_field)(c, local_spline_indices, local_spline_erad, dr, dz);
        local_efield[1] = (*calculate_local_field)(c, local_spline_indices, local_spline_ez,   dr, dz);
        local_efield[2] = (*calculate_local_field)(c, local_spline_indices, local_spline_etor, dr, dz);
        local_efield[0] = get_rad_from_poloidal(R, local_efield[0], local_efield[2], X[0], X[2]);
        local_efield[2] = get_tor_from_poloidal(R, local_efield[0], local_efield[2], X[0], X[2]);
    }
}

__device__ double calculate_local_field_with_splines(TaigaCommons *c, const int *local_spline_indices,
                                                     const double *local_spline, double dr, double dz){
    /* MATLAB CODE:
    sample2(3) =c11(bs1,bs2)*dsx^3*dsy^3 + c12(bs1,bs2)*dsx^3*dsy^2 + c13(bs1,bs2)*dsx^3*dsy + c14(bs1,bs2)*dsx^3 + ...
                c21(bs1,bs2)*dsx^2*dsy^3 + c22(bs1,bs2)*dsx^3*dsy^2 + c23(bs1,bs2)*dsx^2*dsy + c24(bs1,bs2)*dsx^2 + ...
                c31(bs1,bs2)*dsx  *dsy^3 + c32(bs1,bs2)*dsx  *dsy^2 + c33(bs1,bs2)*dsx  *dsy + c34(bs1,bs2)*dsx    + ...
                c41(bs1,bs2)      *dsy^3 + c42(bs1,bs2)      *dsy^2 + c43(bs1,bs2)      *dsy + c44(bs1,bs2);*/
    double local_field = 0.0, local_field_comp[16];
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

__device__ double calculate_local_field_with_bsplines(TaigaCommons *c, const int *local_spline_indices,
                                                      const double *local_spline, double R, double Z){
    const int k = 3;
    const int k_plus_1 = 4;
    double B_R[k_plus_1];
    double B_Z[k_plus_1];
    bspline(B_R, R, k, local_spline_indices[0], c->spline_rgrid);
    bspline(B_Z, Z, k, local_spline_indices[1], c->spline_zgrid);
    double local_field = 0.0;
    for(int i=0; i<k_plus_1; ++i){
        for(int j=0; j<k_plus_1; ++j){
            local_field += local_spline[i*k_plus_1+j] * B_R[i] * B_Z[j];
        }
    }
    return local_field;
}

__device__ double get_dr_with_splines(TaigaCommons *c, const int *local_spline_indices, double R){
    return R - c->spline_rgrid[local_spline_indices[0]];
}
__device__ double get_dz_with_splines(TaigaCommons *c, const int *local_spline_indices, double Z){
    return Z - c->spline_zgrid[local_spline_indices[1]];
}

__device__ double get_dr_with_bsplines(TaigaCommons *c, const int *local_spline_indices, double R){
    return R;
}
__device__ double get_dz_with_bsplines(TaigaCommons *c, const int *local_spline_indices, double Z){
    return Z;
}

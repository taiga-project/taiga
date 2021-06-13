#include "localise_field.cuh"
#include "bspline.cu"

__device__ void (*get_coefficients)(TaigaCommons *c,
                                    const int *local_spline_indices,
                                    double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                    double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                    double *local_polflux, int zgrid_length);

__device__ void get_coefficients_with_bsplines(
                                TaigaCommons *c,
                                const int *local_spline_indices,
                                double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                double *local_polflux, int zgrid_length){
    const int k = 3;
    int i, j, index;
    for(i=local_spline_indices[0]-k; i<=local_spline_indices[0]; ++i){
        for(j=local_spline_indices[1]-k; j<=local_spline_indices[1]; ++j){
            index = i*(zgrid_length-1)+j;
            local_spline_brad[i] = c->brad[0][index];
            local_spline_bz[i] = c->bz[0][index];
            local_spline_btor[i] = c->btor[0][index];
            if (c->is_electric_field_on) {
                local_spline_erad[i] = c->erad[0][index];
                local_spline_ez[i] = c->ez[0][index];
                local_spline_etor[i] = c->etor[0][index];
            }
            if (c->is_magnetic_field_perturbation) {
                local_polflux[i] = c->polflux[0][index];
            }
        }
    }
}

__device__ void get_coefficients_with_splines(
        TaigaCommons *c,
        const int *local_spline_indices,
        double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
        double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
        double *local_polflux, int zgrid_length){
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
            local_polflux[i] = c->polflux[i][i2];
        }
    }
}

__device__ void copy_local_field(TaigaCommons *c,
                                 double position_rad, double position_z,
                                 int *local_spline_indices,
                                 double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                 double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                 double *local_polflux){
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
                                         local_polflux, zgrid_length);
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
    int local_index_r = local_spline_indices[0];
    int local_index_z = local_spline_indices[1];
    double *spline_r_grid = c->spline_rgrid;
    double *spline_z_grid = c->spline_zgrid;
    double local_field = 0.0;
    for(int i=0; i<4; ++i){
        int r_grid_index = local_index_r-k+i;
        for(int j=0; j<4; ++j){
            int z_grid_index = local_index_z-k+j;
            local_field += local_spline[i*4+j] *
                    B(R, k, r_grid_index, spline_r_grid) *
                    B(Z, k, z_grid_index, spline_z_grid);
        }
    }
    return local_field;
}

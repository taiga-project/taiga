#include "lorentz.cu"
#include "localise_field.cuh"

__device__ void (*solve_diffeq)(double *X, double *B, double *E, double *E_prev, double eperm, double timestep);

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

    double local_psi_n[16];

    double local_bfield[9], local_efield[9];
    double local_efield_prev[3]={0, 0, 0};

    double eperm = c->eperm;
    double timestep = c->timestep;
    bool is_electric_field_on = c->is_electric_field_on;
    long max_step_number = c->max_step_number;

    double X_prev[6];

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

    switch(c->field_interpolation_method){
        case CUBIC_SPLINE:
            get_coefficients = &get_coefficients_with_splines;
            calculate_local_field = &calculate_local_field_with_splines;
            get_dr = &get_dr_with_splines;
            get_dz = &get_dz_with_splines;
            break;
        case CUBIC_BSPLINE:
            get_coefficients = &get_coefficients_with_bsplines;
            calculate_local_field = &calculate_local_field_with_bsplines;
            get_dr = &get_dr_with_bsplines;
            get_dz = &get_dz_with_bsplines;
    }

    if (is_electric_field_on){
        get_acceleration_from_lorentz_force = &get_acceleration_from_lorentz_force_with_electric_field;
    }else{
        get_acceleration_from_lorentz_force = &get_acceleration_from_lorentz_force_without_electric_field;
    }

    for (int step_counter=0; (step_counter < max_step_number && (detcellid == CALCULATION_NOT_FINISHED)); ++step_counter){
        for(int i=0; i<6; ++i)  X_prev[i] = X[i];

        get_local_field(X, local_bfield, local_efield, c, is_electric_field_on,
                        local_spline_indices,
                        local_spline_brad, local_spline_bz, local_spline_btor,
                        local_spline_erad, local_spline_ez, local_spline_etor,
                        local_psi_n);

        (*solve_diffeq)(X, local_bfield, local_efield, local_efield_prev, eperm, timestep);

        detcellid = calculate_detection_position(X, X_prev, c->detector_geometry);
    }

    return detcellid;
}

#include "core/policies.cuh"
#include "core/physics/lorentz.cu"
#include "core/localise_field.cuh"
#include "core/physics/ionisation.cuh"

template<typename Solver, typename MagneticField,
         typename SecondaryIonisation, typename Lorentz,
         typename DetectInterp,
         typename Perturb>
__device__ int calculate_trajectory(TaigaCommons *c, double X[X_SIZE], int detcellid) {
    int   local_spline_indices[2] = { SPLINE_INDEX_ERROR, SPLINE_INDEX_ERROR };
    int   local_ts_index[1]       = { SPLINE_INDEX_ERROR };
    double local_spline_brad[16], local_spline_bz[16], local_spline_btor[16];
    double local_spline_erad[16], local_spline_ez[16], local_spline_etor[16];
    double local_spline_psi_n[16];
    double local_ts_psi[2]        = { UNDEFINED_FLOAT, UNDEFINED_FLOAT };

    const double eperm    = c->eperm;
    const double timestep = c->timestep;
    const long   max_step = c->max_step_number;

    double X_prev[X_SIZE];
    double psi_n;

    for (long step = 0; (step < max_step) && (detcellid == CALCULATION_NOT_FINISHED); ++step) {

#pragma unroll
        for (int i = 0; i < X_SIZE; ++i) {
            X_prev[i] = X[i];
        }

        MagneticField::get_coefficients(c, X, local_spline_indices);
        MagneticField::calculate_local_field(c, X, local_spline_indices,
                                             local_spline_brad, local_spline_bz, local_spline_btor,
                                             local_spline_erad, local_spline_ez, local_spline_etor,
                                             local_spline_psi_n);

        // Perturb::apply(X, c, timestep);

        psi_n = Solver::solve(X, eperm, timestep, c, EFieldSwitch::enabled,
                              local_spline_indices,
                              local_spline_brad, local_spline_bz, local_spline_btor,
                              local_spline_erad, local_spline_ez, local_spline_etor,
                              local_spline_psi_n);

        SecondaryIonisation::update(psi_n, X, c, local_ts_index, local_ts_psi, timestep);

        detcellid = DetectInterp::apply(X, X_prev, c->detector_geometry, timestep);
    }
    return detcellid;
}

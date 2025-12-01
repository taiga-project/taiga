#ifndef POLICIES_CUH
#define POLICIES_CUH

// Solver policies
struct SolverRK4 {
    __device__ __forceinline__
    static double solve(double *X, double eperm, double timestep,
                        TaigaCommons *c, bool is_electric_field_on,
                        int *local_spline_indices,
                        double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                        double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                        double *local_spline_psi_n)
    {
        return solve_diffeq_by_rk4(X, eperm, timestep, c, is_electric_field_on,
                                   local_spline_indices,
                                   local_spline_brad, local_spline_bz, local_spline_btor,
                                   local_spline_erad, local_spline_ez, local_spline_etor,
                                   local_spline_psi_n);
    }
};

struct SolverRKN {
    __device__ __forceinline__
    static double solve(double *X, double eperm, double timestep,
                        TaigaCommons *c, bool is_electric_field_on,
                        int *local_spline_indices,
                        double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                        double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                        double *local_spline_psi_n)
    {
        return solve_diffeq_by_rkn(X, eperm, timestep, c, is_electric_field_on,
                                   local_spline_indices,
                                   local_spline_brad, local_spline_bz, local_spline_btor,
                                   local_spline_erad, local_spline_ez, local_spline_etor,
                                   local_spline_psi_n);
    }
};

struct SolverVerlet {
    __device__ __forceinline__
    static double solve(double *X, double eperm, double timestep,
                        TaigaCommons *c, bool is_electric_field_on,
                        int *local_spline_indices,
                        double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                        double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                        double *local_spline_psi_n)
    {
        return solve_diffeq_by_verlet(X, eperm, timestep, c, is_electric_field_on,
                                      local_spline_indices,
                                      local_spline_brad, local_spline_bz, local_spline_btor,
                                      local_spline_erad, local_spline_ez, local_spline_etor,
                                      local_spline_psi_n);
    }
};

struct SolverYoshida {
    __device__ __forceinline__
    static double solve(double *X, double eperm, double timestep,
                        TaigaCommons *c, bool is_electric_field_on,
                        int *local_spline_indices,
                        double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                        double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                        double *local_spline_psi_n)
    {
        return solve_diffeq_by_yoshida(X, eperm, timestep, c, is_electric_field_on,
                                       local_spline_indices,
                                       local_spline_brad, local_spline_bz, local_spline_btor,
                                       local_spline_erad, local_spline_ez, local_spline_etor,
                                       local_spline_psi_n);
    }
};

// Field interpolation policies
struct CubicSplineInterp {
    __device__ __forceinline__
    static void get_coefficients(TaigaCommons *c, const double *X,
                                 int *local_spline_indices)
    {
        get_coefficients_with_splines(c, X, local_spline_indices);
    }

    __device__ __forceinline__
    static void calculate_local_field(TaigaCommons *c, const double *X,
                                      int *local_spline_indices,
                                      double *brad, double *bz, double *btor,
                                      double *erad, double *ez, double *etor,
                                      double *psi_n)
    {
        calculate_local_field_with_splines(c, X, local_spline_indices,
                                           brad, bz, btor,
                                           erad, ez, etor,
                                           psi_n);
    }

    __device__ __forceinline__
    static double get_dr(TaigaCommons *c, int *local_spline_indices) {
        return get_dr_with_splines(c, local_spline_indices);
    }

    __device__ __forceinline__
    static double get_dz(TaigaCommons *c, int *local_spline_indices) {
        return get_dz_with_splines(c, local_spline_indices);
    }
};

struct CubicBSplineInterp {
    __device__ __forceinline__
    static void get_coefficients(TaigaCommons *c, const double *X,
                                 int *local_spline_indices)
    {
        get_coefficients_with_bsplines(c, X, local_spline_indices);
    }

    __device__ __forceinline__
    static void calculate_local_field(TaigaCommons *c, const double *X,
                                      int *local_spline_indices,
                                      double *brad, double *bz, double *btor,
                                      double *erad, double *ez, double *etor,
                                      double *psi_n)
    {
        calculate_local_field_with_bsplines(c, X, local_spline_indices,
                                            brad, bz, btor,
                                            erad, ez, etor,
                                            psi_n);
    }

    __device__ __forceinline__
    static double get_dr(TaigaCommons *c, int *local_spline_indices) {
        return get_dr_with_bsplines(c, local_spline_indices);
    }

    __device__ __forceinline__
    static double get_dz(TaigaCommons *c, int *local_spline_indices) {
        return get_dz_with_bsplines(c, local_spline_indices);
    }
};

// Secondary ionisation policies
struct SecondaryIonisationOn {
    __device__ __forceinline__
    static void update(double psi_n, double *X, TaigaCommons *c,
                       int *local_ts_index, double *local_ts_psi, double timestep)
    {
        calculate_ionisation_loss(psi_n, X, c, local_ts_index, local_ts_psi, timestep);
    }
};

struct SecondaryIonisationOff {
    __device__ __forceinline__
    static void update(double psi_n, double *X, TaigaCommons *c,
                       int *local_ts_index, double *local_ts_psi, double timestep)
    {
        no_ionisation_loss(psi_n, X, c, local_ts_index, local_ts_psi, timestep);
    }
};

// Electric field policies
struct LorentzWithElectric {
    static constexpr bool enabled = true;
    __device__ __forceinline__
    static void acceleration(double *a,
                             const double *v,
                             const double *B,
                             const double *E,
                             double eperm) {
        get_acceleration_from_lorentz_force_with_electric_field(a, v, B, E, eperm);
    }
};

struct LorentzWithoutElectric {
    static constexpr bool enabled = false;
    __device__ __forceinline__
    static void acceleration(double *a,
                             const double *v,
                             const double *B,
                             const double * /*E*/,
                             double eperm) {
        get_acceleration_from_lorentz_force_without_electric_field(a, v, B, eperm);
    }
};

// Field perturbation policies
struct PerturbationOn {
    __device__ __forceinline__
    static void apply(double *X, TaigaCommons *c, double timestep) {
        // apply_field_perturbation(X, c, timestep);
    }
};

struct PerturbationOff {
    __device__ __forceinline__
    static void apply(double *X, TaigaCommons *c, double timestep) {
    }
};



#endif //POLICIES_CUH

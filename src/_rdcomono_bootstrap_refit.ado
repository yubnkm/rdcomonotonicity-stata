*! version 0.1.0 11aug2026
program define _rdcomono_bootstrap_refit, rclass
    version 17.0

    /*
        Internal bootstrap refit for rdcomono.

        This program deliberately does NOT recompute:
            - opposite-side nearest neighbors;
            - distances to opposite-side neighbors;
            - W0 and W1.
        -> Calculate nearest neighbor once! ... it doesn't depend on weights

        Those objects depend only on X, D, and the original first-stage
        bandwidths, so they are fixed across multiplier-bootstrap draws.

        For a bootstrap draw b, the caller supplies

            w_i^b = w_i * e_ib,

        where e_ib ~ Exp(1).

        generate() contains four variables:

            y0_boot
            y1_boot
            q0_grid_boot
            q1_grid_boot
    */

    syntax varlist(min=2 numeric),                         ///
        TREATment(varname numeric)                        ///
        GENerate(string asis)                             ///
        BAND0(real)                                       ///
        BAND1(real)                                       ///
        NEAREST0(varlist numeric)                         ///
        NEAREST1(varlist numeric)                         ///
        W0(varname numeric)                               ///
        W1(varname numeric)                               ///
        WVAR(varname numeric)                             ///
        Q0GRID(varname numeric)                           ///
        Q1GRID(varname numeric)                           ///
        [ KERNEL(string)                                  ///
          FOLDS(integer 5)                                ///
          ORDER(integer 1) ]

    gettoken depvar xvars : varlist

    /******************************************************************
    Parse output names
    ******************************************************************/

    local generate : list retokenize generate
    local n_generate : word count `generate'

    if `n_generate' != 4 {
        display as error ///
            "generate() must contain y0 y1 q0_grid q1_grid"
        exit 198
    }

    gettoken y0_name rest : generate
    gettoken y1_name rest : rest
    gettoken q0_name q1_name : rest

    confirm new variable `y0_name'
    confirm new variable `y1_name'
    confirm new variable `q0_name'
    confirm new variable `q1_name'

    /******************************************************************
    Basic validation
    ******************************************************************/

    marksample touse
    markout `touse' `depvar' `xvars' `treatment' `wvar'

    if `band0' <= 0 | `band1' <= 0 {
        display as error "band0() and band1() must be positive"
        exit 198
    }

    if `folds' < 1 {
        display as error "folds() must be positive"
        exit 198
    }

    if `order' < 0 {
        display as error "order() must be nonnegative"
        exit 198
    }

    if "`kernel'" == "" {
        local kernel "gaussian"
    }

    local kernel = lower("`kernel'")

    if !inlist("`kernel'", "gaussian", "uniform", "triangular") {
        display as error ///
            "kernel() must be gaussian, uniform, or triangular"
        exit 198
    }

    local p : word count `xvars'
    local p_nearest0 : word count `nearest0'
    local p_nearest1 : word count `nearest1'

    if `p_nearest0' != `p' | `p_nearest1' != `p' {
        display as error ///
            "nearest-neighbor vectors must have same dimension as X"
        exit 198
    }

    quietly count if `touse' & `treatment' == 0
    local n0 = r(N)

    quietly count if `touse' & `treatment' == 1
    local n1 = r(N)

    if `n0' == 0 | `n1' == 0 {
        display as error ///
            "bootstrap sample must contain both treatment groups"
        exit 2000
    }

    /******************************************************************
    Evaluation copies of X
    ******************************************************************/

    local evaluation_xvars0
    local evaluation_xvars1

    foreach x of local xvars {
        tempvar evaluation_x0 evaluation_x1

        quietly generate double `evaluation_x0' = `x' if `touse' & `treatment' == 0

        quietly generate double `evaluation_x1' = `x' if `touse' & `treatment' == 1

        local evaluation_xvars0 "`evaluation_xvars0' `evaluation_x0'"

        local evaluation_xvars1 "`evaluation_xvars1' `evaluation_x1'"
    }

    /******************************************************************
    Bootstrap Stage 1:
    re-estimate g0 and g1.

    IMPORTANT:
    band0 and band1 are fixed at the values selected in the original
    rdcomono fit.
    ******************************************************************/

    tempvar g0_factual g1_factual

    quietly _rdcomono_localpoly `depvar' `xvars'             ///
        if `touse' & `treatment' == 0,                       ///
        at(`evaluation_xvars0')                               ///
        generate(`g0_factual')                               ///
        bandwidth(`band0')                                   ///
        wvar(`wvar')                                         ///
        kernel(`kernel')                                     ///
        folds(`folds')                                       ///
        order(`order')

    quietly _rdcomono_localpoly `depvar' `xvars'             ///
        if `touse' & `treatment' == 1,                       ///
        at(`evaluation_xvars1')                               ///
        generate(`g1_factual')                               ///
        bandwidth(`band1')                                   ///
        wvar(`wvar')                                         ///
        kernel(`kernel')                                     ///
        folds(`folds')                                       ///
        order(`order')

    /******************************************************************
    Bootstrap Stage 2:
    generated regressors near the frontier.

    nearest0(), nearest1(), W0(), and W1() are inherited from
    the original fit.
    ******************************************************************/

    tempvar g0_boundary g1_boundary

    local boundary_xvars1
    local boundary_xvars0

    foreach x of local xvars {
        tempvar boundary_x1 boundary_x0

        quietly generate double `boundary_x1' = `x' if `touse' & `treatment' == 1 & `w1' == 1

        quietly generate double `boundary_x0' = `x' if `touse' & `treatment' == 0 & `w0' == 1

        local boundary_xvars1 "`boundary_xvars1' `boundary_x1'"

        local boundary_xvars0 "`boundary_xvars0' `boundary_x0'"
    }

    /*
        g0 evaluated for treated observations, with the kernel centered
        at their nearest untreated observations.
    */
    quietly _rdcomono_localpoly `depvar' `xvars'             ///
        if `touse' & `treatment' == 0,                       ///
        at(`evaluation_xvars1')                               ///
        center(`nearest0')                                   ///
        generate(`g0_boundary')                              ///
        bandwidth(`band0')                                   ///
        wvar(`wvar')                                         ///
        kernel(`kernel')                                     ///
        folds(`folds')                                       ///
        order(`order')

    /*
        g1 evaluated for untreated observations, with the kernel centered
        at their nearest treated observations.
    */
    quietly _rdcomono_localpoly `depvar' `xvars'             ///
        if `touse' & `treatment' == 1,                       ///
        at(`evaluation_xvars0')                               ///
        center(`nearest1')                                   ///
        generate(`g1_boundary')                              ///
        bandwidth(`band1')                                   ///
        wvar(`wvar')                                         ///
        kernel(`kernel')                                     ///
        folds(`folds')                                       ///
        order(`order')

    /******************************************************************
    Bootstrap Stage 3:
    q1 and q0.

    Following MATLAB/R, when band0 != band1 the two original selected
    first-stage bandwidths are the candidate bandwidths for the
    bootstrap q regressions.
    ******************************************************************/

    if abs(`band0' - `band1') < 1e-14 {
        local qbands "`band0'"
    }
    else {
        local qbands "`band0' `band1'"
    }

    local n_qbands : word count `qbands'

    quietly count if ///
        `touse' & `treatment' == 1 & `w1' == 1
    local n_frontier1 = r(N)

    quietly count if ///
        `touse' & `treatment' == 0 & `w0' == 1
    local n_frontier0 = r(N)

    if `n_frontier0' == 0 | `n_frontier1' == 0 {
        display as error ///
            "bootstrap draw has no usable near-frontier observations"
        exit 498
    }

    if `n_qbands' > 1 & ///
       (`folds' > `n_frontier0' | `folds' > `n_frontier1') {

        display as error ///
            "folds() exceeds bootstrap near-frontier sample size"
        exit 198
    }

    tempvar q1_at_g0 q0_at_g1
    tempname q1_band_scalar q0_band_scalar

    /*
        q1:
            E[Y(1)|X] as a function of E[Y(0)|X]
    */
    quietly _rdcomono_localpoly `depvar' `g0_boundary'       ///
        if `touse' & `treatment' == 1 & `w1' == 1,          ///
        at(`g0_factual')                                     ///
        generate(`q1_at_g0')                                 ///
        bandwidth(`qbands')                                  ///
        wvar(`wvar')                                         ///
        kernel(`kernel')                                     ///
        folds(`folds')                                       ///
        order(`order')

    scalar `q1_band_scalar' = r(bandwidth)

    if missing(scalar(`q1_band_scalar')) {
        display as error ///
            "no valid bootstrap q1 bandwidth"
        exit 498
    }

    /*
        q0:
            E[Y(0)|X] as a function of E[Y(1)|X]
    */
    quietly _rdcomono_localpoly `depvar' `g1_boundary'       ///
        if `touse' & `treatment' == 0 & `w0' == 1,          ///
        at(`g1_factual')                                     ///
        generate(`q0_at_g1')                                 ///
        bandwidth(`qbands')                                  ///
        wvar(`wvar')                                         ///
        kernel(`kernel')                                     ///
        folds(`folds')                                       ///
        order(`order')

    scalar `q0_band_scalar' = r(bandwidth)

    if missing(scalar(`q0_band_scalar')) {
        display as error ///
            "no valid bootstrap q0 bandwidth"
        exit 498
    }

    /******************************************************************
    Bootstrap q-curves on FIXED grids.

    We do not allow the x-axis itself to change across bootstrap draws.
    ******************************************************************/

    local q0_band_value = scalar(`q0_band_scalar')
    local q1_band_value = scalar(`q1_band_scalar')

    /*
        q0 grid:
        training regressor = g1_boundary
    */
    quietly _rdcomono_localpoly `depvar' `g1_boundary'       ///
        if `touse' & `treatment' == 0 & `w0' == 1,          ///
        at(`q0grid')                                         ///
        generate(`q0_name')                                  ///
        bandwidth(`q0_band_value')                           ///
        wvar(`wvar')                                         ///
        kernel(`kernel')                                     ///
        folds(`folds')                                       ///
        order(`order')

    /*
        q1 grid:
        training regressor = g0_boundary
    */
    quietly _rdcomono_localpoly `depvar' `g0_boundary'       ///
        if `touse' & `treatment' == 1 & `w1' == 1,          ///
        at(`q1grid')                                         ///
        generate(`q1_name')                                  ///
        bandwidth(`q1_band_value')                           ///
        wvar(`wvar')                                         ///
        kernel(`kernel')                                     ///
        folds(`folds')                                       ///
        order(`order')

    /******************************************************************
    Bootstrap potential-outcome estimates
    ******************************************************************/

    quietly generate double `y0_name' = .
    quietly generate double `y1_name' = .

    /*
        Untreated observations:
            y0 is factual
            y1 = q1(g0)
    */
    quietly replace `y0_name' = `g0_factual' ///
        if `touse' & `treatment' == 0

    quietly replace `y1_name' = `q1_at_g0' ///
        if `touse' & `treatment' == 0

    /*
        Treated observations:
            y1 is factual
            y0 = q0(g1)
    */
    quietly replace `y1_name' = `g1_factual' ///
        if `touse' & `treatment' == 1

    quietly replace `y0_name' = `q0_at_g1' ///
        if `touse' & `treatment' == 1

    /******************************************************************
    Return bootstrap q bandwidths
    ******************************************************************/

    return scalar q0_band = scalar(`q0_band_scalar')
    return scalar q1_band = scalar(`q1_band_scalar')
end
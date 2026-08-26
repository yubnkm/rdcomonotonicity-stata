*! version 0.2.0 11aug2026
program define rdcomono, rclass
    version 17.0

    /*
        Comonotonicity-based extrapolation in a sharp RDD.

        This is based on 
        RDD_extrapolate_CV_band.m and R/RDD_extrapolate_CV_band.R:

          1. estimate g0(x) and g1(x) separately;
          2. locate nearest observations across the treatment frontier;
          3. form generated regressors near the frontier;
          4. estimate q1(g0) and q0(g1);
          5. impute observation-level conditional mean potential outcomes;
          6. construct the support indicator S;
          7. Optionally performs multiplier bootstrap inference.
    */

    syntax varlist(min=2 numeric) [if] [in],                 ///
        TREATment(varname numeric)                           ///
        GENerate(string asis)                                ///
        BANDwidth(numlist min=1)                             ///
        [ WVAR(varname numeric)                              ///
          KERNEL(string)                                    ///
          FOLDS(integer 5)                                  ///
          ORDER(integer 1)                                  ///
          BOOTstrap(integer 0)                              ///
          BOOTFRAME(name)                                   ///
          BOOTPOINTS(integer 100) ]
    /*
        The first variable is the outcome. The remaining variables are the
        assignment variables/covariates X.
    */
    gettoken depvar xvars : varlist

    /*
        generate() must contain exactly three new variable names:

            generate(y0hat y1hat supported)
        
        y0hat       estimated E[Y(0)|X]
        y1hat       estimated E[Y(1)|X]
        supported   estimated support indicator S
        * Treatement effect estimate can be created with
            gen double tauhat = y1hat - y0hat if supported == 1
    */
    local generate : list retokenize generate
    local n_generate : word count `generate'

    if `n_generate' != 3 {
        display as error "generate() must contain exactly three names: y0 y1 support"
        exit 198
    }

    gettoken y0_name generated_rest : generate
    gettoken y1_name support_name : generated_rest

    confirm new variable `y0_name'
    confirm new variable `y1_name'
    confirm new variable `support_name'

    /*
        Define the complete-case estimation sample.
    */
    marksample touse
    markout `touse' `depvar' `xvars' `treatment'

    if "`wvar'" != "" {
        markout `touse' `wvar'
    }

    /*
        Basic option validation.
    */
    if `order' < 0 {
        display as error "order() must be a nonnegative integer"
        exit 198
    }

    if `folds' < 1 {
        display as error "folds() must be a positive integer"
        exit 198
    }

    foreach h of numlist `bandwidth' {
        if `h' <= 0 {
            display as error ///
                "all values in bandwidth() must be strictly positive"
            exit 198
        }
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

    /*
        Bootstrap validation.
    */
    if `bootstrap' < 0 {
        display as error "bootstrap() must be nonnegative"
        exit 198
    }

    if `bootstrap' > 0 & `bootpoints' < 2 {
        display as error ///
            "bootpoints() must be at least 2"
        exit 198
    }

    if `bootstrap' == 0 & "`bootframe'" != "" {
        display as error ///
            "bootframe() requires bootstrap()"
        exit 198
    }

    if `bootstrap' > 0 & "`bootframe'" == "" {
        local bootframe "rdcomono_bootstrap"
    }

    /*
        Use unit numerical weights when wvar() is omitted.
    */
    if "`wvar'" == "" {
        tempvar unit_weight
        quietly generate double `unit_weight' = 1
        local wvar "`unit_weight'"
    }

    quietly count if `touse' & `wvar' < 0

    if r(N) > 0 {
        display as error "wvar() must be nonnegative"
        exit 198
    }

    quietly summarize `wvar' if `touse', meanonly

    if r(sum) <= 0 {
        display as error "wvar() must contain at least one strictly positive weight"
        exit 198
    }

    /*
        Treatment must be binary and both treatment regions must be observed.
    */
    quietly count if `touse' & !inlist(`treatment', 0, 1)

    if r(N) > 0 {
        display as error "treatment() must equal 0 or 1 in the estimation sample"
        exit 198
    }

    quietly count if `touse'
    local n_complete = r(N)

    if `n_complete' == 0 {
        display as error "no complete observations in the estimation sample"
        exit 2000
    }

    quietly count if `touse' & `treatment' == 0
    local n0 = r(N)

    quietly count if `touse' & `treatment' == 1
    local n1 = r(N)

    if `n0' == 0 | `n1' == 0 {
        display as error ///
            "the estimation sample must contain both treated and untreated observations"
        exit 2000
    }

    
    local p : word count `xvars'
    local evaluation_xvars0
    local evaluation_xvars1

    foreach x of local xvars {
        tempvar evaluation_x0 evaluation_x1

        quietly generate double `evaluation_x0' = `x' if `touse' & `treatment' == 0
        quietly generate double `evaluation_x1' = `x' if `touse' & `treatment' == 1

        local evaluation_xvars0 "`evaluation_xvars0' `evaluation_x0'"
        local evaluation_xvars1 "`evaluation_xvars1' `evaluation_x1'"
    }

    /*
        Bandwidth convention.

        One value:
            use the same fixed bandwidth in every local-polynomial fit.

        Two values:
            use the first value for g0 and the second for g1. The two values
            are candidate bandwidths for the second-stage q regressions.

        Three or more values:
            use the full list as cross-validation candidates in every stage.

    */
    local n_bandwidths : word count `bandwidth'

    if `n_bandwidths' == 1 {
        local bands0 "`bandwidth'"
        local bands1 "`bandwidth'"
        local qbands "`bandwidth'"
    }
    else if `n_bandwidths' == 2 {
        local bands0 : word 1 of `bandwidth'
        local bands1 : word 2 of `bandwidth'
        local qbands "`bandwidth'"
    }
    else {
        local bands0 "`bandwidth'"
        local bands1 "`bandwidth'"
        local qbands "`bandwidth'"
    }

    if `n_bandwidths' > 1 & `folds' < 2 {
        display as error ///
            "folds() must be at least 2 when cross-validation is requested"
        exit 198
    }

    /*
        Stage 1: estimate factual conditional mean functions separately by D.
    */
    tempvar g0_factual g1_factual
    tempname band0_scalar band1_scalar

    quietly _rdcomono_localpoly `depvar' `xvars'                 ///
        if `touse' & `treatment' == 0,                          ///
        at(`evaluation_xvars0')                                  ///
        generate(`g0_factual')                                  ///
        bandwidth(`bands0')                                     ///
        wvar(`wvar')                                            ///
        kernel(`kernel')                                        ///
        folds(`folds')                                          ///
        order(`order')

    scalar `band0_scalar' = r(bandwidth)

    if missing(scalar(`band0_scalar')) {
        display as error ///
            "no candidate bandwidth produced a valid untreated first-stage fit"
        exit 498
    }

    quietly _rdcomono_localpoly `depvar' `xvars'                 ///
        if `touse' & `treatment' == 1,                          ///
        at(`evaluation_xvars1')                                  ///
        generate(`g1_factual')                                  ///
        bandwidth(`bands1')                                     ///
        wvar(`wvar')                                            ///
        kernel(`kernel')                                        ///
        folds(`folds')                                          ///
        order(`order')

    scalar `band1_scalar' = r(bandwidth)

    if missing(scalar(`band1_scalar')) {
        display as error ///
            "no candidate bandwidth produced a valid treated first-stage fit"
        exit 498
    }

    local band0_value = scalar(`band0_scalar')
    local band1_value = scalar(`band1_scalar')

    /*
        Stage 2a: for every observation, find its nearest observation on the
        opposite side of the treatment frontier. The Mata routine writes the
        opposite-side covariate vector and Euclidean distance into temporary
        variables.
    */
    local nearest0_vars
    local nearest1_vars

    forvalues j = 1/`p' {
        tempvar nearest0_j nearest1_j
        quietly generate double `nearest0_j' = .
        quietly generate double `nearest1_j' = .

        local nearest0_vars "`nearest0_vars' `nearest0_j'"
        local nearest1_vars "`nearest1_vars' `nearest1_j'"
    }

    tempvar distance0 distance1
    quietly generate double `distance0' = .
    quietly generate double `distance1' = .

    mata: _rdcomono_nn_op_driver(                    ///
        "`xvars'",                                              ///
        "`treatment'",                                         ///
        "`touse'",                                             ///
        "`nearest0_vars'",                                     ///
        "`nearest1_vars'",                                     ///
        "`distance0'",                                         ///
        "`distance1'"                                          ///
    )

    /*
        The frontier-neighborhood radius is omega times the relevant
        first-stage bandwidth, where omega is twice the 0.75 quantile of the
        kernel distribution, as in the MATLAB implementation.
    */
    tempname quartile_multiplier

    if "`kernel'" == "gaussian" {
        scalar `quartile_multiplier' = 2 * 0.67448
    }
    else if "`kernel'" == "uniform" {
        scalar `quartile_multiplier' = 2 * 0.5
    }
    else {
        scalar `quartile_multiplier' = 2 * (1 - 1/sqrt(2))
    }

    tempvar W1 W0
    quietly generate byte `W1' = 0 if `touse'
    quietly generate byte `W0' = 0 if `touse'

    quietly replace `W1' =                                      ///
        !missing(`distance0') &                                 ///
        `distance0' <                                           ///
        scalar(`quartile_multiplier') * scalar(`band0_scalar')  ///
        if `touse' & `treatment' == 1

    quietly replace `W0' =                                      ///
        !missing(`distance1') &                                 ///
        `distance1' <                                           ///
        scalar(`quartile_multiplier') * scalar(`band1_scalar')  ///
        if `touse' & `treatment' == 0

    quietly count if `touse' & `treatment' == 1 & `W1' == 1
    local n_frontier1 = r(N)

    quietly count if `touse' & `treatment' == 0 & `W0' == 1
    local n_frontier0 = r(N)

    if `n_frontier1' == 0 | `n_frontier0' == 0 {
        display as error ///
            "no observations were found sufficiently close to both sides of the frontier"
        exit 498
    }

    local n_qbands : word count `qbands'

    if `n_qbands' > 1 &                                        ///
       (`folds' > `n_frontier1' | `folds' > `n_frontier0') {
        display as error ///
            "folds() exceeds the number of near-frontier observations in a second-stage regression"
        exit 198
    }

    /*
        Stage 2b: evaluate the opposite-side first-stage regression at Xi while
        centering the kernel at Xi's nearest opposite-side neighbor.

        g0_boundary is constructed for treated observations near the frontier;
        g1_boundary is constructed for untreated observations near the frontier.
    */

    local boundary_xvars1
    local boundary_xvars0

    foreach x of local xvars{
        tempvar boundary_x1 boundary_xv0

        quietly generate double `boundary_x1' = `x' if `touse' & `treatment' == 1 & `W1' == 1
        quietly generate double `boundary_x0' = `x' if `touse' & `treatment' == 0 & 'W0' == 0

        local boundary_xvars1 "`boundary_xvars1' `boundary_x1'"
        local boundary_xvars0 "`boundary_xvars0' `boundary_x0'"
    }

    tempvar g0_boundary g1_boundary

    quietly _rdcomono_localpoly `depvar' `xvars'                 ///
        if `touse' & `treatment' == 0,                          ///
        at(`evaluation_xvars1')                                  ///
        center(`nearest0_vars')                                 ///
        generate(`g0_boundary')                                 ///
        bandwidth(`band0_value')                    ///
        wvar(`wvar')                                            ///
        kernel(`kernel')                                        ///
        folds(`folds')                                          ///
        order(`order')

    quietly _rdcomono_localpoly `depvar' `xvars'                 ///
        if `touse' & `treatment' == 1,                          ///
        at(`evaluation_xvars0')                                  ///
        center(`nearest1_vars')                                 ///
        generate(`g1_boundary')                                 ///
        bandwidth(`band1_value')                    ///
        wvar(`wvar')                                            ///
        kernel(`kernel')                                        ///
        folds(`folds')                                          ///
        order(`order')

    /*
        Estimated domains of q1 and q0.

        q1 maps g0 to g1, so its domain is determined by g0_boundary.
        q0 maps g1 to g0, so its domain is determined by g1_boundary.
    */
    tempname y1_min_scalar y1_max_scalar
    tempname y0_min_scalar y0_max_scalar

    quietly summarize `g0_boundary'                             ///
        if `touse' & `treatment' == 1 & `W1' == 1, meanonly

    if r(N) == 0 | missing(r(min)) | missing(r(max)) {
        display as error ///
            "the generated regressor for q1 contains no valid observations"
        exit 498
    }

    scalar `y1_min_scalar' = r(min)
    scalar `y1_max_scalar' = r(max)

    quietly summarize `g1_boundary'                             ///
        if `touse' & `treatment' == 0 & `W0' == 1, meanonly

    if r(N) == 0 | missing(r(min)) | missing(r(max)) {
        display as error ///
            "the generated regressor for q0 contains no valid observations"
        exit 498
    }

    scalar `y0_min_scalar' = r(min)
    scalar `y0_max_scalar' = r(max)

    /*
        Stage 3: estimate q1 and q0 and evaluate them at each observation's
        factual conditional mean.
    */
    tempvar q1_at_g0 q0_at_g1
    tempname q1_band_scalar q0_band_scalar

    quietly _rdcomono_localpoly `depvar' `g0_boundary'          ///
        if `touse' & `treatment' == 1 & `W1' == 1,             ///
        at(`g0_factual')                                        ///
        generate(`q1_at_g0')                                    ///
        bandwidth(`qbands')                                     ///
        wvar(`wvar')                                            ///
        kernel(`kernel')                                        ///
        folds(`folds')                                          ///
        order(`order')

    scalar `q1_band_scalar' = r(bandwidth)

    if missing(scalar(`q1_band_scalar')) {
        display as error ///
            "no candidate bandwidth produced a valid q1 fit"
        exit 498
    }

    quietly _rdcomono_localpoly `depvar' `g1_boundary'          ///
        if `touse' & `treatment' == 0 & `W0' == 1,             ///
        at(`g1_factual')                                        ///
        generate(`q0_at_g1')                                    ///
        bandwidth(`qbands')                                     ///
        wvar(`wvar')                                            ///
        kernel(`kernel')                                        ///
        folds(`folds')                                          ///
        order(`order')

    scalar `q0_band_scalar' = r(bandwidth)

    if missing(scalar(`q0_band_scalar')) {
        display as error ///
            "no candidate bandwidth produced a valid q0 fit"
        exit 498
    }

    /*
        Stage 4: combine factual and imputed conditional mean potential outcomes.

          D = 0: y0 is factual and y1 = q1(g0)
          D = 1: y1 is factual and y0 = q0(g1)
    */
    quietly generate double `y0_name' = .
    quietly generate double `y1_name' = .
    quietly generate byte   `support_name' = .

    quietly replace `y0_name' = `g0_factual'                   ///
        if `touse' & `treatment' == 0
    quietly replace `y0_name' = `q0_at_g1'                     ///
        if `touse' & `treatment' == 1

    quietly replace `y1_name' = `q1_at_g0'                     ///
        if `touse' & `treatment' == 0
    quietly replace `y1_name' = `g1_factual'                   ///
        if `touse' & `treatment' == 1

    quietly replace `support_name' = 0 if `touse'

    quietly replace `support_name' =                            ///
        !missing(`g1_factual') &                                ///
        `g1_factual' > scalar(`y0_min_scalar') &                ///
        `g1_factual' < scalar(`y0_max_scalar')                  ///
        if `touse' & `treatment' == 1

    quietly replace `support_name' =                            ///
        !missing(`g0_factual') &                                ///
        `g0_factual' > scalar(`y1_min_scalar') &                ///
        `g0_factual' < scalar(`y1_max_scalar')                  ///
        if `touse' & `treatment' == 0

    label variable `y0_name' ///
        "Estimated E[Y(0)|X] from rdcomono"
    label variable `y1_name' ///
        "Estimated E[Y(1)|X] from rdcomono"
    label variable `support_name' ///
        "Counterfactual mean supported by rdcomono"

    quietly count if `touse' & `support_name' == 1
    local n_supported = r(N)

    /******************************************************************
    Multiplier bootstrap
    ******************************************************************/

    local bootstrap_frame ""

    tempname bootstrap_q0_bands
    tempname bootstrap_q1_bands

    if `bootstrap' > 0 {

        /*
            Convert temporary scalars to numerical locals before passing
            them to the bootstrap program.
        */
        local band0_value   = scalar(`band0_scalar')
        local band1_value   = scalar(`band1_scalar')

        local q0_band_value = scalar(`q0_band_scalar')
        local q1_band_value = scalar(`q1_band_scalar')

        local y0_min_value  = scalar(`y0_min_scalar')
        local y0_max_value  = scalar(`y0_max_scalar')

        local y1_min_value  = scalar(`y1_min_scalar')
        local y1_max_value  = scalar(`y1_max_scalar')

        /*
            The nearest-neighbor variables and W0/W1 have already been
            constructed by the point estimator. They are reused in every
            bootstrap replication.
        */
        quietly _rdcomono_bootstrap `depvar' `xvars'          ///
            if `touse',                                      ///
            treatment(`treatment')                           ///
            basey0(`y0_name')                                ///
            basey1(`y1_name')                                ///
            support(`support_name')                          ///
            g0boundary(`g0_boundary')                        ///
            g1boundary(`g1_boundary')                        ///
            band0(`band0_value')                             ///
            band1(`band1_value')                             ///
            q0band(`q0_band_value')                          ///
            q1band(`q1_band_value')                          ///
            y0min(`y0_min_value')                            ///
            y0max(`y0_max_value')                            ///
            y1min(`y1_min_value')                            ///
            y1max(`y1_max_value')                            ///
            nearest0(`nearest0_vars')                        ///
            nearest1(`nearest1_vars')                        ///
            w0(`W0')                                         ///
            w1(`W1')                                         ///
            wvar(`wvar')                                     ///
            reps(`bootstrap')                                ///
            frame(`bootframe')                               ///
            points(`bootpoints')                             ///
            kernel(`kernel')                                 ///
            folds(`folds')                                   ///
            order(`order')

        local bootstrap_frame "`r(frame)'"

        matrix `bootstrap_q0_bands' = r(q0_bands)
        matrix `bootstrap_q1_bands' = r(q1_bands)
    }

    /*
        Returned results for testing and later postestimation commands.
    */
    return scalar N = `n_complete'
    return scalar N0 = `n0'
    return scalar N1 = `n1'
    return scalar N_frontier0 = `n_frontier0'
    return scalar N_frontier1 = `n_frontier1'
    return scalar N_supported = `n_supported'

    return scalar band0 = scalar(`band0_scalar')
    return scalar band1 = scalar(`band1_scalar')
    return scalar q0_band = scalar(`q0_band_scalar')
    return scalar q1_band = scalar(`q1_band_scalar')

    return scalar y0_min = scalar(`y0_min_scalar')
    return scalar y0_max = scalar(`y0_max_scalar')
    return scalar y1_min = scalar(`y1_min_scalar')
    return scalar y1_max = scalar(`y1_max_scalar')

    return scalar folds = `folds'
    return scalar order = `order'

    return local depvar "`depvar'"
    return local xvars "`xvars'"
    return local treatment "`treatment'"
    return local kernel "`kernel'"
    return local y0var "`y0_name'"
    return local y1var "`y1_name'"
    return local supportvar "`support_name'"
    return local candidate_bandwidths "`bandwidth'"

    return scalar bootstrap_reps = `bootstrap'

    if `bootstrap' > 0 {

        return scalar bootstrap_points = `bootpoints'

        return matrix bootstrap_q0_bands = ///
            `bootstrap_q0_bands'

        return matrix bootstrap_q1_bands = ///
            `bootstrap_q1_bands'

        return local bootstrap_frame ///
            "`bootstrap_frame'"
    }

    display as text _newline "Comonotonic RD extrapolation"
    display as text "  complete observations: " as result %10.0f `n_complete'
    display as text "  untreated / treated:  " as result ///
        %10.0f `n0' " / " %10.0f `n1'
    display as text "  supported observations:" as result ///
        %10.0f `n_supported'
    display as text "  first-stage bandwidths:" as result ///
        %10.6f scalar(`band0_scalar') " / " ///
        %10.6f scalar(`band1_scalar')
    if `bootstrap' > 0 {

        display as text ///
            "  bootstrap replications:" ///
            as result %10.0f `bootstrap'

        display as text ///
            "  bootstrap frame:       " ///
            as result "`bootstrap_frame'"
    }
end


mata:
mata set matastrict on


/**********************************************************************
Find the nearest observation on the opposite side of the frontier
**********************************************************************/

void _rdcomono_nn_op_driver(
    string scalar x_names,
    string scalar treatment_name,
    string scalar touse_name,
    string scalar nearest0_names,
    string scalar nearest1_names,
    string scalar distance0_name,
    string scalar distance1_name
)
{
    real colvector use
    real colvector sample_index

    real matrix X
    real colvector D

    real colvector index0
    real colvector index1

    real matrix nearest0
    real matrix nearest1
    real colvector distance0
    real colvector distance1

    real scalar i
    real scalar nearest_position

    real colvector squared_distance
    real colvector nearest_candidates

    /*
        Read only the complete estimation sample, but preserve the original
        Stata observation indices so results can be written back correctly.
    */
    use = st_data(., touse_name)
    sample_index = selectindex(use :== 1)

    X = st_data(sample_index, tokens(x_names))
    D = st_data(sample_index, treatment_name)

    index0 = selectindex(D :== 0)
    index1 = selectindex(D :== 1)

    nearest0 = J(rows(X), cols(X), .)
    nearest1 = J(rows(X), cols(X), .)
    distance0 = J(rows(X), 1, .)
    distance1 = J(rows(X), 1, .)

    /*
        For every treated observation, find its nearest untreated neighbor.
        Ties are resolved by taking the first minimum, matching the deterministic
        behavior needed for reproducible tests.
    */
    for (i = 1; i <= rows(index1); i++) {
        squared_distance =
            rowsum(
                (X[index0, .] :- X[index1[i], .]) :^ 2
            )

        nearest_candidates =
            selectindex(
                squared_distance :== min(squared_distance)
            )

        nearest_position = nearest_candidates[1]

        nearest0[index1[i], .] = X[index0[nearest_position], .]
        distance0[index1[i]] = sqrt(squared_distance[nearest_position])
    }

    /*
        For every untreated observation, find its nearest treated neighbor.
    */
    for (i = 1; i <= rows(index0); i++) {
        squared_distance =
            rowsum(
                (X[index1, .] :- X[index0[i], .]) :^ 2
            )

        nearest_candidates =
            selectindex(
                squared_distance :== min(squared_distance)
            )

        nearest_position = nearest_candidates[1]

        nearest1[index0[i], .] = X[index1[nearest_position], .]
        distance1[index0[i]] = sqrt(squared_distance[nearest_position])
    }

    st_store(
        sample_index,
        tokens(nearest0_names),
        nearest0
    )

    st_store(
        sample_index,
        tokens(nearest1_names),
        nearest1
    )

    st_store(
        sample_index,
        distance0_name,
        distance0
    )

    st_store(
        sample_index,
        distance1_name,
        distance1
    )
}

end

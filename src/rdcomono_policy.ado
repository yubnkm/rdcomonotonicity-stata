*! version 0.1.0 17aug2026
program define rdcomono_policy, rclass
    version 17.0

    /*
        Estimate the effect of a counterfactual treatment policy after rdcomono.

        The counterfactual policy is supplied as a variable containing 0 <= p_i <= 1.

        theta_hat = (1 / sum_i S_i) * sum_i { S_i * ( (1 - D_i) * p_i * (y1_i - Y_i) + D_i * (1 - p_i) * (y0_i - Y_i) ) }

        If bootframe() is supplied, CIs are calculated.
    */


    syntax varlist(min=1 max=1 numeric) [if] [in],           ///
        TREATment(varname numeric)                           ///
        POLICY(varname numeric)                             ///
        Y0(varname numeric)                                 ///
        Y1(varname numeric)                                 ///
        SUPPORT(varname numeric)                            ///
        [ BOOTFRAME(name)                                   ///
          LEVEL(real 90) ]

    local depvar `varlist'

    marksample touse

    markout `touse'                                         ///
        `depvar'                                            ///
        `treatment'                                         ///
        `policy'                                            ///
        `y0'                                                ///
        `y1'                                                ///
        `support'


    /******************************************************************
    Validate options and variables
    ******************************************************************/

    if `level' <= 0 | `level' >= 100 {
        display as error "level() must be strictly between 0 and 100"
        exit 198
    }


    /*
        Observed treatment must be binary.
    */
    quietly count if `touse' & !inlist(`treatment', 0, 1)

    if r(N) > 0 {
        display as error "treatment() must contain only 0 and 1"
        exit 198
    }


    /*
        S must be binary.
    */
    quietly count if `touse' & !inlist(`support', 0, 1)

    if r(N) > 0 {
        display as error "support() must contain only 0 and 1"
        exit 198
    }


    /*
        Counterfactual treatment probabilities must lie in [0,1].
    */
    quietly count if `touse' &                      ///
        (`policy' < 0 | `policy' > 1)

    if r(N) > 0 {
        display as error ///
            "policy() must contain treatment probabilities between 0 and 1"
        exit 198
    }


    /*
        Count complete and supported observations.
    */
    quietly count if `touse'
    local n_complete = r(N)

    if `n_complete' == 0 {
        display as error "no complete observations"
        exit 2000
    }


    quietly count if `touse' & `support' == 1
    local n_supported = r(N)

    if `n_supported' == 0 {
        display as error ///
            "no observations have support() == 1; policy effect is not estimable"
        exit 498
    }


    /******************************************************************
    Point estimate
    ******************************************************************/

    /*
        switch_prob is the probability that counterfactual treatment
        differs from factual treatment.

        If D = 0: switch_prob = p(X)
        If D = 1: switch_prob = 1 - p(X)
    */

    tempvar switch_prob contribution

    quietly generate double `switch_prob' =                 ///
        (1 - `treatment') * `policy'                         ///
        +                                                    ///
        `treatment' * (1 - `policy')                         ///
        if `touse'


    /*
        Observation-level contribution to equation (3.4).

        D = 0:
            policy may switch the observation into treatment,
            so use estimated Y(1) - factual Y.

        D = 1:
            policy may switch the observation out of treatment,
            so use estimated Y(0) - factual Y.
    */

    quietly generate double `contribution' =                 ///
        (1 - `treatment') * `policy' *                       ///
            (`y1' - `depvar')                                ///
        +                                                    ///
        `treatment' * (1 - `policy') *                       ///
            (`y0' - `depvar')                                ///
        if `touse' & `support' == 1

    tempname estimate
    tempname num_affected
    tempname affected_share

    quietly summarize `contribution'                         ///
        if `touse' & `support' == 1, meanonly

    scalar `estimate' = r(mean)


    /*
        Expected number whose treatment changes.
    */

    quietly summarize `switch_prob'                          ///
        if `touse' & `support' == 1, meanonly

    scalar `affected_share' = r(mean)

    scalar `num_affected' =                                  ///
        scalar(`affected_share') * `n_supported'


    /******************************************************************
    Initialize confidence interval results
    ******************************************************************/

    tempname conf_rad
    tempname conf_low
    tempname conf_high

    scalar `conf_rad'  = .
    scalar `conf_low'  = .
    scalar `conf_high' = .

    local bootstrap_reps 0
    local bootstrap_valid 0


    /******************************************************************
    Multiplier bootstrap
    ******************************************************************/

    if "`bootframe'" != "" {


        /*
            Check that the requested bootstrap frame exists.
        */

        capture frame `bootframe': describe

        if _rc {
            display as error ///
                "bootstrap frame `bootframe' does not exist"
            exit 111
        }


        /*
            The bootstrap frame should have been created by rdcomono.
        */

        foreach v in _rdcomono_id `depvar' `treatment' `support' {

            capture frame `bootframe': confirm variable `v'

            if _rc {
                display as error ///
                    "variable `v' not found in bootstrap frame `bootframe'"
                display as error ///
                    "bootframe() must be the frame created by the corresponding rdcomono fit"
                exit 111
            }
        }


        /**************************************************************
        Copy the counterfactual policy into the bootstrap frame

        rdcomono was estimated BEFORE policy() was usually created,
        so the bootstrap frame does not contain the policy variable.

        _rdcomono_id identifies the original observation number.
        We temporarily link the bootstrap frame back to the main frame
        and copy the policy and touse variables.
        **************************************************************/

        local mainframe "`c(frame)'"

        tempvar policy_id
        tempvar policy_copy
        tempvar sample_copy
        tempvar bootstrap_link

        quietly generate long `policy_id' = _n

        quietly generate double `policy_copy' = `policy'

        quietly generate byte `sample_copy' = `touse'


        /*
            Create the same observation-ID variable in the bootstrap frame.
        */

        frame `bootframe': quietly generate long `policy_id' = _rdcomono_id


        /*
            Link bootstrap observations to observations in the original data frame.
        */

        frame `bootframe': quietly frlink m:1 `policy_id',    ///
            frame(`mainframe')                               ///
            generate(`bootstrap_link')


        /*
            Copy the policy probabilities and estimation-sample indicator into the bootstrap frame.
        */
        frame `bootframe': quietly frget                     ///
            `policy_copy' = `policy_copy'                    ///
            `sample_copy' = `sample_copy',                   ///
            from(`bootstrap_link')


        /**************************************************************
        Find bootstrap multiplier variables

            _rdm_m1
            _rdm_m2
            ...
        **************************************************************/

        frame `bootframe': quietly ds _rdm_m*

        local multiplier_vars "`r(varlist)'"

        local bootstrap_reps : word count `multiplier_vars'


        if `bootstrap_reps' == 0 {

            frame `bootframe': quietly drop                  ///
                `policy_id'                                  ///
                `bootstrap_link'                             ///
                `policy_copy'                                ///
                `sample_copy'

            display as error ///
                "no multiplier-bootstrap draws were found in `bootframe'"

            exit 498
        }

        tempname bootstrap_estimates

        matrix `bootstrap_estimates' = ///
            J(1, `bootstrap_reps', .)


        tempvar bootstrap_contribution

        local j = 0
        local bad_bootstrap = 0


        foreach multiplier of local multiplier_vars {

            local ++j

            local b : subinstr local multiplier ///
                "_rdm_m" "", all


            local y0_draw "_rdm_y0_`b'"
            local y1_draw "_rdm_y1_`b'"

            capture frame `bootframe': ///
                confirm variable `y0_draw'

            if _rc {

                display as error ///
                    "`y0_draw' not found in bootstrap frame"

                exit 111
            }


            capture frame `bootframe': ///
                confirm variable `y1_draw'

            if _rc {

                display as error ///
                    "`y1_draw' not found in bootstrap frame"

                exit 111
            }

            frame `bootframe': quietly generate double ///
                `bootstrap_contribution' =                   ///
                                                            ///
                (1 - `treatment') * `policy_copy' *          ///
                    (`y1_draw' - `depvar')                   ///
                +                                            ///
                `treatment' * (1 - `policy_copy') *          ///
                    (`y0_draw' - `depvar')                   ///
                                                            ///
                if `sample_copy' == 1 &                      ///
                   `support' == 1 &                          ///
                   !missing(_rdcomono_id)

            frame `bootframe': quietly summarize             ///
                `bootstrap_contribution'                     ///
                [aw = `multiplier']                          ///
                if `sample_copy' == 1 &                      ///
                   `support' == 1 &                          ///
                   !missing(_rdcomono_id),                   ///
                meanonly


            if r(N) == 0 | missing(r(mean)) {

                matrix `bootstrap_estimates'[1, `j'] = .

                local ++bad_bootstrap
            }
            else {

                matrix `bootstrap_estimates'[1, `j'] = ///
                    r(mean)
            }


            frame `bootframe': quietly drop ///
                `bootstrap_contribution'
        }


        /**************************************************************
        Remove temporary link variables from bootstrap frame
        **************************************************************/

        frame `bootframe': quietly drop                      ///
            `policy_id'                                      ///
            `bootstrap_link'                                 ///
            `policy_copy'                                    ///
            `sample_copy'


        local bootstrap_valid = ///
            `bootstrap_reps' - `bad_bootstrap'


        if `bad_bootstrap' > 0 {

            display as text ///
                "warning: `bad_bootstrap' non-finite bootstrap policy estimates were ignored"
        }

        if `bootstrap_valid' > 0 {

            mata: st_numscalar(                               ///
                "`conf_rad'",                                 ///
                _rdcomono_policy_absq7(                       ///
                    "`bootstrap_estimates'",                  ///
                    "`estimate'",                             ///
                    `level'/100                               ///
                )                                            ///
            )


            scalar `conf_low' = ///
                scalar(`estimate') - scalar(`conf_rad')


            scalar `conf_high' = ///
                scalar(`estimate') + scalar(`conf_rad')
        }
    }


    /******************************************************************
    Return results
    ******************************************************************/

    return scalar estimate = scalar(`estimate')
    return scalar N = `n_complete'
    return scalar N_supported = `n_supported'
    return scalar identified_n = `n_supported'
    return scalar num_affected = scalar(`num_affected')
    return scalar affected_share = scalar(`affected_share')
    return scalar level = `level'
    return scalar conf_rad = scalar(`conf_rad')
    return scalar conf_low = scalar(`conf_low')
    return scalar conf_high = scalar(`conf_high')
    return scalar bootstrap_reps = `bootstrap_reps'
    return scalar bootstrap_valid = `bootstrap_valid'


    return local depvar "`depvar'"
    return local treatment "`treatment'"
    return local policy "`policy'"
    return local y0var "`y0'"
    return local y1var "`y1'"
    return local supportvar "`support'"
    return local bootstrap_frame "`bootframe'"

    if "`bootframe'" != "" {
        return matrix bootstrap_estimates = `bootstrap_estimates'
    }

    /******************************************************************
    Display results
    ******************************************************************/

    display as text _newline
    display as text "{hline 64}"
    display as text " Counterfactual Policy Effect"
    display as text "{hline 64}"

    display as text " Estimated policy effect" ///
        _col(39) as result %10.6f scalar(`estimate')

    if "`bootframe'" != "" & !missing(scalar(`conf_rad')) {

        display as text " `level'% bootstrap CI" ///
            _col(39) as result "[" ///
            %9.6f scalar(`conf_low') ///
            ", " ///
            %9.6f scalar(`conf_high') ///
            "]"
    }

    display as text "{hline 64}"

    display as text " Identified observations (S = 1)" ///
        _col(44) as result %8.0f `n_supported'

    display as text " Expected number affected" ///
        _col(44) as result %8.3f scalar(`num_affected')

    display as text " Affected share among identified" ///
        _col(44) as result %8.4f scalar(`affected_share')

    if "`bootframe'" != "" {
        display as text " Bootstrap replications" ///
            _col(44) as result %8.0f `bootstrap_reps'
    }

    display as text "{hline 64}"

end


mata:

mata set matastrict on


real scalar _rdcomono_policy_absq7(
    string scalar matrix_name,
    string scalar center_name,
    real scalar p
)
{
    real colvector x

    real scalar center
    real scalar n
    real scalar h
    real scalar j
    real scalar gamma

    x = st_matrix(matrix_name)'

    x = select(x, x :< .)


    if (rows(x) == 0) {
        return(.)
    }

    center = st_numscalar(center_name)

    x = abs(x :- center)
    x = sort(x, 1)

    n = rows(x)


    if (n == 1) {
        return(x[1])
    }

    h = 1 + (n - 1) * p
    j = floor(h)
    gamma = h - j


    if (j >= n) {
        return(x[n])
    }

    return(
        (1 - gamma) * x[j]
        +
        gamma * x[j + 1]
    )
}

end
*! version 0.1.0 11aug2026
program define _rdcomono_bootstrap, rclass
    version 17.0

    /*
        Multiplier-bootstrap controller for rdcomono.

        For bootstrap replication b:

            e_ib ~ Exp(1)
            w_ib = w_i * e_ib

        Geometry is copied from the original fit and held fixed:
            nearest0
            nearest1
            W0
            W1
    */

    syntax varlist(min=2 numeric) [if] [in],                 ///
        TREATment(varname numeric)                           ///
        BASEY0(varname numeric)                              ///
        BASEY1(varname numeric)                              ///
        SUPPORT(varname numeric)                             ///
        G0BOUNDARY(varname numeric)                          ///
        G1BOUNDARY(varname numeric)                          ///
        BAND0(real)                                          ///
        BAND1(real)                                          ///
        Q0BAND(real)                                         ///
        Q1BAND(real)                                         ///
        Y0MIN(real)                                          ///
        Y0MAX(real)                                          ///
        Y1MIN(real)                                          ///
        Y1MAX(real)                                          ///
        NEAREST0(varlist numeric)                            ///
        NEAREST1(varlist numeric)                            ///
        W0(varname numeric)                                  ///
        W1(varname numeric)                                  ///
        WVAR(varname numeric)                                ///
        REPS(integer)                                        ///
        FRAME(name)                                          ///
        [ POINTS(integer 100)                                ///
          KERNEL(string)                                     ///
          FOLDS(integer 5)                                   ///
          ORDER(integer 1) ]

    gettoken depvar xvars : varlist

    /******************************************************************
    Validate options
    ******************************************************************/

    if `reps' < 1 {
        display as error "reps() must be positive"
        exit 198
    }

    if `points' < 2 {
        display as error "points() must be at least 2"
        exit 198
    }

    if "`kernel'" == "" {
        local kernel "gaussian"
    }

    /*
        Do not silently overwrite an existing frame.
    */
    capture frame `frame': describe

    if _rc == 0 {
        display as error ///
            "frame `frame' already exists"
        display as error ///
            "drop it first or choose another bootframe() name"
        exit 110
    }

    /******************************************************************
    Complete original estimation sample
    ******************************************************************/

    marksample touse

    markout `touse' ///
        `depvar' `xvars' `treatment' `wvar' ///
        `basey0' `basey1' `support'

    /*
        Reserved ID used to map bootstrap results back to the original
        observation number.
    */
    capture confirm variable _rdcomono_id

    if _rc == 0 {
        display as error ///
            "_rdcomono_id is reserved internally by rdcomono"
        exit 110
    }

    quietly generate long _rdcomono_id = _n if `touse'

    /*
        Avoid duplicate variables in frame put.
    */
    local copyvars ///
    "`depvar' `xvars' `treatment' `wvar' `nearest0' `nearest1' `w0' `w1' `g0boundary' `g1boundary' `basey0' `basey1' `support'"

    local copyvars : list uniq copyvars

    frame put _rdcomono_id `copyvars' if `touse', ///
        into(`frame')

    quietly drop _rdcomono_id

    /******************************************************************
    Standardized original observation weights
    ******************************************************************/

    frame `frame': ///
        generate double _rdm_basew = `wvar'

    /******************************************************************
    Make room for the q grids.

    Usually N >> points, but this also works in tiny test datasets.
    ******************************************************************/

    frame `frame': quietly count
    local boot_N = r(N)

    if `points' > `boot_N' {
        frame `frame': quietly set obs `points'
    }

    /*
        q0's argument is E[Y(1)|X], whose estimated frontier range is
            [y0min, y0max].

        q1's argument is E[Y(0)|X], whose estimated frontier range is
            [y1min, y1max].
    */

    frame `frame': ///
        generate double _rdm_q0_grid = ///
        `y0min' +                                ///
        (_n - 1) * (`y0max' - `y0min') / (`points' - 1) ///
        if _n <= `points'

    frame `frame': ///
        generate double _rdm_q1_grid = ///
        `y1min' +                                ///
        (_n - 1) * (`y1max' - `y1min') / (`points' - 1) ///
        if _n <= `points'

    /******************************************************************
    Evaluate ORIGINAL q0 and q1 on those fixed grids
    ******************************************************************/

    /*
        q0 uses untreated observations near the frontier:
            outcome Y
            regressor g1_boundary
    */
    frame `frame': quietly ///
        _rdcomono_localpoly `depvar' `g1boundary'            ///
        if !missing(_rdcomono_id) &                          ///
           `treatment' == 0 & `w0' == 1,                    ///
        at(_rdm_q0_grid)                                     ///
        generate(_rdm_q0)                                    ///
        bandwidth(`q0band')                                  ///
        wvar(_rdm_basew)                                     ///
        kernel(`kernel')                                     ///
        folds(`folds')                                       ///
        order(`order')

    /*
        q1 uses treated observations near the frontier:
            outcome Y
            regressor g0_boundary
    */
    frame `frame': quietly ///
        _rdcomono_localpoly `depvar' `g0boundary'            ///
        if !missing(_rdcomono_id) &                          ///
           `treatment' == 1 & `w1' == 1,                    ///
        at(_rdm_q1_grid)                                     ///
        generate(_rdm_q1)                                    ///
        bandwidth(`q1band')                                  ///
        wvar(_rdm_basew)                                     ///
        kernel(`kernel')                                     ///
        folds(`folds')                                       ///
        order(`order')

    /******************************************************************
    Store bootstrap-selected q bandwidths
    ******************************************************************/

    tempname q0_boot_bands q1_boot_bands

    matrix `q0_boot_bands' = J(1, `reps', .)
    matrix `q1_boot_bands' = J(1, `reps', .)

    /******************************************************************
    Bootstrap loop
    ******************************************************************/

    forvalues b = 1/`reps' {

        /*
            Persistent variables kept in the bootstrap frame.
        */
        local multiplier_name "_rdm_m`b'"
        local y0_name         "_rdm_y0_`b'"
        local y1_name         "_rdm_y1_`b'"
        local q0_name         "_rdm_q0_`b'"
        local q1_name         "_rdm_q1_`b'"

        /*
            Temporary effective bootstrap weight.
        */
        tempvar bootstrap_weight
        local bootstrap_weight "`bootstrap_weight'"

        /*
            If U ~ Uniform(0,1), then -ln(1-U) ~ Exp(1).

            Generate only for actual sample observations. Grid-only rows,
            if any, remain missing.
        */
        frame `frame': ///
            generate double `multiplier_name' = ///
            -ln(1 - runiform()) ///
            if !missing(_rdcomono_id)

        /*
            User weights are multiplied by the bootstrap multiplier.
        */
        frame `frame': ///
            generate double `bootstrap_weight' = ///
            _rdm_basew * `multiplier_name' ///
            if !missing(_rdcomono_id)

        /*
            Re-estimate the full regression part of rdcomono, but reuse
            original nearest-neighbor geometry.
        */
        frame `frame': quietly ///
            _rdcomono_bootstrap_refit `depvar' `xvars',      ///
            treatment(`treatment')                           ///
            generate(                                        ///
                `y0_name'                                    ///
                `y1_name'                                    ///
                `q0_name'                                    ///
                `q1_name'                                    ///
            )                                                ///
            band0(`band0')                                   ///
            band1(`band1')                                   ///
            nearest0(`nearest0')                             ///
            nearest1(`nearest1')                             ///
            w0(`w0')                                         ///
            w1(`w1')                                         ///
            wvar(`bootstrap_weight')                         ///
            q0grid(_rdm_q0_grid)                             ///
            q1grid(_rdm_q1_grid)                             ///
            kernel(`kernel')                                 ///
            folds(`folds')                                   ///
            order(`order')

        matrix `q0_boot_bands'[1, `b'] = r(q0_band)
        matrix `q1_boot_bands'[1, `b'] = r(q1_band)

        /*
            The effective weight can always be reconstructed as

                _rdm_basew * _rdm_m#

            so we do not need to store another N x B set of variables.
        */
        frame `frame': quietly drop `bootstrap_weight'

        display as text ///
            "  bootstrap replication " ///
            as result `b' ///
            as text " of " ///
            as result `reps'
    }

    /******************************************************************
    Return bootstrap information
    ******************************************************************/

    return scalar reps = `reps'
    return scalar points = `points'

    return scalar band0 = `band0'
    return scalar band1 = `band1'

    return matrix q0_bands = `q0_boot_bands'
    return matrix q1_bands = `q1_boot_bands'

    return local frame "`frame'"
end
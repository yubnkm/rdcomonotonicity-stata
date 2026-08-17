*! version 0.1.0 14aug2026

program define rdcomono_qplot, rclass
    version 17.0

    /*
        Post-estimation plotting command for rdcomono.

        This command uses a bootstrap frame created by

            rdcomono ..., bootstrap(#) bootframe(name)

        and produces:

            1. q0 plot:
                   E[Y(1)|X] -> E[Y(0)|X]

            2. q1 plot:
                   E[Y(0)|X] -> E[Y(1)|X]

            3. comparison plot:
                   q1 versus q0^{-1}

    */


    syntax, FRAME(name)                                      ///
        [ LEVEL(real 90)                                     ///
          PREFIX(name) ]


    if `level' <= 0 | `level' >= 100 {

        display as error ///
            "level() must be strictly between 0 and 100"

        exit 198
    }


    /*
        Determines the names of the graphs.
        Default: rdcomono
    */

    if "`prefix'" == "" {
        local prefix "rdcomono"
    }


    /*
        Stata object names may contain at most 32 characters.

        "_compare" is the longest suffix we add here.
    */

    if strlen("`prefix'_compare") > 32 {

        display as error ///
            "prefix() is too long to construct valid graph names"

        exit 198
    }


    local q0_graph      "`prefix'_q0"
    local q1_graph      "`prefix'_q1"
    local compare_graph "`prefix'_compare"



    /******************************************************************
        Verify that the requested bootstrap frame exists
    ******************************************************************/

    capture frame `frame': describe

    if _rc != 0 {

        display as error ///
            "frame `frame' does not exist"

        exit 111
    }



    /******************************************************************
        Verify that this is an rdcomono bootstrap frame
    ******************************************************************/

    /*
        The current rdcomono bootstrap frame must contain:

            _rdm_q0_grid
                Grid on which q0 is evaluated.

            _rdm_q1_grid
                Grid on which q1 is evaluated.

            _rdm_q0
                Original estimated q0 curve.

            _rdm_q1
                Original estimated q1 curve.

        Bootstrap curves have names such as

            _rdm_q0_1
            _rdm_q0_2
            ...

            _rdm_q1_1
            _rdm_q1_2
            ...
    */

    foreach v in                                      ///
        _rdm_q0_grid                                  ///
        _rdm_q1_grid                                  ///
        _rdm_q0                                       ///
        _rdm_q1 {

        capture frame `frame': confirm numeric variable `v'

        if _rc != 0 {

            display as error ///
                "frame `frame' is not a valid rdcomono bootstrap frame"

            display as error ///
                "required variable `v' was not found"

            exit 111
        }
    }



    /******************************************************************
        Find the bootstrap q0 and q1 variables
    ******************************************************************/

    /*
        Find everything beginning with _rdm_q0_.
        Remove _rdm_q0_grid so that only bootstrap draws remain.
    */

    capture frame `frame': ds _rdm_q0_*

    if _rc != 0 {

        display as error ///
            "no q0 bootstrap variables were found in frame `frame'"

        exit 111
    }

    local q0_all "`r(varlist)'"

    local q0_exclude "_rdm_q0_grid"

    local q0_bootvars : list q0_all - q0_exclude

    local q0_reps : word count `q0_bootvars'


    if `q0_reps' == 0 {

        display as error ///
            "no q0 bootstrap draws were found in frame `frame'"

        exit 111
    }



    /*
        Find everything beginning with _rdm_q1_.
        Remove _rdm_q1_grid so that only bootstrap draws remain.
    */

    capture frame `frame': ds _rdm_q1_*

    if _rc != 0 {

        display as error ///
            "no q1 bootstrap variables were found in frame `frame'"

        exit 111
    }

    local q1_all "`r(varlist)'"

    local q1_exclude "_rdm_q1_grid"

    local q1_bootvars : list q1_all - q1_exclude

    local q1_reps : word count `q1_bootvars'


    if `q1_reps' == 0 {

        display as error ///
            "no q1 bootstrap draws were found in frame `frame'"

        exit 111
    }
    /*
        q0 and q1 should have the same number of
        bootstrap replications.
    */

    if `q0_reps' != `q1_reps' {

        display as error ///
            "q0 and q1 have different numbers of bootstrap replications"

        display as error ///
            "q0 replications: `q0_reps'; q1 replications: `q1_reps'"

        exit 498
    }

    local reps = `q0_reps'


    /******************************************************************
        Check that valid q grids exist
    ******************************************************************/

    frame `frame': quietly count if        ///
        !missing(_rdm_q0_grid) &           ///
        !missing(_rdm_q0)

    local n_q0_points = r(N)


    frame `frame': quietly count if        ///
        !missing(_rdm_q1_grid) &           ///
        !missing(_rdm_q1)

    local n_q1_points = r(N)


    if `n_q0_points' < 2 {

        display as error ///
            "q0 requires at least two valid evaluation points"

        exit 498
    }


    if `n_q1_points' < 2 {

        display as error ///
            "q1 requires at least two valid evaluation points"

        exit 498
    }



    /******************************************************************
        Construct temporary frames
    ******************************************************************/

    /*
        temporary copies: q0frame, q1frame
    */

    tempname q0frame q1frame compareframe

    frame copy `frame' `q0frame'
    frame copy `frame' `q1frame'

    tempfile q0_compare_data
    tempfile q1_compare_data



    /******************************************************************
        Construct q0 plot data
    ******************************************************************/

    frame `q0frame': keep if                ///
        !missing(_rdm_q0_grid) &            ///
        !missing(_rdm_q0)

    frame `q0frame': rename _rdm_q0_grid q0_input
    frame `q0frame': rename _rdm_q0      q0_estimate

    frame `q0frame': keep                  ///
        q0_input                           ///
        q0_estimate                        ///
        `q0_bootvars'

    frame `q0frame': generate long q0_point = _n

    frame `q0frame': reshape long _rdm_q0_, ///
        i(q0_point)                          ///
        j(rep)



    /******************************************************************
        Calculate q0 pointwise bootstrap band
    ******************************************************************/

    frame `q0frame': generate double q0_absdev = ///
        abs(_rdm_q0_ - q0_estimate)

    frame `q0frame': bysort q0_point: egen double q0_conf = ///
        pctile(q0_absdev), p(`level')

    frame `q0frame': bysort q0_point: keep if _n == 1

    /*
        Pointwise lower and upper confidence-band endpoints.
    */

    frame `q0frame': generate double q0_lower = ///
        q0_estimate - q0_conf

    frame `q0frame': generate double q0_upper = ///
        q0_estimate + q0_conf


    /*
        45-degree line.
    */

    frame `q0frame': generate double q0_identity = ///
        q0_input


    /*
        Sort by horizontal-axis variable before graphing.
    */

    frame `q0frame': sort q0_input



    /******************************************************************
        Obtain q0 support endpoints
    ******************************************************************/

    frame `q0frame': quietly summarize q0_input, meanonly

    local q0_x_min = r(min)
    local q0_x_max = r(max)


    if `q0_x_min' >= `q0_x_max' {

        display as error ///
            "q0 evaluation grid has no nondegenerate support"

        capture frame drop `q0frame'
        capture frame drop `q1frame'

        exit 498
    }



    /******************************************************************
        Plot q0
    ******************************************************************/

    /*
        Plot order:

            1. bootstrap confidence region
            2. estimated q0
            3. 45-degree line
    */

    frame `q0frame': twoway                              ///
        (rarea q0_lower q0_upper q0_input,              ///
            sort                                        ///
            fcolor(gs10)                                ///
            fintensity(100)                             ///
            lcolor(gs10)                                ///
            lwidth(vthin))                              ///
        (line q0_estimate q0_input,                     ///
            sort                                        ///
            lcolor(black)                               ///
            lwidth(medthick))                           ///
        (line q0_identity q0_input,                     ///
            sort                                        ///
            lcolor(gs6)                                 ///
            lpattern(dash)                              ///
            lwidth(thin)),                              ///
        xline(`q0_x_min' `q0_x_max',                    ///
            lcolor(gs6)                                 ///
            lpattern(dot)                               ///
            lwidth(thin))                               ///
        title("Estimated q0 function")                   ///
        subtitle("`level'% pointwise multiplier-bootstrap band") ///
        xtitle("E[Y(1)|X = x]")                         ///
        ytitle("E[Y(0)|X = x]")                         ///
        legend(order(                                   ///
            1 "`level'% bootstrap band"                 ///
            2 "Estimated q0"                            ///
            3 "45-degree line"))                        ///
        name(`q0_graph', replace)



    /******************************************************************
        Save q0 data needed for q0^{-1}
    ******************************************************************/

    frame `q0frame': generate double EY0 = q0_estimate
    frame `q0frame': generate double EY1 = q0_input

    frame `q0frame': generate byte function_id = 2

    frame `q0frame': keep EY0 EY1 function_id

    frame `q0frame': save "`q0_compare_data'", replace



    /******************************************************************
        Construct q1 plot data
    ******************************************************************/

    frame `q1frame': keep if                ///
        !missing(_rdm_q1_grid) &            ///
        !missing(_rdm_q1)

    frame `q1frame': rename _rdm_q1_grid q1_input
    frame `q1frame': rename _rdm_q1      q1_estimate

    frame `q1frame': keep                  ///
        q1_input                           ///
        q1_estimate                        ///
        `q1_bootvars'

    frame `q1frame': generate long q1_point = _n

    frame `q1frame': reshape long _rdm_q1_, ///
        i(q1_point)                          ///
        j(rep)



    /******************************************************************
        Calculate q1 pointwise bootstrap band
    ******************************************************************/

    frame `q1frame': generate double q1_absdev = ///
        abs(_rdm_q1_ - q1_estimate)


    frame `q1frame': bysort q1_point: egen double q1_conf = ///
        pctile(q1_absdev), p(`level')


    frame `q1frame': bysort q1_point: keep if _n == 1


    frame `q1frame': generate double q1_lower = ///
        q1_estimate - q1_conf

    frame `q1frame': generate double q1_upper = ///
        q1_estimate + q1_conf


    frame `q1frame': generate double q1_identity = ///
        q1_input


    frame `q1frame': sort q1_input



    /******************************************************************
        Obtain q1 support endpoints
    ******************************************************************/

    frame `q1frame': quietly summarize q1_input, meanonly

    local q1_x_min = r(min)
    local q1_x_max = r(max)


    if `q1_x_min' >= `q1_x_max' {

        display as error ///
            "q1 evaluation grid has no nondegenerate support"

        capture frame drop `q0frame'
        capture frame drop `q1frame'

        exit 498
    }



    /******************************************************************
        Plot q1
    ******************************************************************/

    frame `q1frame': twoway                              ///
        (rarea q1_lower q1_upper q1_input,              ///
            sort                                        ///
            fcolor(gs10)                                ///
            fintensity(100)                             ///
            lcolor(gs10)                                ///
            lwidth(vthin))                              ///
        (line q1_estimate q1_input,                     ///
            sort                                        ///
            lcolor(black)                               ///
            lwidth(medthick))                           ///
        (line q1_identity q1_input,                     ///
            sort                                        ///
            lcolor(gs6)                                 ///
            lpattern(dash)                              ///
            lwidth(thin)),                              ///
        xline(`q1_x_min' `q1_x_max',                    ///
            lcolor(gs6)                                 ///
            lpattern(dot)                               ///
            lwidth(thin))                               ///
        title("Estimated q1 function")                   ///
        subtitle("`level'% pointwise multiplier-bootstrap band") ///
        xtitle("E[Y(0)|X = x]")                         ///
        ytitle("E[Y(1)|X = x]")                         ///
        legend(order(                                   ///
            1 "`level'% bootstrap band"                 ///
            2 "Estimated q1"                            ///
            3 "45-degree line"))                        ///
        name(`q1_graph', replace)



    /******************************************************************
        Save q1 data for comparison graph
    ******************************************************************/

    frame `q1frame': generate double EY0 = q1_input
    frame `q1frame': generate double EY1 = q1_estimate

    frame `q1frame': generate byte function_id = 1

    frame `q1frame': keep EY0 EY1 function_id

    frame `q1frame': save "`q1_compare_data'", replace



    /******************************************************************
        Construct comparison dataset
    ******************************************************************/

    frame create `compareframe'

    frame `compareframe': use ///
        "`q1_compare_data'", clear

    frame `compareframe': append using ///
        "`q0_compare_data'"



    /******************************************************************
        Find common support for q1 and q0^{-1}
    ******************************************************************/

    frame `compareframe': quietly summarize EY0 ///
        if function_id == 1, meanonly

    local q1_comp_min = r(min)
    local q1_comp_max = r(max)

    frame `compareframe': quietly summarize EY0 ///
        if function_id == 2, meanonly

    local q0inv_comp_min = r(min)
    local q0inv_comp_max = r(max)

    local comp_x_min = max(                ///
        `q1_comp_min',                     ///
        `q0inv_comp_min'                   ///
    )

    local comp_x_max = min(                ///
        `q1_comp_max',                     ///
        `q0inv_comp_max'                   ///
    )


    if `comp_x_min' >= `comp_x_max' {

        display as error ///
            "q1 and q0 inverse have no overlapping horizontal support"

        capture frame drop `q0frame'
        capture frame drop `q1frame'
        capture frame drop `compareframe'

        exit 498
    }



    /******************************************************************
        Plot q1 versus q0^{-1}
    ******************************************************************/

    frame `compareframe': twoway                         ///
        (line EY1 EY0                                    ///
            if function_id == 1 &                        ///
            inrange(EY0, `comp_x_min', `comp_x_max'),    ///
            sort                                         ///
            lcolor(blue)                                 ///
            lwidth(medthick))                            ///
        (line EY1 EY0                                    ///
            if function_id == 2 &                        ///
            inrange(EY0, `comp_x_min', `comp_x_max'),    ///
            sort                                         ///
            lcolor(red)                                  ///
            lwidth(medthick))                            ///
        (function y = x,                                 ///
            range(`comp_x_min' `comp_x_max')             ///
            lcolor(gs6)                                  ///
            lpattern(dash)                               ///
            lwidth(thin)),                               ///
        title("Comparison of q1 and q0^{-1}")          ///
        xtitle("E[Y(0)|X = x]")                          ///
        ytitle("E[Y(1)|X = x]")                          ///
        legend(order(                                    ///
            1 "Estimated q1"                             ///
            2 "Estimated q0 inverse"))                   ///
        name(`compare_graph', replace)


    /******************************************************************
        Remove temporary frames
    ******************************************************************/

    capture frame drop `q0frame'
    capture frame drop `q1frame'
    capture frame drop `compareframe'


    /******************************************************************
        Return results
    ******************************************************************/

    return scalar level = `level'
    return scalar bootstrap_reps = `reps'

    return scalar q0_points = `n_q0_points'
    return scalar q1_points = `n_q1_points'

    return scalar comparison_xmin = `comp_x_min'
    return scalar comparison_xmax = `comp_x_max'

    return local bootstrap_frame "`frame'"

    return local q0_graph "`q0_graph'"

    return local q1_graph "`q1_graph'"

    return local comparison_graph "`compare_graph'"



    /******************************************************************
        Display summary
    ******************************************************************/

    display as text _newline  "rdcomono q-function plots"

    display as text ///
        "  bootstrap frame:      " ///
        as result "`frame'"

    display as text ///
        "  bootstrap replications:" ///
        as result %9.0f `reps'

    display as text ///
        "  confidence level:     " ///
        as result %9.1f `level' "%"

    display as text ///
        "  q0 graph:             " ///
        as result "`q0_graph'"

    display as text ///
        "  q1 graph:             " ///
        as result "`q1_graph'"

    display as text ///
        "  comparison graph:     " ///
        as result "`compare_graph'"

end
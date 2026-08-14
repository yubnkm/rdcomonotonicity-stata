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

        The bootstrap confidence bands are pointwise bands:

            c_d(y)
              = level-th percentile of
                |q_d^b(y) - q_d(y)|

        and the plotted interval is

            q_d(y) +/- c_d(y).

        IMPORTANT:
        This first version requires bootstrap results.  It therefore
        requires frame(), where frame() is the bootstrap frame returned
        by rdcomono.

        showpoints is not implemented yet because the second-stage
        training data are not currently stored with stable public names
        in the bootstrap frame.
    */


    /******************************************************************
    1. Parse syntax
    ******************************************************************/

    syntax, FRAME(name)                                      ///
        [ LEVEL(real 90)                                     ///
          PREFIX(name)                                       ///
          NOCOMBINE ]


    /*
        Stata convention:

            level(90)

        means a 90 percent band.

        This differs from the R function, where level = 0.90.
    */

    if `level' <= 0 | `level' >= 100 {

        display as error ///
            "level() must be strictly between 0 and 100"

        exit 198
    }


    /*
        prefix() determines the names of the graphs.

        Default:

            rdcomono_q0
            rdcomono_q1
            rdcomono_compare
            rdcomono_qplots
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
    local combined_graph "`prefix'_qplots"



    /******************************************************************
    2. Verify that the requested bootstrap frame exists
    ******************************************************************/

    capture frame `frame': describe

    if _rc != 0 {

        display as error ///
            "frame `frame' does not exist"

        exit 111
    }



    /******************************************************************
    3. Verify that this is an rdcomono bootstrap frame
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
    4. Find the bootstrap q0 and q1 variables
    ******************************************************************/

    /*
        First find everything beginning with _rdm_q0_.

        This returns

            _rdm_q0_grid
            _rdm_q0_1
            _rdm_q0_2
            ...

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
        Do the same thing for q1.
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
        rdcomono should generate one q0 curve and one q1 curve for
        every bootstrap replication.

        If these numbers differ, something is wrong with the bootstrap
        frame.
    */

    if `q0_reps' != `q1_reps' {

        display as error ///
            "q0 and q1 have different numbers of bootstrap draws"

        exit 498
    }


    local reps = `q0_reps'



    /******************************************************************
    5. Check that valid q grids exist
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
    6. Construct temporary frames
    ******************************************************************/

    /*
        We never modify the user's bootstrap frame.

        Instead, make temporary copies:

            q0frame
            q1frame

        These copies will be reshaped and transformed for plotting.
    */

    tempname q0frame q1frame compareframe

    frame copy `frame' `q0frame'
    frame copy `frame' `q1frame'


    /*
        Temporary files are used only to build the comparison plot.

        They disappear automatically when this command finishes.
    */

    tempfile q0_compare_data
    tempfile q1_compare_data



    /******************************************************************
    7. Construct q0 plot data
    ******************************************************************/

    /*
        q0 maps

            E[Y(1)|X]  ->  E[Y(0)|X].

        Therefore:

            horizontal axis = q0 input
            vertical axis   = q0 estimate
    */


    /*
        Keep only valid grid rows.
    */

    frame `q0frame': keep if                ///
        !missing(_rdm_q0_grid) &            ///
        !missing(_rdm_q0)


    /*
        Give the original grid and estimator readable names.
    */

    frame `q0frame': rename _rdm_q0_grid q0_input
    frame `q0frame': rename _rdm_q0      q0_estimate


    /*
        Retain only the variables needed for q0 inference.
    */

    frame `q0frame': keep                  ///
        q0_input                           ///
        q0_estimate                        ///
        `q0_bootvars'


    /*
        Each q0 grid location needs an identifier before reshape.
    */

    frame `q0frame': generate long q0_point = _n



    /******************************************************************
    8. Reshape q0 bootstrap curves
    ******************************************************************/

    /*
        Before reshape, we have

            q0_point   q0_input  q0_estimate   _rdm_q0_1 ...
                1          y1       q0(y1)        q0^1(y1)
                2          y2       q0(y2)        q0^1(y2)
                ...

        reshape long turns this into

            q0_point   rep   q0_input   q0_estimate   _rdm_q0_
                1       1       y1         q0(y1)       q0^1(y1)
                1       2       y1         q0(y1)       q0^2(y1)
                ...
    */

    frame `q0frame': reshape long _rdm_q0_, ///
        i(q0_point)                          ///
        j(rep)



    /******************************************************************
    9. Calculate q0 pointwise bootstrap band
    ******************************************************************/

    /*
        For bootstrap draw b calculate

            |q0_b(y) - q0_hat(y)|.
    */

    frame `q0frame': generate double q0_absdev = ///
        abs(_rdm_q0_ - q0_estimate)


    /*
        At each grid point find the requested percentile of the
        absolute bootstrap deviations.

        For level(90), for example,

            q0_conf(y)

        is the 90th percentile of the bootstrap deviations.
    */

    frame `q0frame': bysort q0_point: egen double q0_conf = ///
        pctile(q0_absdev), p(`level')


    /*
        After computing the percentile we no longer need one row per
        bootstrap replication.

        Keep one row for each q0 evaluation point.
    */

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
    10. Obtain q0 support endpoints
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
    11. Plot q0
    ******************************************************************/

    /*
        Plot order:

            1. bootstrap confidence region
            2. estimated q0
            3. 45-degree line

        gs10 is used for the bootstrap region.
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
    12. Save q0 data needed for q0^{-1}
    ******************************************************************/

    /*
        We do NOT numerically invert q0.

        If

            y0 = q0(y1),

        then the inverse relationship can be plotted parametrically as

            horizontal coordinate = q0(y1)
            vertical coordinate   = y1.

        Therefore:

            EY0 = q0_estimate
            EY1 = q0_input

        gives the q0^{-1} curve in the comparison graph.
    */

    frame `q0frame': generate double EY0 = q0_estimate
    frame `q0frame': generate double EY1 = q0_input

    frame `q0frame': generate byte function_id = 2

    frame `q0frame': keep EY0 EY1 function_id

    frame `q0frame': save "`q0_compare_data'", replace



    /******************************************************************
    13. Construct q1 plot data
    ******************************************************************/

    /*
        q1 maps

            E[Y(0)|X]  ->  E[Y(1)|X].
    */

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



    /******************************************************************
    14. Reshape q1 bootstrap curves
    ******************************************************************/

    frame `q1frame': reshape long _rdm_q1_, ///
        i(q1_point)                          ///
        j(rep)



    /******************************************************************
    15. Calculate q1 pointwise bootstrap band
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
    16. Obtain q1 support endpoints
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
    17. Plot q1
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
    18. Save q1 data for comparison graph
    ******************************************************************/

    /*
        q1 already has the orientation needed for the comparison:

            horizontal = E[Y(0)|X]
            vertical   = E[Y(1)|X].
    */

    frame `q1frame': generate double EY0 = q1_input
    frame `q1frame': generate double EY1 = q1_estimate

    frame `q1frame': generate byte function_id = 1

    frame `q1frame': keep EY0 EY1 function_id

    frame `q1frame': save "`q1_compare_data'", replace



    /******************************************************************
    19. Construct comparison dataset
    ******************************************************************/

    /*
        function_id = 1:
            q1

        function_id = 2:
            q0^{-1}
    */

    frame create `compareframe'

    frame `compareframe': use ///
        "`q1_compare_data'", clear

    frame `compareframe': append using ///
        "`q0_compare_data'"



    /******************************************************************
    20. Find common support for q1 and q0^{-1}
    ******************************************************************/

    /*
        This reproduces the R calculation

            comp_x_min =
                max(
                    min(q1 x),
                    min(q0 inverse x)
                )

            comp_x_max =
                min(
                    max(q1 x),
                    max(q0 inverse x)
                )
    */


    /*
        q1 support.
    */

    frame `compareframe': quietly summarize EY0 ///
        if function_id == 1, meanonly

    local q1_comp_min = r(min)
    local q1_comp_max = r(max)


    /*
        q0 inverse support.
    */

    frame `compareframe': quietly summarize EY0 ///
        if function_id == 2, meanonly

    local q0inv_comp_min = r(min)
    local q0inv_comp_max = r(max)


    /*
        Intersection.
    */

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
    21. Plot q1 versus q0^{-1}
    ******************************************************************/

    /*
        Under comonotonicity,

            q1(y)

        and

            q0^{-1}(y)

        should coincide over their common support.

        This is the Stata counterpart of comp_plot in the R package.
    */

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
        title("Comparison of q1 and q0 inverse")          ///
        xtitle("E[Y(0)|X = x]")                          ///
        ytitle("E[Y(1)|X = x]")                          ///
        legend(order(                                    ///
            1 "Estimated q1"                             ///
            2 "Estimated q0 inverse"))                   ///
        name(`compare_graph', replace)



    /******************************************************************
    22. Optionally combine q0 and q1 plots
    ******************************************************************/

    /*
        By default, also construct one graph containing the q0 and q1
        panels.

        The user can suppress this with

            nocombine
    */

    if "`nocombine'" == "" {

        graph combine                              ///
            `q0_graph'                             ///
            `q1_graph',                            ///
            cols(2)                                ///
            title("Estimated comonotonic mappings") ///
            name(`combined_graph', replace)
    }



    /******************************************************************
    23. Remove temporary frames
    ******************************************************************/

    /*
        The source bootstrap frame is NOT modified or dropped.
    */

    capture frame drop `q0frame'
    capture frame drop `q1frame'
    capture frame drop `compareframe'



    /******************************************************************
    24. Return results
    ******************************************************************/

    return scalar level = `level'
    return scalar bootstrap_reps = `reps'

    return scalar q0_points = `n_q0_points'
    return scalar q1_points = `n_q1_points'

    return scalar comparison_xmin = `comp_x_min'
    return scalar comparison_xmax = `comp_x_max'

    return local bootstrap_frame "`frame'"

    return local q0_graph ///
        "`q0_graph'"

    return local q1_graph ///
        "`q1_graph'"

    return local comparison_graph ///
        "`compare_graph'"


    if "`nocombine'" == "" {

        return local combined_graph ///
            "`combined_graph'"
    }



    /******************************************************************
    25. Display summary
    ******************************************************************/

    display as text _newline ///
        "rdcomono q-function plots"

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


    if "`nocombine'" == "" {

        display as text ///
            "  combined graph:       " ///
            as result "`combined_graph'"
    }

end
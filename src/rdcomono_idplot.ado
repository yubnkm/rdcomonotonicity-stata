*! version 0.1.0 14aug2026
program define rdcomono_idplot, rclass
    version 17.0

    /*
        Plot the region with identified conditional average treatment effects.

        This is the Stata counterpart of the R function plot_idf_region().

        The command is available only for two-dimensional assignment variables.

        Treatment status is represented by marker color.

        Identification status is represented by marker shape:

            support = 0   open circle
            support = 1   x

        Unlike the R function, this command does not require a
        threshold_function(). rdcomono already uses an observed treatment
        indicator, so treatment() supplies the treated/untreated region labels
        directly.
    */

    syntax varlist(min=2 max=2 numeric) [if] [in],          ///
        TREATment(varname numeric)                           ///
        SUPPORT(varname numeric)                             ///
        [ NAME(name)                                        ///
          TITLE(string asis) ]


    /******************************************************************
    1. Parse the two assignment variables
    ******************************************************************/

    gettoken x1 x2 : varlist


    /******************************************************************
    2. Define the plotting sample
    ******************************************************************/

    marksample touse

    markout `touse' ///
        `x1' ///
        `x2' ///
        `treatment' ///
        `support'


    quietly count if `touse'

    local N = r(N)


    if `N' == 0 {

        display as error ///
            "no complete observations are available for plotting"

        exit 2000
    }


    /******************************************************************
    3. Validate treatment and support indicators
    ******************************************************************/

    /*
        Treatment must be binary.
    */

    quietly count if ///
        `touse' & ///
        !inlist(`treatment', 0, 1)


    if r(N) > 0 {

        display as error ///
            "treatment() must equal 0 or 1 in the plotting sample"

        exit 198
    }


    /*
        The rdcomono support indicator must also be binary.
    */

    quietly count if ///
        `touse' & ///
        !inlist(`support', 0, 1)


    if r(N) > 0 {

        display as error ///
            "support() must equal 0 or 1 in the plotting sample"

        exit 198
    }


    /******************************************************************
    4. Graph labels and graph name
    ******************************************************************/

    /*
        Use Stata variable labels for the two axes when available.

        Otherwise use the variable names themselves.
    */

    local x1_label : variable label `x1'
    local x2_label : variable label `x2'


    if `"`x1_label'"' == "" {
        local x1_label "`x1'"
    }


    if `"`x2_label'"' == "" {
        local x2_label "`x2'"
    }


    /*
        Default graph name.
    */

    if "`name'" == "" {
        local name "rdcomono_idplot"
    }


    /*
        Default graph title.
    */

    if `"`title'"' == "" {

        local title ///
            `"Estimated extrapolation region"'
    }


    /******************************************************************
    5. Count observations by identification and treatment status
    ******************************************************************/

    quietly count if ///
        `touse' & ///
        `support' == 1

    local N_identified = r(N)


    quietly count if ///
        `touse' & ///
        `support' == 0

    local N_notidentified = r(N)


    quietly count if ///
        `touse' & ///
        `treatment' == 0

    local N0 = r(N)


    quietly count if ///
        `touse' & ///
        `treatment' == 1

    local N1 = r(N)


    /******************************************************************
    6. Plot identified and unidentified observations
    ******************************************************************/

    /*
        The R implementation maps:

            color -> treatment region
            shape -> CATE identification

        Stata legends are plot-based rather than aesthetic-based.

        Therefore, the four treatment-by-identification combinations
        are plotted separately:

            untreated, not identified
            untreated, identified
            treated, not identified
            treated, identified

        Marker conventions follow the R function:

            Not identified: open circle
            Identified:     x
    */

    twoway                                                        ///
    (scatter `x2' `x1' if                                   ///
        `touse' &                                           ///
        `treatment' == 0 &                                  ///
        `support' == 0,                                     ///
        msymbol(Oh)                                         ///
        msize(vsmall)                                       ///
        mcolor(blue))                                       ///
    (scatter `x2' `x1' if                                   ///
        `touse' &                                           ///
        `treatment' == 0 &                                  ///
        `support' == 1,                                     ///
        msymbol(x)                                          ///
        msize(vsmall)                                       ///
        mcolor(blue))                                       ///
    (scatter `x2' `x1' if                                   ///
        `touse' &                                           ///
        `treatment' == 1 &                                  ///
        `support' == 0,                                     ///
        msymbol(Oh)                                         ///
        msize(vsmall)                                       ///
        mcolor(red))                                        ///
    (scatter `x2' `x1' if                                   ///
        `touse' &                                           ///
        `treatment' == 1 &                                  ///
        `support' == 1,                                     ///
        msymbol(x)                                          ///
        msize(vsmall)                                       ///
        mcolor(red)),                                       ///
    title(`title', size(medsmall))                           ///
    xtitle(`"`x1_label'"', size(small))                      ///
    ytitle(`"`x2_label'"', size(small))                      ///
    xlabel(, labsize(small))                                ///
    ylabel(, labsize(small))                                ///
legend(order(                                               ///
    1 "Not identified"                                      ///
    2 "Identified")                                         ///
    cols(2)                                                 ///
    position(6)                                             ///
    ring(1)                                                 ///
    size(small)                                             ///
    region(lstyle(none)))                                    ///
    graphregion(color(white))                               ///
    name(`name', replace)


    /******************************************************************
    7. Return plotting information
    ******************************************************************/

    return scalar N = `N'

    return scalar N0 = `N0'
    return scalar N1 = `N1'

    return scalar N_identified = ///
        `N_identified'

    return scalar N_notidentified = ///
        `N_notidentified'


    return local xvars ///
        "`x1' `x2'"

    return local treatment ///
        "`treatment'"

    return local support ///
        "`support'"

    return local graph ///
        "`name'"


    /******************************************************************
    8. Display summary
    ******************************************************************/

    display as text _newline ///
        "rdcomono identified-region plot"

    display as text ///
        "  observations:          " ///
        as result %9.0f `N'

    display as text ///
        "  identified:            " ///
        as result %9.0f `N_identified'

    display as text ///
        "  not identified:        " ///
        as result %9.0f `N_notidentified'

    display as text ///
        "  untreated / treated:   " ///
        as result %9.0f `N0' " / " %9.0f `N1'

    display as text ///
        "  graph:                 " ///
        as result "`name'"

end
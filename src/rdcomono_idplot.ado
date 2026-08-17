*! version 0.1.0 14aug2026
program define rdcomono_idplot, rclass
    version 17.0

    /*
        Plot the region with identified conditional average treatment effects.

        The command is available only for two-dimensional assignment variables.
    */

    syntax varlist(min=2 max=2 numeric) [if] [in],          ///
        TREATment(varname numeric)                           ///
        SUPPORT(varname numeric)                             ///
        [ NAME(name)                                        ///
          TITLE(string asis) ]


    gettoken x1 x2 : varlist

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
        Graph labels and graph name
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
        Default 
    */

    if "`name'" == "" {
        local name "rdcomono_idplot"
    }

    if `"`title'"' == "" {

        local title ///
            `"Estimated extrapolation region"'
    }


    /******************************************************************
        Count observations by identification and treatment status
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
        Plot identified and unidentified observations
    ******************************************************************/

    twoway                                                        ///
    (scatter `x2' `x1' if                                   ///
        `touse' &                                           ///
        `treatment' == 0 &                                  ///
        `support' == 0,                                     ///
        msymbol(x)                                         ///
        msize(medsmall)                                       ///
        mcolor(blue))                                       ///
    (scatter `x2' `x1' if                                   ///
        `touse' &                                           ///
        `treatment' == 0 &                                  ///
        `support' == 1,                                     ///
        msymbol(Oh)                                          ///
        msize(small)                                       ///
        mcolor(blue))                                       ///
    (scatter `x2' `x1' if                                   ///
        `touse' &                                           ///
        `treatment' == 1 &                                  ///
        `support' == 0,                                     ///
        msymbol(x)                                         ///
        msize(medsmall)                                       ///
        mcolor(red))                                        ///
    (scatter `x2' `x1' if                                   ///
        `touse' &                                           ///
        `treatment' == 1 &                                  ///
        `support' == 1,                                     ///
        msymbol(Oh)                                          ///
        msize(small)                                       ///
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
        Return plotting information
    ******************************************************************/

    return scalar N = `N'
    return scalar N0 = `N0'
    return scalar N1 = `N1'

    return scalar N_identified = `N_identified'
    return scalar N_notidentified = `N_notidentified'

    return local xvars "`x1' `x2'"
    return local treatment "`treatment'"
    return local support "`support'"
    return local graph "`name'"


    /******************************************************************
        Display summary
    ******************************************************************/

    display as text _newline "rdcomono identified-region plot"

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
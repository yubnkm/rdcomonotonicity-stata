version 17.0
clear all
set more off
set linesize 120

adopath ++ "`c(pwd)'/src"

capture mkdir "output"
capture mkdir "output/summer_school"

local datafile "examples/ssextract_karthik.csv"

local B      = 10
local points = 100
local level  = 90
local bands = 0.2

set seed 1

import delimited using "`datafile'", ///
    clear varnames(1)

* Imbens and Wager only use this subsample:
keep if inrange(mdcut01, -40, 40) & inrange(rdcut01, -40, 40)

quietly count

display as text _newline "Summer-school restricted sample"
display as text "  observations after +/-40 restriction: " as result %10.0fc r(N)

* Treatment
generate byte D = !(mdcut01 > 0 & rdcut01 > 0)
label variable D "Mandatory summer-school treatment"
label define treatment_label ///
    0 "Untreated" ///
    1 "Treated"
label values D treatment_label

* Outcome
generate double Y_math = zmscr02
generate double Y_read = zrscr02

label variable Y_math "Next-year math score"
label variable Y_read "Next-year reading score"


* Normalize the covariates

quietly summarize mdcut01, meanonly

scalar min_math = r(min)
scalar max_math = r(max)

generate double x_math = (mdcut01 - scalar(min_math)) / (scalar(max_math) - scalar(min_math))
scalar cutoff_math = -scalar(min_math) / (scalar(max_math) - scalar(min_math))

quietly summarize rdcut01, meanonly

scalar min_read = r(min)
scalar max_read = r(max)

generate double x_read = (rdcut01 - scalar(min_read)) / (scalar(max_read) - scalar(min_read))

scalar cutoff_read = -scalar(min_read) / (scalar(max_read) - scalar(min_read))

label variable x_math "Math Score"
label variable x_read "Reading Score"


/*
Put the scalar values into locals so they are easy to use in graphsand policy definitions.
*/

local cm = scalar(cutoff_math)
local cr = scalar(cutoff_read)



/**********************************************************************
    Plot the empirical assignment-variable distribution
**********************************************************************/


twoway                                                     ///
    (scatter x_read x_math if D == 0,                      ///
        msymbol(Oh)                                        ///
        msize(tiny)                                        ///
        mcolor(magenta))                                   ///
    (scatter x_read x_math if D == 1,                      ///
        msymbol(+)                                         ///
        msize(tiny)                                        ///
        mcolor(cyan))                                      ///
    (pci 1 `cm' `cr' `cm',                                 ///
        lcolor(blue)                                       ///
        lwidth(medthick))                                  ///
    (pci `cr' `cm' `cr' 1,                                 ///
        lcolor(blue)                                       ///
        lwidth(medthick)),                                 ///
    title("Summer-school assignment")                      ///
    xtitle("Math Score")                                   ///
    ytitle("Reading Score")                                ///
    xscale(range(0 1))                                     ///
    yscale(range(0 1))                                     ///
    xlabel(0 1)                                            ///
    ylabel(0 1)                                            ///
    legend(order(                                          ///
        1 "Untreated"                                      ///
        2 "Treated"                                        ///
        3 "Treatment frontier"))                           ///
    name(summer_design, replace)

graph export ///
    "output/summer_school/figSummerSchool_dist.png", ///
    name(summer_design) ///
    replace width(1800)


/**********************************************************************
    Estimation
**********************************************************************/

** Math

capture frame drop math_boot

rdcomono Y_math x_math x_read,                    ///
    treatment(D)                                  ///
    generate(y0_math y1_math S_math)              ///
    bandwidth(`bands')                            ///
    kernel(gaussian)                              ///
    folds(5)                                      ///
    order(1)                                      ///
    bootstrap(`B')                                ///
    bootframe(math_boot)                          ///
    bootpoints(`points')

scalar math_band0   = r(band0)
scalar math_band1   = r(band1)
scalar math_q0_band = r(q0_band)
scalar math_q1_band = r(q1_band)
scalar math_N_id    = r(N_supported)


display as text _newline  "Math outcome results"
display as text "  g0 bandwidth: " as result %8.5f scalar(math_band0)
display as text "  g1 bandwidth: " as result %8.5f scalar(math_band1)
display as text "  q0 bandwidth: " as result %8.5f scalar(math_q0_band)
display as text "  q1 bandwidth: " as result %8.5f scalar(math_q1_band)
display as text "  identified observations: " as result %10.0fc scalar(math_N_id)


* Math - tau
generate double tau_math = y1_math - y0_math if S_math == 1
label variable tau_math "Estimated math CATE"

* Math - q plots
rdcomono_qplot,                         ///
    frame(math_boot)                    ///
    level(`level')                      ///
    prefix(math)                        ///
    nocombine

graph export ///
    "output/summer_school/q0_math.png", ///
    name(math_q0) ///
    replace width(1800)

graph export ///
    "output/summer_school/q1_math.png", ///
    name(math_q1) ///
    replace width(1800)

graph export ///
    "output/summer_school/qcompare_math.png", ///
    name(math_compare) ///
    replace width(1800)

* Math - Identified regions
rdcomono_idplot x_math x_read,              ///
    treatment(D)                            ///
    support(S_math)                         ///
    name(math_id)                           ///
    title("Math outcome")

graph export ///
    "output/summer_school/extrapolation_math.png", ///
    name(math_id) ///
    replace width(1800)


** Reading

capture frame drop read_boot

rdcomono Y_read x_math x_read,                    ///
    treatment(D)                                  ///
    generate(y0_read y1_read S_read)              ///
    bandwidth(`bands')                            ///
    kernel(gaussian)                              ///
    folds(5)                                      ///
    order(1)                                      ///
    bootstrap(`B')                                ///
    bootframe(read_boot)                          ///
    bootpoints(`points')


scalar read_band0   = r(band0)
scalar read_band1   = r(band1)
scalar read_q0_band = r(q0_band)
scalar read_q1_band = r(q1_band)
scalar read_N_id    = r(N_supported)


display as text _newline "Reading outcome results"
display as text "  g0 bandwidth: " as result %8.5f scalar(read_band0)
display as text "  g1 bandwidth: " as result %8.5f scalar(read_band1)
display as text "  q0 bandwidth: " as result %8.5f scalar(read_q0_band)
display as text "  q1 bandwidth: " as result %8.5f scalar(read_q1_band)
display as text "  identified observations: " as result %10.0fc scalar(read_N_id)


* Reading - tau
generate double tau_read = y1_read - y0_read if S_read == 1
label variable tau_read "Estimated reading CATE"

* Reading - q plots
rdcomono_qplot,                         ///
    frame(read_boot)                    ///
    level(`level')                      ///
    prefix(read)                        ///
    nocombine

graph export ///
    "output/summer_school/q0_read.png", ///
    name(read_q0) ///
    replace width(1800)

graph export ///
    "output/summer_school/q1_read.png", ///
    name(read_q1) ///
    replace width(1800)

graph export ///
    "output/summer_school/qcompare_read.png", ///
    name(read_compare) ///
    replace width(1800)

* Reading - identified regions
rdcomono_idplot x_math x_read,              ///
    treatment(D)                            ///
    support(S_read)                         ///
    name(read_id)                           ///
    title("Reading outcome")

graph export ///
    "output/summer_school/extrapolation_read.png", ///
    name(read_id) ///
    replace width(1800)


/**********************************************************************
    COUNTERFACTUAL MATH-THRESHOLD POLICIES
**********************************************************************/

quietly levelsof x_math ///
    if x_math >= scalar(cutoff_math), ///
    local(math_thresholds)


tempfile math_policy_results
tempname mathpost

postfile `mathpost'                         ///
    double threshold                        ///
    double effect                           ///
    double number_affected                  ///
    using `math_policy_results', replace



* Temporary policy variables.
generate byte D_cf_math = .
generate double effect_component_math = .

quietly count if S_math == 1
scalar denom_math = r(N)

foreach c of local math_thresholds {

    * Counterfactual policy assignment
    quietly replace D_cf_math = ///
        !(x_math > `c' & ///
          x_read > scalar(cutoff_read))


    * Individual contribution to policy effect.
    quietly replace effect_component_math = ///
          (1-D) * D_cf_math * ///
              (y1_math - Y_math)            ///
        + D * (1-D_cf_math) * ///
              (y0_math - Y_math)            ///
        if S_math == 1


    quietly summarize effect_component_math if S_math == 1, meanonly

    local policy_effect = r(sum) / scalar(denom_math)

    * number whose treatment status changes
    quietly count if S_math == 1 & D_cf_math != D

    local affected = r(N)

    post `mathpost' ///
        (`c') ///
        (`policy_effect') ///
        (`affected')
}

postclose `mathpost'

capture frame drop math_policy
frame create math_policy
frame math_policy: use `math_policy_results', clear
frame math_policy: sort threshold


* Plot the counterfactual policies - math
frame math_policy: twoway                              ///
    (line effect threshold,                            ///
        yaxis(1)                                       ///
        lcolor(blue)                                   ///
        lwidth(medthick))                              ///
    (line number_affected threshold,                   ///
        yaxis(2)                                       ///
        lcolor(red)                                    ///
        lpattern(dash)                                 ///
        lwidth(medthick)),                             ///
    title("Counterfactual math threshold")             ///
    xtitle("Math Score Threshold")                     ///
    ytitle("Estimated policy effect", axis(1))         ///
    ytitle("Number affected", axis(2))                 ///
    yscale(range(0 .025) axis(1))                      ///
    yscale(range(0 15000) axis(2))                     ///
    xscale(range(`cm' 1))                              ///
    legend(order(                                      ///
        1 "Estimated policy effect"                    ///
        2 "Number affected"))                          ///
    name(counterfactual_math, replace)

graph export ///
    "output/summer_school/counterfactual_math.png", ///
    name(counterfactual_math) ///
    replace width(1800)

/**********************************************************************
    COUNTERFACTUAL READING-THRESHOLD POLICIES
**********************************************************************/

quietly levelsof x_read if x_read >= scalar(cutoff_read), local(read_thresholds)

tempfile read_policy_results
tempname readpost

postfile `readpost'                         ///
    double threshold                        ///
    double effect                           ///
    double number_affected                  ///
    using `read_policy_results', replace


generate byte D_cf_read = .
generate double effect_component_read = .


quietly count if S_read == 1
scalar denom_read = r(N)


foreach c of local read_thresholds {

    quietly replace D_cf_read = ///
        !(x_math > scalar(cutoff_math) & ///
          x_read > `c')


    quietly replace effect_component_read = ///
          (1-D) * D_cf_read * ///
              (y1_read - Y_read)            ///
        + D * (1-D_cf_read) * ///
              (y0_read - Y_read)            ///
        if S_read == 1


    quietly summarize ///
        effect_component_read ///
        if S_read == 1, meanonly

    local policy_effect = ///
        r(sum) / scalar(denom_read)


    quietly count if ///
        S_read == 1 & ///
        D_cf_read != D

    local affected = r(N)


    post `readpost' ///
        (`c') ///
        (`policy_effect') ///
        (`affected')
}

postclose `readpost'

capture frame drop read_policy
frame create read_policy
frame read_policy: use `read_policy_results', clear
frame read_policy: sort threshold

* Plot the counterfactual policies - reading
frame read_policy: twoway                              ///
    (line effect threshold,                            ///
        yaxis(1)                                       ///
        lcolor(blue)                                   ///
        lwidth(medthick))                              ///
    (line number_affected threshold,                   ///
        yaxis(2)                                       ///
        lcolor(red)                                    ///
        lpattern(dash)                                 ///
        lwidth(medthick)),                             ///
    title("Counterfactual reading threshold")          ///
    xtitle("Reading Score Threshold")                  ///
    ytitle("Estimated policy effect", axis(1))         ///
    ytitle("Number affected", axis(2))                 ///
    yscale(range(0 .025) axis(1))                      ///
    yscale(range(0 10000) axis(2))                     ///
    xscale(range(`cr' 1))                              ///
    legend(order(                                      ///
        1 "Estimated policy effect"                    ///
        2 "Number affected"))                          ///
    name(counterfactual_read, replace)

graph export ///
    "output/summer_school/counterfactual_reading.png", ///
    name(counterfactual_read) ///
    replace width(1800)


/**********************************************************************
    Final summary
**********************************************************************/
queitly{
noisily display as text _newline "============================================================"
noisily display as text "Summer-school empirical application completed"
noisily display as text "============================================================"
noisily display as text "Math:"
noisily display as text "  identified observations = " as result %10.0fc scalar(math_N_id)
noisily display as text "  first-stage h0 / h1     = " as result %7.4f scalar(math_band0) " / " %7.4f scalar(math_band1)
noisily display as text "Reading:"
noisily display as text "  identified observations = " as result %10.0fc scalar(read_N_id)
noisily display as text "  first-stage h0 / h1     = " as result %7.4f scalar(read_band0) " / " %7.4f scalar(read_band1)
}
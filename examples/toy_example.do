version 17.0
clear all
set more off
set linesize 100

/**********************************************************************
Toy example for rdcomono

This reproduces the data-generating process used in the R package README:

    g1(x1,x2) = sin(x1) + 0.3 sin(x2)
    g0(x1,x2) = 0.8 - 0.8 cos(g1(x1,x2))
    D          = 1{x2 < 0.7 - 0.4 x1}

The example:
    1. simulates the toy data;
    2. estimates the model with rdcomono;
    3. plots the treatment frontier and identified region;
    4. plots q0 and q1 with bootstrap confidence bands.

Run this do-file from the root of the Stata package repository, where the
internal ado-files are stored in the src/ directory.
**********************************************************************/

adopath ++ "`c(pwd)'/src"

capture mkdir "output"


/**********************************************************************
1. Construct the toy data
**********************************************************************/

set seed 1
set obs 1000

generate double x1 = runiform()
generate double x2 = runiform()

generate double mu1 = sin(x1) + 0.3*sin(x2)
generate double mu0 = 0.8 - 0.8*cos(mu1)

generate double frontier = 0.7 - 0.4*x1
generate byte D = x2 < frontier

generate double factual_mean = cond(D == 1, mu1, mu0)
generate double y = factual_mean + rnormal(0, 0.02)

label variable x1 "X1"
label variable x2 "X2"
label variable D  "Treatment"
label variable y  "Observed outcome"
label define treatment_label 0 "Untreated" 1 "Treated"
label values D treatment_label


/**********************************************************************
2. Plot the simulated RDD design
**********************************************************************/

twoway                                                      ///
    (scatter x2 x1 if D == 0,                               ///
        msymbol(Oh) msize(vsmall))                          ///
    (scatter x2 x1 if D == 1,                               ///
        msymbol(+) msize(vsmall))                           ///
    (function y = 0.7 - 0.4*x, range(0 1)                  ///
        lwidth(medthick)),                                  ///
    title("Toy multivariate RDD")                           ///
    subtitle("Treatment rule: X2 < 0.7 - 0.4 X1")          ///
    xtitle("X1")                                           ///
    ytitle("X2")                                           ///
    xscale(range(0 1)) yscale(range(0 1))                  ///
    xlabel(0(.2)1) ylabel(0(.2)1)                          ///
    legend(order(1 "Untreated" 2 "Treated" 3 "Frontier")) ///
    name(toy_design, replace)

graph display toy_design
graph export "output/toy_example/toy_design.png", replace width(1800)


/**********************************************************************
3. Estimate the comonotonic RD model with multiplier bootstrap

The bootstrap uses Exp(1) multiplier weights.  The toy example uses
100 replications to keep execution time manageable.  Increase this to
100 for the number of draws used in the paper.
**********************************************************************/

local boot_reps   100
local boot_points 100
local ci_pct      90

/* Drop the bootstrap frame if this do-file was already run. */
capture frame drop toy_bootstrap

/*
    Separate seed for bootstrap draws so that the multiplier draws are
    reproducible independently of the simulated-data seed.
*/
set seed 24680

rdcomono y x1 x2,                         ///
    treatment(D)                          ///
    generate(y0_hat y1_hat supported)     ///
    bandwidth(0.2 0.3 0.4 0.5 0.6)       ///
    kernel(gaussian)                      ///
    folds(5)                              ///
    order(1)                              ///
    bootstrap(`boot_reps')                ///
    bootframe(toy_bootstrap)              ///
    bootpoints(`boot_points')


/**********************************************************************
4. q0, q1, and combinded plot
**********************************************************************/

rdcomono_qplot, frame(toy_bootstrap)

/**********************************************************************
5. Identified region 
**********************************************************************/

rdcomono_idplot x1 x2,                  ///
    treatment(D)                        ///
    support(supported)                  ///
    name(idplot)

/**********************************************************************
8. Informal accuracy summaries
**********************************************************************/

generate double y0_error = y0_hat - mu0 if supported == 1
generate double y1_error = y1_hat - mu1 if supported == 1

generate double y0_sq_error = y0_error^2 if supported == 1
generate double y1_sq_error = y1_error^2 if supported == 1

quietly summarize y0_sq_error if supported == 1, meanonly
scalar rmse_y0 = sqrt(r(mean))

quietly summarize y1_sq_error if supported == 1, meanonly
scalar rmse_y1 = sqrt(r(mean))

display as text _newline "Informal supported-region accuracy"
display as text "  RMSE for E[Y(0)|X]: " as result %9.5f rmse_y0
display as text "  RMSE for E[Y(1)|X]: " as result %9.5f rmse_y1

display as result _newline "Toy example completed."
display as text "Graphs were saved in the output/ directory:"
display as text "  output/toy_example/toy_design.png"
display as text "  output/toy_example/toy_support.png"
display as text "  output/toy_example/toy_q0.png"
display as text "  output/toy_example/toy_q1.png"
display as text "  output/toy_example/toy_q0_q1.png"
*/
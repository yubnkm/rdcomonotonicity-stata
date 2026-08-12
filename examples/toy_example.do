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

/* Store returned results before another command overwrites r(). */
scalar band0_used = r(band0)
scalar band1_used = r(band1)
scalar q0_band_used = r(q0_band)
scalar q1_band_used = r(q1_band)
scalar N_supported = r(N_supported)

scalar bootstrap_reps_used = r(bootstrap_reps)
scalar bootstrap_points_used = r(bootstrap_points)

local bootstrap_frame "`r(bootstrap_frame)'"

display as text _newline "Selected bandwidths"
display as text "  g0 bandwidth: " as result %6.3f band0_used
display as text "  g1 bandwidth: " as result %6.3f band1_used
display as text "  q0 bandwidth: " as result %6.3f q0_band_used
display as text "  q1 bandwidth: " as result %6.3f q1_band_used
display as text "  supported observations: " ///
    as result %9.0f N_supported

display as text _newline "Bootstrap"
display as text "  replications: " ///
    as result %9.0f bootstrap_reps_used
display as text "  q-grid points: " ///
    as result %9.0f bootstrap_points_used
display as text "  results frame: " ///
    as result "`bootstrap_frame'"

/* Estimated CATE where extrapolation is supported. */
generate double tau_hat = y1_hat - y0_hat if supported == 1
label variable tau_hat "Estimated E[Y(1)-Y(0)|X]"


/**********************************************************************
4. Plot the estimated identification region
**********************************************************************/

twoway                                                          ///
    (scatter x2 x1 if supported == 0,                           ///
        msymbol(Oh) msize(vsmall))                              ///
    (scatter x2 x1 if supported == 1,                           ///
        msymbol(x) msize(vsmall))                               ///
    (function y = 0.7 - 0.4*x, range(0 1)                      ///
        lwidth(medthick)),                                      ///
    title("Estimated extrapolation region")                    ///
    subtitle("Crosses indicate supported observations")        ///
    xtitle("X1")                                               ///
    ytitle("X2")                                               ///
    xscale(range(0 1)) yscale(range(0 1))                      ///
    xlabel(0(.2)1) ylabel(0(.2)1)                              ///
    legend(order(1 "Not supported" 2 "Supported" 3 "Frontier")) ///
    name(toy_support, replace)

graph display toy_support
graph export "output/toy_example/toy_support.png", replace width(1800)


/**********************************************************************
5. Construct q0 plot with 90% pointwise multiplier-bootstrap band

The bootstrap frame contains:

    _rdm_q0_grid       evaluation points y
    _rdm_q0            original q0 estimate
    _rdm_q0_1 ...      bootstrap q0 estimates

At each grid point we calculate

    c(y) = 90th percentile of |q0_boot(y) - q0_hat(y)|

and plot

    q0_hat(y) +/- c(y).
**********************************************************************/

local mainframe "`c(frame)'"

capture frame drop toy_q0_ci
frame copy `bootstrap_frame' toy_q0_ci
frame change toy_q0_ci


/* Keep only rows corresponding to the q0 grid. */
keep if !missing(_rdm_q0_grid)


/*
    Rename the original curve before reshape so that only the
    bootstrap variables match the _rdm_q0_ stub.
*/
rename _rdm_q0_grid q0_input
rename _rdm_q0      q0_estimate

keep q0_input q0_estimate _rdm_q0_*

generate long q0_point = _n


/*
    Convert

        _rdm_q0_1
        _rdm_q0_2
        ...
        _rdm_q0_B

    from wide bootstrap draws to one bootstrap observation per row.
*/
reshape long _rdm_q0_, i(q0_point) j(rep)


/* Absolute bootstrap deviation from the original estimator. */
generate double q0_absdev = ///
    abs(_rdm_q0_ - q0_estimate)


/*
    90th percentile across bootstrap draws at each evaluation point.
*/
bysort q0_point: egen double q0_conf = ///
    pctile(q0_absdev), p(`ci_pct')


/*
    q0_input, q0_estimate, and q0_conf are identical within each
    q0_point after reshape, so retain one observation per point.
*/
bysort q0_point: keep if _n == 1


/* Pointwise bootstrap confidence band. */
generate double q0_lower = q0_estimate - q0_conf
generate double q0_upper = q0_estimate + q0_conf


/* Population q0(y1) from the toy DGP. */
generate double q0_true_curve = ///
    0.8 - 0.8*cos(q0_input)

generate double q0_45_degree = q0_input


label variable q0_input      "E[Y(1)|X]"
label variable q0_estimate   "Estimated q0"
label variable q0_lower      "90% lower band"
label variable q0_upper      "90% upper band"
label variable q0_true_curve "True q0"


twoway                                                        ///
    (rarea q0_lower q0_upper q0_input, sort                  ///
        fcolor(gs10)                                         ///
        fintensity(100)                                      ///
        lcolor(gs10)                                         ///
        lwidth(vthin))                                       ///
    (line q0_estimate q0_input, sort                         ///
        lwidth(medthick))                                    ///
    (line q0_true_curve q0_input, sort                       ///
        lpattern(dash) lwidth(medthick))                     ///
    (line q0_45_degree q0_input, sort                        ///
        lpattern(shortdash)),                                ///
    title("Estimated q0 function")                            ///
    subtitle("90% pointwise multiplier-bootstrap band")       ///
    xtitle("E[Y(1)|X = x]")                                  ///
    ytitle("E[Y(0)|X = x]")                                  ///
    legend(order(                                             ///
        1 "90% bootstrap band"                               ///
        2 "Estimated q0"                                     ///
        3 "True q0"                                          ///
        4 "45-degree line"))                                 ///
    name(q0_plot, replace)

graph display q0_plot
graph export "output/toy_example/toy_q0.png", replace width(1800)


/* Return to the original toy-data frame. */
frame change `mainframe'

/**********************************************************************
6. Construct q1 plot with 90% pointwise multiplier-bootstrap band
**********************************************************************/

capture frame drop toy_q1_ci
frame copy `bootstrap_frame' toy_q1_ci
frame change toy_q1_ci


/* Keep only rows corresponding to the q1 grid. */
keep if !missing(_rdm_q1_grid)


rename _rdm_q1_grid q1_input
rename _rdm_q1      q1_estimate

keep q1_input q1_estimate _rdm_q1_*

generate long q1_point = _n


/* Put bootstrap replications into long format. */
reshape long _rdm_q1_, i(q1_point) j(rep)


/* Absolute bootstrap deviation. */
generate double q1_absdev = ///
    abs(_rdm_q1_ - q1_estimate)


/* 90% pointwise critical value. */
bysort q1_point: egen double q1_conf = ///
    pctile(q1_absdev), p(`ci_pct')


/* Retain one row for each q1 evaluation point. */
bysort q1_point: keep if _n == 1


/* Pointwise bootstrap confidence band. */
generate double q1_lower = q1_estimate - q1_conf
generate double q1_upper = q1_estimate + q1_conf


/*
    Population q1(y0), which is the inverse of

        q0(y1) = 0.8 - 0.8*cos(y1).
*/
generate double q1_acos_argument = ///
    max(-1, min(1, 1 - q1_input/0.8))

generate double q1_true_curve = ///
    acos(q1_acos_argument)

generate double q1_45_degree = q1_input


label variable q1_input      "E[Y(0)|X]"
label variable q1_estimate   "Estimated q1"
label variable q1_lower      "90% lower band"
label variable q1_upper      "90% upper band"
label variable q1_true_curve "True q1"


twoway                                                        ///
    (rarea q1_lower q1_upper q1_input, sort                  ///
        fcolor(gs10)                                         ///
        fintensity(100)                                      ///
        lcolor(gs10)                                         ///
        lwidth(vthin))                                       ///
    (line q1_estimate q1_input, sort                         ///
        lwidth(medthick))                                    ///
    (line q1_true_curve q1_input, sort                       ///
        lpattern(dash) lwidth(medthick))                     ///
    (line q1_45_degree q1_input, sort                        ///
        lpattern(shortdash)),                                ///
    title("Estimated q1 function")                            ///
    subtitle("90% pointwise multiplier-bootstrap band")       ///
    xtitle("E[Y(0)|X = x]")                                  ///
    ytitle("E[Y(1)|X = x]")                                  ///
    legend(order(                                             ///
        1 "90% bootstrap band"                               ///
        2 "Estimated q1"                                     ///
        3 "True q1"                                          ///
        4 "45-degree line"))                                 ///
    name(q1_plot, replace)

graph display q1_plot
graph export "output/toy_example/toy_q1.png", replace width(1800)


frame change `mainframe'


/**********************************************************************
7. Show q0 and q1 together
**********************************************************************/

graph combine q0_plot q1_plot,                               ///
    cols(2)                                                   ///
    title("Comonotonic mappings in the toy example")          ///
    name(q_plots, replace)

graph display q_plots
graph export "output/toy_example/toy_q0_q1.png", replace width(2400)


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

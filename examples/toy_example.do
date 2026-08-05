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
    4. plots q0 and q1 without bootstrap confidence bands.

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
graph export "output/toy_design.png", replace width(1800)


/**********************************************************************
3. Estimate the comonotonic RD model

The bandwidth candidates match the R README:
    seq(0.2, 0.6, length.out = 5)
**********************************************************************/

rdcomono y x1 x2,                         ///
    treatment(D)                          ///
    generate(y0_hat y1_hat supported)     ///
    bandwidth(0.2 0.3 0.4 0.5 0.6)       ///
    kernel(gaussian)                      ///
    folds(5)                              ///
    order(1)

/* Store returned results before another command overwrites r(). */
scalar band0_used = r(band0)
scalar band1_used = r(band1)
scalar q0_band_used = r(q0_band)
scalar q1_band_used = r(q1_band)
scalar N_supported = r(N_supported)

display as text _newline "Selected bandwidths"
display as text "  g0 bandwidth: " as result %6.3f band0_used
display as text "  g1 bandwidth: " as result %6.3f band1_used
display as text "  q0 bandwidth: " as result %6.3f q0_band_used
display as text "  q1 bandwidth: " as result %6.3f q1_band_used
display as text "  supported observations: " as result %9.0f N_supported

/* Estimated conditional average treatment effect where it is supported. */
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
graph export "output/toy_support.png", replace width(1800)


/**********************************************************************
5. Construct the q0 plot

For treated observations:

    y1_hat = estimated factual E[Y(1)|X]
    y0_hat = q0(y1_hat)

Therefore, plotting y0_hat against y1_hat over supported treated
observations displays the estimated q0 function.
**********************************************************************/

generate double q0_input = y1_hat if D == 1 & supported == 1
generate double q0_estimate = y0_hat if D == 1 & supported == 1

/* Population q0(y1) used by the toy DGP. */
generate double q0_true_curve = 0.8 - 0.8*cos(q0_input)          ///
    if !missing(q0_input)
generate double q0_45_degree = q0_input if !missing(q0_input)

label variable q0_input "E[Y(1)|X]"
label variable q0_estimate "Estimated q0"
label variable q0_true_curve "True q0"

twoway                                                         ///
    (line q0_estimate q0_input, sort lwidth(medthick))         ///
    (line q0_true_curve q0_input, sort lpattern(dash)          ///
        lwidth(medthick))                                      ///
    (line q0_45_degree q0_input, sort lpattern(shortdash)),    ///
    title("Estimated q0 function")                            ///
    subtitle("No bootstrap confidence band")                  ///
    xtitle("E[Y(1)|X = x]")                                  ///
    ytitle("E[Y(0)|X = x]")                                  ///
    legend(order(1 "Estimated q0" 2 "True q0"               ///
        3 "45-degree line"))                                  ///
    name(q0_plot, replace)

graph display q0_plot
graph export "output/toy_q0.png", replace width(1800)


/**********************************************************************
6. Construct the q1 plot

For untreated observations:

    y0_hat = estimated factual E[Y(0)|X]
    y1_hat = q1(y0_hat)

Therefore, plotting y1_hat against y0_hat over supported untreated
observations displays the estimated q1 function.
**********************************************************************/

generate double q1_input = y0_hat if D == 0 & supported == 1
generate double q1_estimate = y1_hat if D == 0 & supported == 1

/*
Population q1(y0), the inverse of q0. The inner max/min operation prevents
small numerical errors from sending the argument of acos() outside [-1,1].
*/
generate double q1_acos_argument =                         ///
    max(-1, min(1, 1 - q1_input/0.8))                     ///
    if !missing(q1_input)
generate double q1_true_curve = acos(q1_acos_argument)    ///
    if !missing(q1_acos_argument)
generate double q1_45_degree = q1_input if !missing(q1_input)

label variable q1_input "E[Y(0)|X]"
label variable q1_estimate "Estimated q1"
label variable q1_true_curve "True q1"

twoway                                                         ///
    (line q1_estimate q1_input, sort lwidth(medthick))         ///
    (line q1_true_curve q1_input, sort lpattern(dash)          ///
        lwidth(medthick))                                      ///
    (line q1_45_degree q1_input, sort lpattern(shortdash)),    ///
    title("Estimated q1 function")                            ///
    subtitle("No bootstrap confidence band")                  ///
    xtitle("E[Y(0)|X = x]")                                  ///
    ytitle("E[Y(1)|X = x]")                                  ///
    legend(order(1 "Estimated q1" 2 "True q1"               ///
        3 "45-degree line"))                                  ///
    name(q1_plot, replace)

graph display q1_plot
graph export "output/toy_q1.png", replace width(1800)


/**********************************************************************
7. Show q0 and q1 together
**********************************************************************/

graph combine q0_plot q1_plot,                               ///
    cols(2)                                                   ///
    title("Comonotonic mappings in the toy example")          ///
    name(q_plots, replace)

graph display q_plots
graph export "output/toy_q0_q1.png", replace width(2400)


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
display as text "  output/toy_design.png"
display as text "  output/toy_support.png"
display as text "  output/toy_q0.png"
display as text "  output/toy_q1.png"
display as text "  output/toy_q0_q1.png"

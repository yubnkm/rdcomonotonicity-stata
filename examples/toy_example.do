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
    Construct the toy example data
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
    Plot the simulated RDD design
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
    Estimate the comonotonic RD model with multiplier bootstrap
**********************************************************************/

local boot_reps   100
local boot_points 100
local ci_pct      90

/* Drop the bootstrap frame if this do-file was already run. */
capture frame drop toy_bootstrap

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
    q0, q1, and combinded plot
**********************************************************************/

rdcomono_qplot, frame(toy_bootstrap)

/**********************************************************************
    Identified region 
**********************************************************************/

rdcomono_idplot x1 x2,                  ///
    treatment(D)                        ///
    support(supported)                  ///
    name(idplot)


/**********************************************************************
    Counterfactual policy effect
    - Original policy: D = 1{x2 < 0.7 - 0.4*x1}
    - Counterfactual policy: D_cf = 1{x2 < 0.7 - 0.4*x1 + 0.05}
**********************************************************************/

local delta = 0.05
generate byte D_counterfactual = x2 < (0.7 - 0.4*x1 + `delta')
label variable D_counterfactual "Treatment under counterfactual policy"

rdcomono_policy y,                            ///
    treatment(D)                              ///
    policy(D_counterfactual)                  ///
    y0(y0_hat)                                ///
    y1(y1_hat)                                ///
    support(supported)                        ///
    bootframe(toy_bootstrap)                  ///
    level(90)
    
scalar policy_estimate = r(estimate)
scalar policy_ci_low = r(conf_low)
scalar policy_ci_high = r(conf_high)
scalar policy_identified_n = r(N_supported)
scalar policy_num_affected = r(num_affected)
scalar policy_affected_share = r(affected_share)

* True counterfactual policy effect to compare with the estimate
generate double true_counterfactual_mean = cond(D_counterfactual == 1, mu1, mu0)
generate double true_policy_change = true_counterfactual_mean - factual_mean if supported == 1
quietly summarize true_policy_change if supported == 1, meanonly
scalar true_policy_effect = r(mean)

* Display comparison
quietly {
	noisily display as text "{text}{hline 62}"
	noisily display as text "Counterfactual policy: frontier shifted upward by 0.05"
	noisily display as text "{text}{hline 62}"
	
    noisily display as text "True policy effect"                 _col(35) ": " as result %10.6f true_policy_effect
    noisily display as text "Estimated policy effect"            _col(35) ": " as result %10.6f policy_estimate
    noisily display as text "90% bootstrap confidence interval"  _col(35) ": [" as result %10.6f policy_ci_low as text ", " as result %10.6f policy_ci_high as text "]"
    
    noisily display as text _newline "Observations with S = 1"     _col(35) ": " as result %10.0f policy_identified_n
    noisily display as text "Number affected by policy"          _col(35) ": " as result %10.0f policy_num_affected
    noisily display as text "Affected share among S = 1"         _col(35) ": " as result %10.4f policy_affected_share
	
	noisily display as text "{text}{hline 62}"

}

version 17.0
clear all
set more off

/*
    Run this file from the repository root. The repository is assumed to have

        src/_rdcomono_localpoly.ado
        src/rdcomono.ado
*/
adopath ++ "`c(pwd)'/src"


/********************************************************************
Test 1: exact linear comonotonic DGP
********************************************************************/

clear
set seed 20260805
set obs 800

generate double x1 = runiform()
generate double x2 = runiform()

/* Sharp multivariate treatment rule. */
generate byte D = x2 < 0.70 - 0.40*x1

/*
    Both potential-outcome means are increasing functions of the same index.
    Therefore q1(y0) = 1 + 2*y0 and q0(y1) = (y1 - 1)/2.
*/
generate double y0_true = 0.5 + x1 + 0.5*x2
generate double y1_true = 1 + 2*y0_true
generate double y = D*y1_true + (1-D)*y0_true

rdcomono y x1 x2,                    ///
    treatment(D)                     ///
    generate(y0_hat y1_hat supported) ///
    bandwidth(0.30)                  ///
    kernel(gaussian)                 ///
    folds(5)                         ///
    order(1)

scalar rd_N = r(N)
scalar rd_N0 = r(N0)
scalar rd_N1 = r(N1)
scalar rd_N_frontier0 = r(N_frontier0)
scalar rd_N_frontier1 = r(N_frontier1)
scalar rd_N_supported = r(N_supported)

assert rd_N == 800
assert rd_N0 > 0
assert rd_N1 > 0
assert rd_N_frontier0 > 0
assert rd_N_frontier1 > 0
assert rd_N_supported > 0

assert inlist(supported, 0, 1)

/* Factual conditional means must be reproduced everywhere. */
assert abs(y0_hat - y0_true) < 1e-7 if D == 0
assert abs(y1_hat - y1_true) < 1e-7 if D == 1

/* Counterfactual means must be reproduced on the estimated support. */
assert abs(y0_hat - y0_true) < 1e-6 if D == 1 & supported == 1
assert abs(y1_hat - y1_true) < 1e-6 if D == 0 & supported == 1

/* Hence the supported CATE is also recovered. */
generate double tau_true = y1_true - y0_true
generate double tau_hat = y1_hat - y0_hat
assert abs(tau_hat - tau_true) < 1e-6 if supported == 1

display as result ///
    "Test 1 passed: exact linear comonotonic DGP is recovered."


/********************************************************************
Test 2: unit numerical weights reproduce the unweighted estimator
********************************************************************/

generate double unit_weight = 1

rdcomono y x1 x2,                           ///
    treatment(D)                            ///
    generate(y0_weighted y1_weighted S_w)   ///
    bandwidth(0.30)                         ///
    wvar(unit_weight)                       ///
    kernel(gaussian)                        ///
    folds(5)                                ///
    order(1)

assert abs(y0_weighted - y0_hat) < 1e-10 ///
    if !missing(y0_weighted, y0_hat)
assert abs(y1_weighted - y1_hat) < 1e-10 ///
    if !missing(y1_weighted, y1_hat)
assert S_w == supported

display as result ///
    "Test 2 passed: unit weights reproduce the unweighted estimator."


/********************************************************************
Test 3: candidate bandwidths invoke cross-validation
********************************************************************/

set seed 13579

rdcomono y x1 x2 if _n <= 300,          ///
    treatment(D)                        ///
    generate(y0_cv y1_cv S_cv)          ///
    bandwidth(0.20 0.30 0.40)           ///
    kernel(gaussian)                    ///
    folds(3)                            ///
    order(1)

scalar cv_band0 = r(band0)
scalar cv_band1 = r(band1)
scalar cv_q0_band = r(q0_band)
scalar cv_q1_band = r(q1_band)

assert inlist(cv_band0, 0.20, 0.30, 0.40)
assert inlist(cv_band1, 0.20, 0.30, 0.40)
assert inlist(cv_q0_band, 0.20, 0.30, 0.40)
assert inlist(cv_q1_band, 0.20, 0.30, 0.40)
assert inlist(S_cv, 0, 1) if _n <= 300
assert missing(y0_cv, y1_cv, S_cv) if _n > 300

display as result ///
    "Test 3 passed: bandwidth cross-validation returns candidates."


/********************************************************************
Test 4: incomplete observations are excluded and outputs remain missing
********************************************************************/

replace y = . in 1/5
replace x1 = . in 6/10

rdcomono y x1 x2,                              ///
    treatment(D)                               ///
    generate(y0_missing y1_missing S_missing)  ///
    bandwidth(0.30)                            ///
    kernel(gaussian)                           ///
    order(1)

scalar missing_N = r(N)
assert missing_N == 790
assert missing(y0_missing, y1_missing, S_missing) in 1/10
assert !missing(S_missing) in 11/L

display as result ///
    "Test 4 passed: complete-case handling works."


/********************************************************************
Test 5: invalid treatment values are rejected
********************************************************************/

replace D = 2 in 11

capture noisily rdcomono y x1 x2,              ///
    treatment(D)                               ///
    generate(should_not_y0 should_not_y1 should_not_S) ///
    bandwidth(0.30)

assert _rc == 198

display as result ///
    "Test 5 passed: nonbinary treatment is rejected."


display as result "All basic rdcomono tests passed."

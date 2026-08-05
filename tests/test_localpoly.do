version 17.0
clear all
set more off

/*
    Run this file from the repository root.
*/
adopath ++ "`c(pwd)'/src"


/********************************************************************
Test 1: order(1) reproduces an exactly linear function
********************************************************************/

clear
set seed 12345
set obs 500

generate double x1 = runiform()
generate double x2 = runiform()
generate double y  = 1 + 2*x1 - 3*x2

_rdcomono_localpoly y x1 x2,    ///
    at(x1 x2)                   ///
    generate(yhat)              ///
    bandwidth(0.30)             ///
    kernel(gaussian)            ///
    order(1)

generate double error_linear = abs(yhat - y)

summarize error_linear, detail
assert error_linear < 1e-8 if !missing(error_linear)

display as result ///
    "Test 1 passed: order(1) reproduces a linear function."


/********************************************************************
Test 2: unit observation weights reproduce the unweighted result
********************************************************************/

generate double user_weight = 1

_rdcomono_localpoly y x1 x2,    ///
    at(x1 x2)                   ///
    generate(yhat_weighted)     ///
    bandwidth(0.30)             ///
    kernel(gaussian)            ///
    order(1)                    ///
    wvar(user_weight)

generate double error_unit_weight = ///
    abs(yhat_weighted - yhat)

assert error_unit_weight < 1e-10 if !missing(error_unit_weight)

display as result ///
    "Test 2 passed: unit weights reproduce the unweighted fit."


/********************************************************************
Test 3: separate evaluation and center variables
********************************************************************/

generate double x1_fit = x1 + 0.01
generate double x2_fit = x2 - 0.01
generate double y_shifted_true = 1 + 2*x1_fit - 3*x2_fit

_rdcomono_localpoly y x1 x2,    ///
    at(x1_fit x2_fit)           ///
    center(x1 x2)               ///
    generate(yhat_shifted)      ///
    bandwidth(0.30)             ///
    kernel(gaussian)            ///
    order(1)

generate double error_shifted = ///
    abs(yhat_shifted - y_shifted_true)

assert error_shifted < 1e-8 if !missing(error_shifted)

display as result ///
    "Test 3 passed: at() and center() can differ."


/********************************************************************
Test 4: all three kernels work for order(1)
********************************************************************/

foreach kernel in gaussian uniform triangular {

    _rdcomono_localpoly y x1 x2,         ///
        at(x1 x2)                        ///
        generate(fit_`kernel')           ///
        bandwidth(0.30)                  ///
        kernel(`kernel')                 ///
        order(1)

    generate double error_`kernel' = ///
        abs(fit_`kernel' - y)

    assert error_`kernel' < 1e-8 ///
        if !missing(error_`kernel')

    display as result ///
        "Test 4 passed for kernel: `kernel'"
}


/********************************************************************
Test 5: one candidate bandwidth is returned unchanged
********************************************************************/

_rdcomono_localpoly y x1 x2,         ///
    at(x1 x2)                        ///
    generate(yhat_single_band)       ///
    bandwidth(0.30)                  ///
    folds(5)                         ///
    kernel(gaussian)                 ///
    order(1)

assert abs(r(bandwidth) - 0.30) < 1e-12

display as result ///
    "Test 5 passed: one bandwidth is used directly."


/********************************************************************
Test 6: CV selects one candidate and is reproducible
********************************************************************/

set seed 987654

_rdcomono_localpoly y x1 x2,         ///
    at(x1 x2)                        ///
    generate(yhat_cv)                ///
    bandwidth(0.10 0.20 0.30 0.40)  ///
    folds(5)                         ///
    kernel(gaussian)                 ///
    order(1)

scalar selected_h = r(bandwidth)

assert ///
    selected_h == 0.10 | ///
    selected_h == 0.20 | ///
    selected_h == 0.30 | ///
    selected_h == 0.40

set seed 987654

_rdcomono_localpoly y x1 x2,         ///
    at(x1 x2)                        ///
    generate(yhat_cv_repeat)         ///
    bandwidth(0.10 0.20 0.30 0.40)  ///
    folds(5)                         ///
    kernel(gaussian)                 ///
    order(1)

assert abs(r(bandwidth) - selected_h) < 1e-12
assert abs(yhat_cv_repeat - yhat_cv) < 1e-10 ///
    if !missing(yhat_cv_repeat, yhat_cv)

display as result ///
    "Test 6 passed: CV selects a candidate reproducibly."


/********************************************************************
Test 7: weighted cross-validation
********************************************************************/

generate double nonconstant_weight = 0.5 + x1 + x2

set seed 246810

_rdcomono_localpoly y x1 x2,         ///
    at(x1 x2)                        ///
    generate(yhat_cv_weighted)       ///
    bandwidth(0.10 0.20 0.30 0.40)  ///
    folds(5)                         ///
    kernel(gaussian)                 ///
    order(1)                         ///
    wvar(nonconstant_weight)

generate double error_cv_weighted = ///
    abs(yhat_cv_weighted - y)

assert error_cv_weighted < 1e-8 ///
    if !missing(error_cv_weighted)

display as result ///
    "Test 7 passed: weighted CV works."


/********************************************************************
Test 8: order(0) reproduces a constant function
********************************************************************/

generate double y_constant = 7.25

_rdcomono_localpoly y_constant x1 x2, ///
    at(x1 x2)                         ///
    generate(yhat_constant)           ///
    bandwidth(0.30)                   ///
    kernel(gaussian)                  ///
    order(0)

generate double error_constant = ///
    abs(yhat_constant - y_constant)

assert error_constant < 1e-10 ///
    if !missing(error_constant)

display as result ///
    "Test 8 passed: order(0) reproduces a constant function."


/********************************************************************
Test 9: order(2) reproduces a total-degree quadratic function
********************************************************************/

generate double y_quadratic = ///
      1                         ///
    + 2*x1                      ///
    - 3*x2                      ///
    + 0.5*x1^2                  ///
    - 0.75*x2^2                 ///
    + 1.25*x1*x2

_rdcomono_localpoly y_quadratic x1 x2, ///
    at(x1 x2)                          ///
    generate(yhat_quadratic)           ///
    bandwidth(0.35)                    ///
    kernel(gaussian)                   ///
    order(2)

generate double error_quadratic = ///
    abs(yhat_quadratic - y_quadratic)

summarize error_quadratic, detail
assert error_quadratic < 1e-8 ///
    if !missing(error_quadratic)

display as result ///
    "Test 9 passed: order(2) reproduces a quadratic with interaction."


/********************************************************************
Test 10: order(2) works with separate evaluation and center points
********************************************************************/

generate double y_quadratic_shifted = ///
      1                                      ///
    + 2*x1_fit                               ///
    - 3*x2_fit                               ///
    + 0.5*x1_fit^2                           ///
    - 0.75*x2_fit^2                          ///
    + 1.25*x1_fit*x2_fit

_rdcomono_localpoly y_quadratic x1 x2, ///
    at(x1_fit x2_fit)                   ///
    center(x1 x2)                       ///
    generate(yhat_quadratic_shifted)    ///
    bandwidth(0.35)                     ///
    kernel(gaussian)                    ///
    order(2)

generate double error_quadratic_shifted = ///
    abs(yhat_quadratic_shifted - y_quadratic_shifted)

assert error_quadratic_shifted < 1e-8 ///
    if !missing(error_quadratic_shifted)

display as result ///
    "Test 10 passed: higher-order prediction works when at() differs from center()."


/********************************************************************
Test 11: order(3) reproduces a cubic function
********************************************************************/

generate double y_cubic = ///
      0.5                       ///
    + x1                        ///
    - 2*x2                      ///
    + 0.75*x1*x2                ///
    + 0.4*x1^2                  ///
    - 0.2*x2^2                  ///
    + 0.3*x1^3                  ///
    - 0.6*x1^2*x2               ///
    + 0.25*x1*x2^2              ///
    + 0.1*x2^3

_rdcomono_localpoly y_cubic x1 x2, ///
    at(x1 x2)                      ///
    generate(yhat_cubic)           ///
    bandwidth(0.45)                ///
    kernel(gaussian)               ///
    order(3)

generate double error_cubic = ///
    abs(yhat_cubic - y_cubic)

summarize error_cubic, detail
assert error_cubic < 1e-7 ///
    if !missing(error_cubic)

display as result ///
    "Test 11 passed: order(3) reproduces a cubic function."


/********************************************************************
Test 12: cross-validation works with order(2)
********************************************************************/

set seed 13579

_rdcomono_localpoly y_quadratic x1 x2, ///
    at(x1 x2)                          ///
    generate(yhat_quadratic_cv)        ///
    bandwidth(0.20 0.30 0.40 0.50)    ///
    folds(5)                           ///
    kernel(gaussian)                   ///
    order(2)

scalar selected_h_quadratic = r(bandwidth)

assert ///
    selected_h_quadratic == 0.20 | ///
    selected_h_quadratic == 0.30 | ///
    selected_h_quadratic == 0.40 | ///
    selected_h_quadratic == 0.50

generate double error_quadratic_cv = ///
    abs(yhat_quadratic_cv - y_quadratic)

assert error_quadratic_cv < 1e-7 ///
    if !missing(error_quadratic_cv)

display as result ///
    "Test 12 passed: CV works with a higher-order basis."


/********************************************************************
Test 13: incomplete training observations are excluded
********************************************************************/

generate double y_missing = y
replace y_missing = . in 1/10

_rdcomono_localpoly y_missing x1 x2, ///
    at(x1 x2)                        ///
    generate(yhat_missing_training)  ///
    bandwidth(0.30)                  ///
    kernel(gaussian)                 ///
    order(1)

scalar n_training = r(n)

generate double error_missing_training = ///
    abs(yhat_missing_training - y)

assert error_missing_training < 1e-8 ///
    if !missing(error_missing_training)

assert n_training == 490

display as result ///
    "Test 13 passed: incomplete training observations are excluded."


/********************************************************************
Test 14: incomplete evaluation rows remain missing
********************************************************************/

generate double x1_eval_missing = x1
replace x1_eval_missing = . in 1/5

_rdcomono_localpoly y x1 x2,          ///
    at(x1_eval_missing x2)            ///
    generate(yhat_eval_missing)       ///
    bandwidth(0.30)                   ///
    kernel(gaussian)                  ///
    order(1)

assert missing(yhat_eval_missing) in 1/5
assert !missing(yhat_eval_missing) in 6/L

display as result ///
    "Test 14 passed: incomplete evaluation rows remain missing."


/********************************************************************
Test 15: invalid order and invalid weights are rejected
********************************************************************/

capture noisily _rdcomono_localpoly y x1 x2, ///
    at(x1 x2)                               ///
    generate(should_not_exist_order)        ///
    bandwidth(0.30)                         ///
    order(-1)

assert _rc == 198

generate double negative_weight = 1
replace negative_weight = -1 in 1

capture noisily _rdcomono_localpoly y x1 x2, ///
    at(x1 x2)                               ///
    generate(should_not_exist_weight)       ///
    bandwidth(0.30)                         ///
    wvar(negative_weight)

assert _rc == 198

display as result ///
    "Test 15 passed: invalid inputs are rejected."


display as result "All local-polynomial tests passed."

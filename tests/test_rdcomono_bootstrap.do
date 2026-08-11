version 17.0
clear all
set more off

/*
    Run from repository root.
*/
adopath ++ "`c(pwd)'/src"


/********************************************************************
Test 1: bootstrap runs and returns expected objects
********************************************************************/

clear
capture frame drop boot_exact

set seed 20260811
set obs 500

generate double x1 = runiform()
generate double x2 = runiform()

generate byte D = x2 < 0.70 - 0.40*x1

/*
    Exact linear comonotonic DGP:

        y0 = 0.5 + x1 + 0.5*x2
        y1 = 1 + 2*y0

    Therefore

        q1(y0) = 1 + 2*y0
        q0(y1) = (y1 - 1)/2
*/
generate double y0_true = ///
    0.5 + x1 + 0.5*x2

generate double y1_true = ///
    1 + 2*y0_true

generate double y = ///
    D*y1_true + (1-D)*y0_true


set seed 24680

rdcomono y x1 x2,                         ///
    treatment(D)                          ///
    generate(y0_hat y1_hat supported)     ///
    bandwidth(0.30)                       ///
    kernel(gaussian)                      ///
    folds(5)                              ///
    order(1)                              ///
    bootstrap(3)                          ///
    bootframe(boot_exact)                 ///
    bootpoints(7)


scalar boot_reps   = r(bootstrap_reps)
scalar boot_points = r(bootstrap_points)

matrix boot_q0_bands = r(bootstrap_q0_bands)
matrix boot_q1_bands = r(bootstrap_q1_bands)

assert boot_reps == 3
assert boot_points == 7

if colsof(boot_q0_bands) != 3 {
    display as error ///
        "q0 bootstrap bandwidth matrix has wrong size"
    exit 9
}

if colsof(boot_q1_bands) != 3 {
    display as error ///
        "q1 bootstrap bandwidth matrix has wrong size"
    exit 9
}


/*
    Bootstrap frame should contain one row per original observation
    because N = 500 > bootpoints = 7.
*/
frame boot_exact: quietly count

if r(N) != 500 {
    display as error ///
        "bootstrap frame has incorrect number of rows"
    exit 9
}


/*
    Exactly seven q-grid points.
*/
frame boot_exact: ///
    assert !missing(_rdm_q0_grid) in 1/7

frame boot_exact: ///
    assert missing(_rdm_q0_grid) in 8/L

frame boot_exact: ///
    assert !missing(_rdm_q1_grid) in 1/7

frame boot_exact: ///
    assert missing(_rdm_q1_grid) in 8/L


display as result ///
    "Bootstrap Test 1 passed: expected objects were created."


/********************************************************************
Test 2: exponential multipliers are valid
********************************************************************/

frame boot_exact: ///
    assert _rdm_m1 >= 0 ///
    if !missing(_rdcomono_id)

frame boot_exact: ///
    assert _rdm_m2 >= 0 ///
    if !missing(_rdcomono_id)

frame boot_exact: ///
    assert _rdm_m3 >= 0 ///
    if !missing(_rdcomono_id)


frame boot_exact: ///
    quietly summarize _rdm_m1, meanonly

if r(max) <= r(min) {
    display as error ///
        "bootstrap multipliers do not vary"
    exit 9
}


display as result ///
    "Bootstrap Test 2 passed: Exp(1) multipliers are valid."


/********************************************************************
Test 3: exact linear DGP remains exact under multiplier weighting
********************************************************************/

/*
    Reconstruct true functions in bootstrap frame.
*/
frame boot_exact: ///
    generate double truth_y0 = ///
    0.5 + x1 + 0.5*x2

frame boot_exact: ///
    generate double truth_y1 = ///
    1 + 2*truth_y0


/*
    Multiplier weights should not change an exactly specified local-linear
    regression, apart from numerical error.
*/
forvalues b = 1/3 {

    frame boot_exact: ///
        assert abs(_rdm_y0_`b' - truth_y0) < 1e-6 ///
        if !missing(_rdcomono_id)

    frame boot_exact: ///
        assert abs(_rdm_y1_`b' - truth_y1) < 1e-6 ///
        if !missing(_rdcomono_id)
}


/*
    q0(y1) = (y1 - 1)/2
*/
frame boot_exact: ///
    assert abs( ///
        _rdm_q0 - ///
        ((_rdm_q0_grid - 1)/2) ///
    ) < 1e-6 ///
    if !missing(_rdm_q0_grid)


forvalues b = 1/3 {

    frame boot_exact: ///
        assert abs( ///
            _rdm_q0_`b' - ///
            ((_rdm_q0_grid - 1)/2) ///
        ) < 1e-6 ///
        if !missing(_rdm_q0_grid)
}


/*
    q1(y0) = 1 + 2*y0
*/
frame boot_exact: ///
    assert abs( ///
        _rdm_q1 - ///
        (1 + 2*_rdm_q1_grid) ///
    ) < 1e-6 ///
    if !missing(_rdm_q1_grid)


forvalues b = 1/3 {

    frame boot_exact: ///
        assert abs( ///
            _rdm_q1_`b' - ///
            (1 + 2*_rdm_q1_grid) ///
        ) < 1e-6 ///
        if !missing(_rdm_q1_grid)
}


display as result ///
    "Bootstrap Test 3 passed: exact linear DGP is preserved."


/********************************************************************
Test 4: bootstrap does NOT alter the original S
********************************************************************/

/*
    Compare original support vector with the copy stored in the bootstrap
    result.
*/
mkmat supported, matrix(S_original)

frame boot_exact: ///
    mkmat supported ///
    if !missing(_rdcomono_id), ///
    matrix(S_bootstrap)

scalar S_difference = ///
    mreldif(S_original, S_bootstrap)

if S_difference > 1e-14 {
    display as error ///
        "bootstrap modified the original support indicator"
    exit 9
}


display as result ///
    "Bootstrap Test 4 passed: original S is held fixed."


/********************************************************************
Test 5: reproducibility under set seed
********************************************************************/

capture frame drop boot_repeat

set seed 24680

rdcomono y x1 x2,                         ///
    treatment(D)                          ///
    generate(y0_hat_repeat                ///
             y1_hat_repeat                ///
             supported_repeat)            ///
    bandwidth(0.30)                       ///
    kernel(gaussian)                      ///
    folds(5)                              ///
    order(1)                              ///
    bootstrap(3)                          ///
    bootframe(boot_repeat)                ///
    bootpoints(7)


frame boot_exact: ///
    mkmat ///
    _rdm_m1 ///
    _rdm_y0_1 ///
    _rdm_y1_1 ///
    if !missing(_rdcomono_id), ///
    matrix(Boot_A)


frame boot_repeat: ///
    mkmat ///
    _rdm_m1 ///
    _rdm_y0_1 ///
    _rdm_y1_1 ///
    if !missing(_rdcomono_id), ///
    matrix(Boot_B)


scalar reproducibility_error = ///
    mreldif(Boot_A, Boot_B)

if reproducibility_error > 1e-12 {
    display as error ///
        "bootstrap is not reproducible under the same set seed"
    exit 9
}


display as result ///
    "Bootstrap Test 5 passed: set seed gives reproducible draws."


/********************************************************************
Test 6: noisy DGP produces non-degenerate bootstrap variation
********************************************************************/

clear
capture frame drop boot_noisy

set seed 98765
set obs 600

generate double x1 = runiform()
generate double x2 = runiform()

generate byte D = ///
    x2 < 0.70 - 0.40*x1

generate double y0_true = ///
    0.5 + x1 + 0.5*x2

generate double y1_true = ///
    1 + 2*y0_true

generate double y = ///
    D*y1_true + ///
    (1-D)*y0_true + ///
    rnormal(0, 0.10)


set seed 13579

rdcomono y x1 x2,                         ///
    treatment(D)                          ///
    generate(y0_noise y1_noise S_noise)   ///
    bandwidth(0.30)                       ///
    kernel(gaussian)                      ///
    folds(5)                              ///
    order(1)                              ///
    bootstrap(3)                          ///
    bootframe(boot_noisy)                 ///
    bootpoints(7)


frame boot_noisy: ///
    generate double bootstrap_difference = ///
    abs(_rdm_y0_1 - _rdm_y0_2) + ///
    abs(_rdm_y1_1 - _rdm_y1_2) ///
    if !missing(_rdcomono_id)


frame boot_noisy: ///
    quietly summarize bootstrap_difference, meanonly

scalar max_bootstrap_difference = r(max)

if max_bootstrap_difference <= 1e-10 {
    display as error ///
        "noisy bootstrap estimates do not vary across draws"
    exit 9
}


display as result ///
    "Bootstrap Test 6 passed: noisy estimates vary across draws."


/********************************************************************
Test 7: original user weights are retained and multiplied by multipliers
********************************************************************/

clear
capture frame drop boot_weighted

set seed 112233
set obs 500

generate double x1 = runiform()
generate double x2 = runiform()

generate byte D = ///
    x2 < 0.70 - 0.40*x1

generate double y0_true = ///
    0.5 + x1 + 0.5*x2

generate double y1_true = ///
    1 + 2*y0_true

generate double y = ///
    D*y1_true + ///
    (1-D)*y0_true + ///
    rnormal(0, 0.05)

/*
    Nonconstant original weights.
*/
generate double user_weight = ///
    0.5 + x1


set seed 445566

rdcomono y x1 x2,                         ///
    treatment(D)                          ///
    generate(y0_wboot y1_wboot S_wboot)   ///
    bandwidth(0.30)                       ///
    wvar(user_weight)                     ///
    kernel(gaussian)                      ///
    folds(5)                              ///
    order(1)                              ///
    bootstrap(2)                          ///
    bootframe(boot_weighted)              ///
    bootpoints(7)


/*
    The bootstrap frame preserves w_i as _rdm_basew.
*/
frame boot_weighted: ///
    assert abs(_rdm_basew - user_weight) < 1e-14 ///
    if !missing(_rdcomono_id)

/*
    Hence the actual weight used by draw 1 is reconstructible as
        user_weight * _rdm_m1.
*/
frame boot_weighted: ///
    generate double reconstructed_weight = ///
    _rdm_basew * _rdm_m1 ///
    if !missing(_rdcomono_id)

frame boot_weighted: ///
    assert reconstructed_weight >= 0 ///
    if !missing(_rdcomono_id)


display as result ///
    "Bootstrap Test 7 passed: original weights are combined with multipliers."


/********************************************************************
Test 8: invalid bootstrap repetitions are rejected
********************************************************************/

capture noisily rdcomono y x1 x2,         ///
    treatment(D)                          ///
    generate(bad_y0 bad_y1 bad_S)         ///
    bandwidth(0.30)                       ///
    bootstrap(-1)

assert _rc == 198


display as result ///
    "Bootstrap Test 8 passed: invalid bootstrap() is rejected."


display as result _newline ///
    "All rdcomono multiplier-bootstrap tests passed."
*! version 0.3.0 05aug2026
program define _rdcomono_localpoly, rclass
    version 17.0

    /*
        Internal multivariate local-polynomial regression.

        Features:
          - polynomial orders 0, 1, 2, ...
          - total-degree multivariate polynomial basis
          - fixed or cross-validated bandwidth
          - Gaussian, uniform, and triangular kernels
          - optional nonnegative numerical observation weights
          - separate evaluation points and kernel centers
          - missing-value handling
    */

    syntax varlist(min=2 numeric) [if] [in],                 ///
        AT(varlist numeric)                                 ///
        GENerate(name)                                      ///
        BANDwidth(numlist min=1)                            ///
        [ CENTER(varlist numeric)                           ///
          WVAR(varname numeric)                             ///
          KERNEL(string)                                    ///
          FOLDS(integer 5)                                  ///
          ORDER(integer 1) ]

    /*
        The first variable is the outcome. The remaining variables are
        the training covariates.
        e.g. varlist = "y x1 x2" then 
            depvar = "y"
            xvars  = "x1 x2"
    */
    gettoken depvar xvars : varlist

    /*
        Define the estimation sample and remove incomplete training rows.
        markout: changes `touse' to 0 when an obs has missing values in the listed variables
    */
    marksample touse
    markout `touse' `depvar' `xvars'

    if "`wvar'" != "" {
        markout `touse' `wvar'
    }

    /*
        Validate polynomial order.
    */
    if `order' < 0 {
        display as error "order() must be a nonnegative integer"
        exit 198
    }

    /*
        Validate all candidate bandwidths.
    */
    foreach h of numlist `bandwidth' {
        if `h' <= 0 {
            display as error "all values in bandwidth() must be strictly positive"
            exit 198
        }
    }

    /*
        Cross-validation requires at least two folds whenever more than
        one bandwidth is supplied.
    */
    local n_bandwidths : word count `bandwidth'

    if `n_bandwidths' > 1 & `folds' < 2 {
        display as error "folds() must be at least 2 when bandwidth() contains multiple values"
        exit 198
    }

    /*
        Set and validate the kernel.
    */
    if "`kernel'" == "" {
        local kernel "gaussian"
    }

    local kernel = lower("`kernel'")

    if !inlist("`kernel'", "gaussian", "uniform", "triangular") {
        display as error "kernel() must be gaussian, uniform, or triangular"
        exit 198
    }

    /*
        Training, evaluation, and centering points must have the same
        covariate dimension.
    */
    local p_train : word count `xvars'
    local p_at    : word count `at'

    if `p_train' != `p_at' {
        display as error "at() must contain the same number of variables as the training covariates"
        exit 198
    }

    /*
        If center() is omitted, evaluate and center at the same points.
    */
    if "`center'" == "" {
        local center "`at'"
    }

    local p_center : word count `center'

    if `p_train' != `p_center' {
        display as error ///
            "center() must contain the same number of variables as the training covariates"
        exit 198
    }

    /*
        The output variable must not already exist.
    */
    confirm new variable `generate'
    quietly generate double `generate' = .

    /*
        If no weights are supplied, use unit weights.
    */
    if "`wvar'" == "" {
        tempvar unitweight
        quietly generate double `unitweight' = 1
        local wvar "`unitweight'"
    }

    /*
        Numerical weights must be nonnegative and cannot all be zero.
    */
    quietly count if `touse' & `wvar' < 0

    if r(N) > 0 {
        display as error "wvar() must be nonnegative"
        exit 198
    }

    quietly summarize `wvar' if `touse', meanonly

    if r(sum) <= 0 {
        display as error "wvar() must contain at least one strictly positive weight"
        exit 198
    }

    /*
        Check the estimation-sample size.
    */
    quietly count if `touse'
    local n_train = r(N)

    if `n_train' == 0 {
        display as error "no complete observations in the estimation sample"
        exit 2000
    }

    if `n_bandwidths' > 1 & `folds' > `n_train' {
        display as error ///
            "folds() cannot exceed the number of estimation observations"
        exit 198
    }

    /*
        Create a temporary scalar for the selected bandwidth
        A polynomial of total degree order in p variables contains
        choose(p + order, order) terms. Mata constructs this basis.
    */
    tempname selected_bandwidth

    scalar `selected_bandwidth' = .

    mata: _rdcomono_localpoly_driver(                       ///
        "`depvar'",                                         ///
        "`xvars'",                                          ///
        "`at'",                                             ///
        "`center'",                                         ///
        "`wvar'",                                           ///
        "`touse'",                                          ///
        "`generate'",                                       ///
        "`bandwidth'",                                      ///
        `folds',                                            ///
        `order',                                            ///
        "`kernel'",                                         ///
        "`selected_bandwidth'"                              ///
    )

    /*
        Return results needed by the later RDD estimator.
    */
    return scalar bandwidth = scalar(`selected_bandwidth')
    return scalar folds = `folds'
    return scalar order = `order'
    return scalar n = `n_train'

    return local kernel "`kernel'"
    return local predictvar "`generate'"
    return local candidate_bandwidths "`bandwidth'"
end


mata:
mata set matastrict on


/**********************************************************************
    Construct exponent vectors for a total-degree polynomial basis
**********************************************************************/

real matrix _rdcomono_exponent_matrix(
    real scalar p,
    real scalar polynomial_order
)
{
    real scalar base
    real scalar total_rows
    real scalar code
    real scalar j
    real scalar remainder

    real matrix exponents
    real colvector keep
    real colvector total_degree
    real colvector sorting_index

    /*
        Enumerate every p-vector whose entries lie in
        {0, ..., polynomial_order}. We then retain only vectors whose
        entries sum to at most polynomial_order.

        Example: p = 2 and order = 2 gives

            (0,0), (1,0), (0,1), (2,0), (1,1), (0,2)

        up to column ordering.
    */
    base = polynomial_order + 1
    total_rows = base^p

    exponents = J(total_rows, p, 0)

    for (code = 0; code < total_rows; code++) {
        remainder = code

        for (j = 1; j <= p; j++) {
            exponents[code + 1, j] = mod(remainder, base)
            remainder = floor(remainder / base)
        }
    }

    keep = rowsum(exponents) :<= polynomial_order
    exponents = select(exponents, keep)

    /*
        Put lower-degree terms first. This ensures that the intercept is
        the first basis term and makes the matrix easier to inspect.
    */
    total_degree = rowsum(exponents)
    sorting_index = order((total_degree, exponents), 1..(p + 1))
    exponents = exponents[sorting_index, .]

    return(exponents)
}


/**********************************************************************
    Evaluate a multivariate polynomial basis
**********************************************************************/

real matrix _rdcomono_polynomial_basis(
    real matrix X,
    real matrix exponents
)
{
    real scalar term
    real scalar variable
    real scalar exponent_value

    real matrix basis

    basis = J(rows(X), rows(exponents), 1)

    /*
        For exponent vector a = (a1,...,ap), construct

            X1^a1 * ... * Xp^ap.

        The exponent matrix contains every monomial with total degree no
        greater than order().
    */
    for (term = 1; term <= rows(exponents); term++) {
        for (variable = 1; variable <= cols(exponents); variable++) {
            exponent_value = exponents[term, variable]

            if (exponent_value != 0) {
                basis[, term] =
                    basis[, term] :*
                    (X[, variable] :^ exponent_value)
            }
        }
    }

    return(basis)
}


/**********************************************************************
    Kernel-weight calculation
**********************************************************************/

real colvector _rdcomono_kernel_weights(
    real colvector distance,
    real scalar h,
    string scalar kernel
)
{
    real colvector kernel_weight

    if (kernel == "gaussian") {
        /*
            The common Gaussian normalizing constant is omitted because
            it cancels from weighted least squares at a fixed center.
        */
        kernel_weight =
            exp(-0.5 :* (distance :/ h) :^ 2)
    }
    else if (kernel == "uniform") {
        kernel_weight =
            distance :<= h
    }
    else {
        /*
            The ado layer has already verified that the remaining case
            is the triangular kernel.
        */
        kernel_weight =
            (distance :<= h) :*
            (1 :- distance :/ h)
    }

    return(kernel_weight)
}


/**********************************************************************
    Core multivariate local-polynomial predictor
    Identical (X_at, X_center) pairs are evaluated only once.
**********************************************************************/

real colvector _rdcomono_poly_predict(
    real colvector y,
    real matrix X,
    real colvector sampling_weight,
    real matrix X_at,
    real matrix X_center,
    real scalar h,
    real matrix exponents,
    string scalar kernel
)
{
    real scalar i
    real scalar n_eval
    real scalar n_unique
    real scalar p

    real matrix training_basis
    real matrix evaluation_basis

    real matrix pair_matrix
    real matrix pair_sorted

    real matrix X_at_unique
    real matrix X_center_unique

    real colvector sorting_index
    real colvector group_sorted
    real colvector unique_row_position

    real colvector predictions
    real colvector unique_predictions

    real colvector distance
    real colvector kernel_weight
    real colvector total_weight

    real matrix weighted_crossproduct
    real colvector weighted_outcome_crossproduct
    real colvector beta


    /** 1. Identify unique (X_at, X_center) pairs **/

    n_eval = rows(X_at)
    p = cols(X_at)

    pair_matrix =
        (
            X_at,
            X_center
        )

    sorting_index =
        order(
            pair_matrix,
            1..cols(pair_matrix)
        )

    pair_sorted =
        pair_matrix[sorting_index, .]


    /*
        group_sorted tells us which unique pair each sorted row belongs to.

        unique_row_position stores the first row of every unique pair.
    */

    group_sorted =
        J(n_eval, 1, 1)

    unique_row_position =
        J(n_eval, 1, .)

    n_unique = 1

    unique_row_position[1] = 1


    for (i = 2; i <= n_eval; i++) {

        if (
            any(
                pair_sorted[i, .]
                :!=
                pair_sorted[i - 1, .]
            )
        ) {

            n_unique = n_unique + 1

            unique_row_position[n_unique] = i
        }

        group_sorted[i] = n_unique
    }


    unique_row_position =
        unique_row_position[1..n_unique]


    X_at_unique =
        pair_sorted[
            unique_row_position,
            1..p
        ]


    X_center_unique =
        pair_sorted[
            unique_row_position,
            (p + 1)..(2 * p)
        ]


    /** 2. Construct polynomial bases

        Training basis is unchanged.

        Evaluation basis now only needs to be constructed at unique
        evaluation points. **/

    training_basis =
        _rdcomono_polynomial_basis(
            X,
            exponents
        )


    evaluation_basis =
        _rdcomono_polynomial_basis(
            X_at_unique,
            exponents
        )


    /** 3. Run local polynomial regression only for unique pairs **/

    unique_predictions =
        J(n_unique, 1, .)


    for (i = 1; i <= n_unique; i++) {

        /*
            Kernel weights depend on X_center.
        */

        distance =
            sqrt(
                rowsum(
                    (
                        X :-
                        X_center_unique[i, .]
                    ) :^ 2
                )
            )


        kernel_weight =
            _rdcomono_kernel_weights(
                distance,
                h,
                kernel
            )


        total_weight =
            sampling_weight :*
            kernel_weight


        /*
            No observations receive positive weight.
        */

        if (sum(total_weight) <= 0) {
            continue
        }


        /*
            Weighted local-polynomial regression.
        */

        weighted_crossproduct =
            cross(
                training_basis,
                total_weight,
                training_basis
            )


        weighted_outcome_crossproduct =
            cross(
                training_basis,
                total_weight,
                y
            )


        beta =
            qrsolve(
                weighted_crossproduct,
                weighted_outcome_crossproduct
            )


        /*
            Evaluate at X_at.
        */

        unique_predictions[i] =
            evaluation_basis[i, .] *
            beta
    }


    /** 4. Expand predictions back to the original observations **/

    predictions =
        J(n_eval, 1, .)


    predictions[sorting_index] =
        unique_predictions[group_sorted]


    return(predictions)
}


/**********************************************************************
K-fold bandwidth selection
**********************************************************************/

real scalar _rdcomono_select_bandwidth(
    real colvector y,
    real matrix X,
    real colvector sampling_weight,
    real rowvector bandwidths,
    real scalar num_folds,
    real matrix exponents,
    string scalar kernel
)
{
    real scalar n
    real scalar fold_size
    real scalar bandwidth_index
    real scalar fold

    real colvector permutation
    real colvector validation_index
    real colvector training_index
    real colvector fold_marker

    real colvector cv_criterion
    real colvector fold_prediction

    real colvector valid_bandwidth_index
    real colvector minimum_index

    n = rows(y)
    fold_size = floor(n / num_folds)

    /*
        Suffle observations, then form consecutive equal-size folds.
        Any remainder observations are never validation observations
        but remain in complementray training samples.
    */
    permutation = order(runiform(n, 1), 1)

    cv_criterion = J(cols(bandwidths), 1, .)

    for (
        bandwidth_index = 1;
        bandwidth_index <= cols(bandwidths);
        bandwidth_index++
    ) {
        cv_criterion[bandwidth_index] = 0

        for (fold = 1; fold <= num_folds; fold++) {

            validation_index =
                permutation[
                    ((fold - 1) * fold_size + 1)
                    ..
                    (fold * fold_size)
                ]

            fold_marker = J(n, 1, 0)
            fold_marker[validation_index] = J(rows(validation_index), 1, 1)

            training_index =
                selectindex(
                    fold_marker :== 0
                )

            fold_prediction =
                _rdcomono_poly_predict(
                    y[training_index],
                    X[training_index, .],
                    sampling_weight[training_index],
                    X[validation_index, .],
                    X[validation_index, .],
                    bandwidths[bandwidth_index],
                    exponents,
                    kernel
                )

            /*
                Missing predictions indicate that the candidate bandwidth
                cannot support every validation fit. Mark that candidate
                as invalid.
            */
            if (any(fold_prediction :>= .)) {
                cv_criterion[bandwidth_index] = .
                break
            }

            cv_criterion[bandwidth_index] =
                cv_criterion[bandwidth_index] +
                sum(
                    sampling_weight[validation_index] :*
                    (
                        y[validation_index] :-
                        fold_prediction
                    ) :^ 2
                )
        }

        if (cv_criterion[bandwidth_index] < .) {
            cv_criterion[bandwidth_index] =
                cv_criterion[bandwidth_index] / n
        }
    }

    valid_bandwidth_index =
        selectindex(
            cv_criterion :< .
        )

    /*
        No candidate produced a valid fit in every fold.
    */
    if (rows(valid_bandwidth_index) == 0) {
        return(.)
    }

    minimum_index =
        selectindex(
            cv_criterion[valid_bandwidth_index] :==
            min(cv_criterion[valid_bandwidth_index])
        )

    return(
        bandwidths[
            valid_bandwidth_index[
                minimum_index[1]
            ]
        ]
    )
}


/**********************************************************************
    Main driver called by the ado program
**********************************************************************/

void _rdcomono_localpoly_driver(
    string scalar y_name,
    string scalar x_names,
    string scalar at_names,
    string scalar center_names,
    string scalar weight_name,
    string scalar touse_name,
    string scalar prediction_name,
    string scalar bandwidth_string,
    real scalar num_folds,
    real scalar polynomial_order,
    string scalar kernel,
    string scalar selected_bandwidth_scalar
)
{
    real colvector y
    real matrix X
    real colvector sampling_weight

    real matrix X_at
    real matrix X_center
    real matrix exponents

    real colvector evaluation_complete
    real colvector evaluation_index

    real rowvector bandwidths
    real scalar selected_bandwidth

    real colvector predictions

    /*
        Read the complete estimation sample from Stata.
    */
    y =
        st_data(
            .,
            y_name,
            touse_name
        )

    X =
        st_data(
            .,
            tokens(x_names),
            touse_name
        )

    sampling_weight =
        st_data(
            .,
            weight_name,
            touse_name
        )

    bandwidths =
        strtoreal(
            tokens(bandwidth_string)
        )

    /*
        Generate every monomial exponent vector with total degree less
        than or equal to order().
    */
    exponents =
        _rdcomono_exponent_matrix(
            cols(X),
            polynomial_order
        )

    /*
        One candidate means fixed bandwidth. Multiple candidates invoke
        K-fold cross-validation.
    */
    if (cols(bandwidths) == 1) {
        selected_bandwidth =
            bandwidths[1]
    }
    else {
        selected_bandwidth =
            _rdcomono_select_bandwidth(
                y,
                X,
                sampling_weight,
                bandwidths,
                num_folds,
                exponents,
                kernel
            )
    }

    /*
        If no candidate bandwidth can produce valid CV predictions,
        leave the returned scalar missing and do not calculate the final
        fit.
    */
    if (selected_bandwidth >= .) {
        st_numscalar(
            selected_bandwidth_scalar,
            .
        )
        return
    }

    /*
        Read evaluation and center variables for all observations.
    */
    X_at =
        st_data(
            .,
            tokens(at_names)
        )

    X_center =
        st_data(
            .,
            tokens(center_names)
        )

    evaluation_complete =
        (rowmissing(X_at) :== 0) :&
        (rowmissing(X_center) :== 0)

    evaluation_index =
        selectindex(
            evaluation_complete
        )

    /*
        Nothing needs to be written if every evaluation row is incomplete.
    */
    if (rows(evaluation_index) == 0) {
        st_numscalar(
            selected_bandwidth_scalar,
            selected_bandwidth
        )
        return
    }

    X_at =
        X_at[evaluation_index, .]

    X_center =
        X_center[evaluation_index, .]

    predictions =
        _rdcomono_poly_predict(
            y,
            X,
            sampling_weight,
            X_at,
            X_center,
            selected_bandwidth,
            exponents,
            kernel
        )

    st_store(
        evaluation_index,
        prediction_name,
        predictions
    )

    st_numscalar(
        selected_bandwidth_scalar,
        selected_bandwidth
    )
}

end

# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

discretizeCpp <- function(a, b, w, mask, ij, Dx, Dy) {
    .Call(`_MultiFit_discretizeCpp`, a, b, w, mask, ij, Dx, Dy)
}

single_Fisher_test <- function(t, correct, ret_all_probs) {
    .Call(`_MultiFit_single_Fisher_test`, t, correct, ret_all_probs)
}


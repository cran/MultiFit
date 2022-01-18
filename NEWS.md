# MultiFit 0.1.1

* Added a `NEWS.md` file to track changes to the package.
Fixed issues with std::pow usage, compilation warnings regarding comparison between unsigned and signed int, shortened vigentter execution time by setting "eval=F" for longest running chunk.

# MultiFit 0.1.2
* Corrected saving of temporary files by using "file.path()".

# MultiFit 1.0.0
* New method to select tests.
* Moved to STRICT_R_HEADERS

# MultiFit 1.0.1
* Updates to vignette.

# MultiFIT 1.1.0
* Adding optional stopping rule.
* Removing legacy code: univariate permutation nulls. Modified Holm correction.

# MultiFIT 1.1.1
* Conforming to Rcpp 1.0.8
* Change DOUBLE_XMIN to DBL_MIN
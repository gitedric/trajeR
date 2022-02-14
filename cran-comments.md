## Resubmission

This is a resubmission. I believe I have fixed all the links apart from one mentioned below.

## R CMD check results

> checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-Werror=format-security’ ‘-Wp,-D_FORTIFY_SOURCE=2’
    ‘-Wp,-D_GLIBCXX_ASSERTIONS’

This is a local configuration issue.

## Test environments (via Github Actions)

* windows-latest (release)
* macOS-latest (release)
* ubuntu-20.04 (release)


## R CMD check results (via Github Actions)

There were no ERRORs or WARNINGs.

There are 1 NOTE, only for Windows:

> checking installed package size ... NOTE
    installed size is  6.2Mb
    sub-directories of 1Mb or more:
      libs   5.4Mb

Installed package size exceeding 5 MB is mainly caused by use of the packages Rcpp, RcppArmadillo and RcppEigen.


## R CMD check results (peronnal computer fedora 35)

> checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-Werror=format-security’ ‘-Wp,-D_FORTIFY_SOURCE=2’
    ‘-Wp,-D_GLIBCXX_ASSERTIONS’


This is a local configuration issue


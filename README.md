# SIMULATION OF GAUSSIAN RANDOM FIELDS ON CLOSED SURFACES

This repository contains the companion code for the paper **Non-stationary Gaussian random fields on hypersurfaces: Sampling and strong error analysis"** by E. Jansson, A. Lang and M. Pereira.

It is composed of two folders:

* The folder `package` contains the C++ and R source code for the R package `SFEMsim`, created by the authors. This package offers routines to simulate Gaussian random fields on meshed 1D or 2D surfaces using the Galerkin--Chebyshev approach presented in the paper. This package can be built by compiling these source files using Rstudio and Rcpp. An alternative is to [download directly the following archive folder containing the package](https://cloud.minesparis.psl.eu/index.php/s/1tkyLuYw7N5nKr2), and installing it using Rstudio or the R command `install.packages`. Note that the `SFEsim` packages relies on the R packages `Rcpp`, `RcppEigen`, `Matrix`, `fields` and `parallel`.

* The folder `scripts` contains  the R script used ro generate the simulation figures that appear in the paper.




% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ke.R
\name{r_array_to_d_array}
\alias{r_array_to_d_array}
\title{Calculate dipole-dipole interaction tensors from internuclear vectors}
\usage{
r_array_to_d_array(r_array, gradient = FALSE)
}
\arguments{
\item{r_array}{3D array (pairs, models, xyz) with internuclear vectors}

\item{gradient}{a logical value indicating whether to calculate the derivative}
}
\value{
3D array (pairs, models, tensor elements) with interaction tensors

The optional derivative is contained in the \code{"gradient"} attribute. It is a 4D array (pairs, models, tensor elements, xyz).
}
\description{
Calculate dipole-dipole interaction tensors from internuclear vectors
}
\details{
Given an internuclear vector with components (\code{x}, \code{y}, \code{z}) and distance (\code{r}), the vector form (\code{d}) of the dipole-dipole interaction tensor is defined:

\eqn{d = 1/r^5 [ z^2 - 1/2 (x^2+y^2), \sqrt{3}/2 (x^2-y^2), \sqrt{3} xz, \sqrt{3} yz, \sqrt{3} xy) ]}

Note that the coefficients of the last four elements are different than what was given in \href{https://doi.org/10.1007/s10858-019-00288-8}{Smith 2020}. Those published coefficients do not account for the 3/2 factor in the cartesian tensor.
}

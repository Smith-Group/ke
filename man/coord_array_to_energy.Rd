% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ke.R
\name{coord_array_to_energy}
\alias{coord_array_to_energy}
\title{Calculate restraint energy from atomic coordinates}
\usage{
coord_array_to_energy(
  coord_array,
  atom_pairs,
  grouping_list,
  g0,
  k,
  gradient = FALSE
)
}
\arguments{
\item{coord_array}{3D array (atoms, xyz, models) with atomic coordinates}

\item{atom_pairs}{matrix with each row having the names or indices of an atom pair (first dimension in \code{coord_array})}

\item{grouping_list}{list of lists of integer vectors giving groupings of models to average interaction tensors}

\item{g0}{target group norm squared values}

\item{k}{force constant}

\item{gradient}{a logical value indicating whether to calculate the derivative}
}
\value{
total restraint energy calculated using \eqn{k*(g-g0)^2}

The optional derivative is contained in the \code{"gradient"} attribute. It is a 3D array (atoms, xyz, models).
}
\description{
Calculate restraint energy from atomic coordinates
}

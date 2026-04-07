#' Internal physical constants
#'
#' Lookup table for gyromagnetic ratios and short frequency labels used by the
#' spectral-density term helpers. The physical constants satisfy
#' \deqn{\frac{\mu_0}{4\pi} = 10^{-7}}
#' in SI units, and `hbar` is the reduced Planck constant
#' \eqn{\hbar = h/(2\pi)}.
#'
#' @return A list with named numeric vector `gamma` in rad s^-1 T^-1, named
#'   character vector `label`, scalar `mu0_over_4pi` in T m A^-1, and scalar
#'   `hbar` in J s
#'
#' @noRd
.physical_constants <- list(
	gamma = c(
		"1H" = 267.52218744e6,
		"13C" = 67.28284000e6,
		"15N" = -27.11600000e6,
		"19F" = 251.66200000e6,
		"31P" = 108.29100000e6,
		"2H" = 41.06600000e6
	),
	label = c(
		"1H" = "H",
		"13C" = "C",
		"15N" = "N",
		"19F" = "F",
		"31P" = "P",
		"2H" = "D"
	),
	mu0_over_4pi = 1e-7,
	hbar = 1.054571817e-34
)

#' Get the gyromagnetic ratio of a nucleus
#'
#' The gyromagnetic ratio \eqn{\gamma_X} relates the Larmor angular frequency
#' of nucleus \eqn{X} to the static magnetic field through
#' \deqn{\omega_X = |\gamma_X| B_0.}
#'
#' @param nucleus Character scalar nucleus identifier such as `"1H"` or
#'   `"15N"`
#'
#' @return Numeric scalar gyromagnetic ratio in rad s^-1 T^-1
#'
#' @noRd
.nucleus_gamma <- function(nucleus) {
	tab <- .physical_constants
	if (!nucleus %in% names(tab$gamma)) {
		stop("Unsupported nucleus: `", nucleus, "`")
	}
	tab$gamma[[nucleus]]
}

#' Get a short frequency label for a nucleus
#'
#' @param nucleus Character scalar nucleus identifier such as `"1H"` or
#'   `"15N"`
#'
#' @return Character scalar short label such as `"H"` or `"N"`
#'
#' @noRd
.nucleus_label <- function(nucleus) {
	tab <- .physical_constants
	if (!nucleus %in% names(tab$label)) {
		stop("Unsupported nucleus: `", nucleus, "`")
	}
	tab$label[[nucleus]]
}

#' Convert proton spectrometer frequency to magnetic field strength
#'
#' If the spectrometer proton frequency is \eqn{\nu_H} in MHz, this helper
#' first converts to angular frequency
#' \deqn{\omega_H = 2\pi \nu_H \times 10^6}
#' and then evaluates
#' \deqn{B_0 = \omega_H / |\gamma_H|.}
#'
#' @param proton_mhz Spectrometer proton frequency in MHz
#'
#' @return Numeric scalar magnetic field strength in Tesla
#'
#' @noRd
.proton_mhz_to_b0 <- function(proton_mhz) {
	omega_h <- proton_mhz * 1e6 * 2 * pi
	omega_h / .nucleus_gamma("1H")
}

#' Calculate the Larmor angular frequency of a nucleus
#'
#' This helper evaluates
#' \deqn{\omega_X = |\gamma_X| B_0}
#' for nucleus \eqn{X}, using `proton_mhz` to infer \eqn{B_0}.
#'
#' @param proton_mhz Spectrometer proton frequency in MHz
#' @param nucleus Character scalar nucleus identifier
#'
#' @return Numeric scalar angular frequency in rad s^-1
#'
#' @noRd
.larmor_omega <- function(proton_mhz, nucleus) {
	abs(.proton_mhz_to_b0(proton_mhz) * .nucleus_gamma(nucleus))
}

#' Calculate the squared dipolar prefactor for two nuclei
#'
#' For nuclei \eqn{I} and \eqn{S}, define the distance-independent dipolar
#' prefactor
#' \deqn{d_{IS} = \frac{\mu_0}{4\pi}\hbar \gamma_I \gamma_S.}
#' This helper returns
#' \deqn{d_{IS}^2.}
#' When an explicit internuclear distance is needed, the term-array helpers
#' combine this quantity with \eqn{r_{IS}^{-6}} separately.
#'
#' @param nucleus_i Character scalar first nucleus
#' @param nucleus_s Character scalar second nucleus
#'
#' @return Numeric scalar squared dipolar prefactor
#'
#' @noRd
.dipolar_prefactor_sq <- function(nucleus_i, nucleus_s) {
	tab <- .physical_constants
	(tab$mu0_over_4pi * tab$hbar * .nucleus_gamma(nucleus_i) * .nucleus_gamma(nucleus_s))^2
}

#' Evaluate the second Legendre polynomial
#'
#' This helper evaluates
#' \deqn{P_2(x) = \frac{3x^2 - 1}{2}.}
#'
#' @param x Numeric vector
#'
#' @return Numeric vector with `P2(x)`
#'
#' @noRd
.p2 <- function(x) {
	0.5 * (3 * x^2 - 1)
}

#' Calculate the squared CSA prefactor for a nucleus
#'
#' For observed nucleus \eqn{I} and CSA anisotropy
#' \eqn{\Delta\sigma_I}, this helper evaluates
#' \deqn{
#' c_I^2 =
#' \frac{2}{15}\omega_I^2 \Delta\sigma_I^2.
#' }
#' Here `delta_sigma_ppm` is converted to fractional units via
#' \deqn{\Delta\sigma_I = \texttt{delta\_sigma\_ppm} \times 10^{-6}.}
#'
#' @param proton_mhz Spectrometer proton frequency in MHz
#' @param nucleus Character scalar observed nucleus
#' @param delta_sigma_ppm CSA anisotropy in ppm
#'
#' @return Numeric scalar squared CSA prefactor
#'
#' @noRd
.csa_prefactor_sq <- function(proton_mhz, nucleus, delta_sigma_ppm) {
	omega <- .larmor_omega(proton_mhz, nucleus)
	delta_sigma <- delta_sigma_ppm * 1e-6
	(2 / 15) * omega^2 * delta_sigma^2
}

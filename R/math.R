#' Construct a ZYZ Euler rotation matrix
#'
#' Return the `3 x 3` rotation matrix corresponding to the active
#' `Rz(alpha) Ry(beta) Rz(gamma)` Euler-angle convention.
#'
#' This convention matches the standard parameterization used for the
#' rank-2 Wigner rotation matrix `D^(2)(alpha, beta, gamma)`.
#'
#' @param alpha first rotation angle in radians about the `z` axis
#' @param beta second rotation angle in radians about the `y` axis
#' @param gamma third rotation angle in radians about the `z` axis
#'
#' @return `3 x 3` numeric rotation matrix
#'
#' @examples
#' euler_zyz_matrix(0, 0, 0)
#' euler_zyz_matrix(pi / 2, 0, 0)
#'
#' @export
euler_zyz_matrix <- function(alpha, beta, gamma) {
	ca <- cos(alpha)
	sa <- sin(alpha)
	cb <- cos(beta)
	sb <- sin(beta)
	cg <- cos(gamma)
	sg <- sin(gamma)

	rz_alpha <- matrix(c(
		ca, -sa, 0,
		sa, ca, 0,
		0, 0, 1
	), nrow = 3, byrow = TRUE)

	ry_beta <- matrix(c(
		cb, 0, sb,
		0, 1, 0,
		-sb, 0, cb
	), nrow = 3, byrow = TRUE)

	rz_gamma <- matrix(c(
		cg, -sg, 0,
		sg, cg, 0,
		0, 0, 1
	), nrow = 3, byrow = TRUE)

	rz_alpha %*% ry_beta %*% rz_gamma
}

# Orthonormal rank-2 Cartesian basis matching r_array_to_d_array(unit=TRUE).
.rank2_real_basis <- function() {
	sr3 <- sqrt(3)

	list(
		matrix(c(
			-0.5, 0, 0,
			0, -0.5, 0,
			0, 0, 1
		), nrow = 3, byrow = TRUE),
		matrix(c(
			sr3 / 2, 0, 0,
			0, -sr3 / 2, 0,
			0, 0, 0
		), nrow = 3, byrow = TRUE),
		matrix(c(
			0, 0, sr3 / 2,
			0, 0, 0,
			sr3 / 2, 0, 0
		), nrow = 3, byrow = TRUE),
		matrix(c(
			0, 0, 0,
			0, 0, sr3 / 2,
			0, sr3 / 2, 0
		), nrow = 3, byrow = TRUE),
		matrix(c(
			0, sr3 / 2, 0,
			sr3 / 2, 0, 0,
			0, 0, 0
		), nrow = 3, byrow = TRUE)
	)
}

#' Construct the real rank-2 rotation matrix
#'
#' Return the `5 x 5` real rotation matrix acting on the rank-2 tesseral basis
#' used by [r_array_to_d_array()] with `unit = TRUE`.
#'
#' The basis ordering is
#' `(Y20, Y22c, Y21c, Y21s, Y22s)`, equivalent to the five components
#' `(d1, d2, d3, d4, d5)` returned by [r_array_to_d_array()].
#'
#' The matrix is constructed by rotating the corresponding orthonormal
#' Cartesian rank-2 basis tensors with the `3 x 3` rotation matrix from
#' [euler_zyz_matrix()]. This produces a real orthogonal representation of the
#' same rotation described by the complex Wigner matrix `D^(2)(alpha, beta, gamma)`.
#'
#' @param alpha first rotation angle in radians about the `z` axis
#' @param beta second rotation angle in radians about the `y` axis
#' @param gamma third rotation angle in radians about the `z` axis
#'
#' @return `5 x 5` numeric rotation matrix
#'
#' @examples
#' d_real_rank2(0, 0, 0)
#'
#' x <- c(1, 0, 0)
#' d_x <- drop(r_array_to_d_array(matrix(x, nrow = 1), dist = FALSE, unit = TRUE))
#' y <- drop(euler_zyz_matrix(pi / 2, 0, 0) %*% x)
#' d_y <- drop(r_array_to_d_array(matrix(y, nrow = 1), dist = FALSE, unit = TRUE))
#' all.equal(drop(d_real_rank2(pi / 2, 0, 0) %*% d_x), d_y, tolerance = 1e-12)
#'
#' @export
# Generated with: wolframscript -code 'mvals={2,1,0,-1,-2};U={{0,0,-1,0,0},{-1/Sqrt[2],0,0,0,-1/Sqrt[2]},{0,-1/Sqrt[2],0,1/Sqrt[2],0},{0,I/Sqrt[2],0,I/Sqrt[2],0},{I/Sqrt[2],0,0,0,-I/Sqrt[2]}};dmat=FullSimplify[ComplexExpand[U.Table[WignerD[{2,m,mp},alpha,beta,gamma],{m,mvals},{mp,mvals}].ConjugateTranspose[U],TargetFunctions->{Re,Im}],Assumptions->Element[{alpha,beta,gamma},Reals]];toR[s_]:=StringReplace[ToString[InputForm[s]],{"Cos["->"cos(","Sin["->"sin(","Sqrt["->"sqrt(","["->"(","]"->")"}];Print["d_real_rank2 <- function(alpha, beta, gamma) {"];Print["\tdmat <- matrix(0, 5, 5)"];Do[Print["\tdmat[",i,", ",j,"] <- ",toR[dmat[[i,j]]]],{i,5},{j,5}];Print["\tdmat"];Print["}"]'
d_real_rank2 <- function(alpha, beta, gamma) {
	dmat <- matrix(0, 5, 5)
	ca <- cos(alpha)
	sa <- sin(alpha)
	cb <- cos(beta)
	sb <- sin(beta)
	cg <- cos(gamma)
	sg <- sin(gamma)
	c2a <- cos(2 * alpha)
	s2a <- sin(2 * alpha)
	c2b <- cos(2 * beta)
	s2b <- sin(2 * beta)
	c2g <- cos(2 * gamma)
	s2g <- sin(2 * gamma)
	sr3 <- sqrt(3)
	sb_sq <- sb^2
	three_plus_c2b <- 3 + c2b
	cb_cg <- cb * cg
	cb_sg <- cb * sg
	ca_cb <- ca * cb
	sa_cb <- sa * cb
	ca_c2b <- ca * c2b
	sa_c2b <- sa * c2b
	ca_sa <- ca * sa

	dmat[1, 1] <- (1 + 3 * c2b) / 4
	dmat[1, 2] <- (sr3 * c2g * sb_sq) / 2
	dmat[1, 3] <- -(sr3 * cb_cg * sb)
	dmat[1, 4] <- sr3 * cb * sb * sg
	dmat[1, 5] <- -(sr3 * cg * sb_sq * sg)
	dmat[2, 1] <- (sr3 * c2a * sb_sq) / 2
	dmat[2, 2] <- (c2a * three_plus_c2b * c2g) / 4 - cb * s2a * s2g
	dmat[2, 3] <- (c2a * cg * s2b) / 2 - s2a * sb * sg
	dmat[2, 4] <- -(sb * (cg * s2a + c2a * cb_sg))
	dmat[2, 5] <- -(cb * c2g * s2a) - (c2a * three_plus_c2b * s2g) / 4
	dmat[3, 1] <- sr3 * ca_cb * sb
	dmat[3, 2] <- -(ca * c2g * s2b) / 2 + sa * sb * s2g
	dmat[3, 3] <- ca_c2b * cg - cb_sg * sa
	dmat[3, 4] <- -(cb_cg * sa) - ca_c2b * sg
	dmat[3, 5] <- sb * (c2g * sa + ca_cb * s2g)
	dmat[4, 1] <- sr3 * sa_cb * sb
	dmat[4, 2] <- -(sb * (cb * c2g * sa + ca * s2g))
	dmat[4, 3] <- sa_c2b * cg + ca * cb_sg
	dmat[4, 4] <- ca_cb * cg - sa_c2b * sg
	dmat[4, 5] <- sb * (-(ca * c2g) + cb * sa * s2g)
	dmat[5, 1] <- sr3 * ca_sa * sb_sq
	dmat[5, 2] <- (three_plus_c2b * c2g * s2a) / 4 + c2a * cb * s2g
	dmat[5, 3] <- sb * (cb_cg * s2a + c2a * sg)
	dmat[5, 4] <- sb * (c2a * cg - s2a * cb_sg)
	dmat[5, 5] <- c2a * cb * c2g - ca_sa * three_plus_c2b * cg * sg
	dmat
}

# Reference version
.d_real_rank2 <- function(alpha, beta, gamma) {
	r <- euler_zyz_matrix(alpha, beta, gamma)
	basis <- .rank2_real_basis()
	basis_norm_sq <- sum(basis[[1]] * basis[[1]])

	value <- matrix(0, nrow = 5, ncol = 5)

	for (j in seq_len(5)) {
		basis_rot <- r %*% basis[[j]] %*% t(r)
		for (i in seq_len(5)) {
			value[i, j] <- sum(basis[[i]] * basis_rot) / basis_norm_sq
		}
	}

	value
}

library(ke)

pdb_path <- system.file("extdata", "gb3", "2lum_subset.pdb.gz", package = "ke")
pdb2lum <- read_ensemble(pdb_path, proton_only = TRUE)

# coordinates against which derivatives will be calculated
coord_array <- pdb2lum[, , c(1, 3, 5)]

# coordinates for creating synthetic relaxation rates
coord_synthetic <- pdb2lum[, , c(2, 4, 6)]

# ensemble member interconversion eigenvalue
base_rates <- c(kens = 1 / 2e-9)

# internal permutation rates
kmethyl <- 1 / 1e-12
karo <- 1 / 100e-6

# isotropic rotational diffusion constant used for both forward and gradient checks
d_rot <- 1 / 4e-9

# proton field strength in MHz
proton_mhz <- 700

# transition rate matrix for ensemble members
base_rate_mat <- rate_mat_simple(base_rates, dimnames(coord_synthetic)[[3]])

# data structure to hold data related to kinetic ensemble calculation
ke_data <- make_ke_data(
	coord_synthetic,
	base_rate_mat,
	base_rates,
	kc = d_rot,
	proton_mhz = proton_mhz,
	mix_times = numeric()
)

# calculate an array with interatomic distances for three ensemble members
dist_array <- simplify2array(apply(coord_synthetic, 3, function(x) as.matrix(dist(t(x))), simplify = FALSE))

# calculate matrix with minimum distances across three ensemble members
min_dist_mat <- apply(dist_array, 1:2, min)

# create empty matrix with minimum distances between equivalent atom groups
equiv_min_dist_mat <- matrix(
	NA_real_,
	nrow = length(ke_data$equiv_list),
	ncol = length(ke_data$equiv_list),
	dimnames = list(names(ke_data$equiv_list), names(ke_data$equiv_list))
)

# create a named vector mapping for equivalent atom groups with just a single atom
equiv_single <- unlist(ke_data$equiv_list[sapply(ke_data$equiv_list, length) == 1])

# copy distances from equivalent groups with just a single atom
equiv_min_dist_mat[names(equiv_single), names(equiv_single)] <- min_dist_mat[equiv_single, equiv_single]

# find equivalent groups with multiple atoms
equiv_multi <- ke_data$equiv_list[sapply(ke_data$equiv_list, length) > 1]

# create pairs of all possible combinations of all groups with multi groups
equiv_multi_pairs <- expand.grid(names(ke_data$equiv_list), names(equiv_multi), stringsAsFactors = FALSE)

# compute minimum distances for multi group pairs
for (i in seq_len(nrow(equiv_multi_pairs))) {
	min_dist <- min(min_dist_mat[
		ke_data$equiv_list[[equiv_multi_pairs[i, 1]]],
		ke_data$equiv_list[[equiv_multi_pairs[i, 2]]]
	])
	equiv_min_dist_mat[equiv_multi_pairs[i, 1], equiv_multi_pairs[i, 2]] <-
		equiv_min_dist_mat[equiv_multi_pairs[i, 2], equiv_multi_pairs[i, 1]] <- min_dist
}

# select pairs of equivalent atom groups whose minimum distances are less than 5 angstroms
idx <- which(equiv_min_dist_mat < 5 & upper.tri(equiv_min_dist_mat), arr.ind = TRUE)
idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
sigma_pairs <- matrix(rownames(equiv_min_dist_mat)[idx], ncol = 2)

coord_array_aperm <- aperm(coord_array, c(2, 1, 3))
coord_synthetic_aperm <- aperm(coord_synthetic, c(2, 1, 3))
rates <- c(base_rates, kmethyl = kmethyl, karo = karo, Dx = d_rot, Dy = d_rot, Dz = d_rot)

perm_rates_by_type <- list(
	"1-1" = c(NA_real_, NA_real_),
	"1-2" = c(NA_real_, karo = karo),
	"1-3" = c(NA_real_, kmethyl = kmethyl),
	"2-2" = c(karo = karo, karo = karo),
	"2-3" = c(karo = karo, kmethyl = kmethyl),
	"3-3" = c(kmethyl = kmethyl, kmethyl = kmethyl),
	"2 Internal" = c(NA_real_, karo = karo),
	"3 Internal" = c(NA_real_, kmethyl = kmethyl)
)

spec_den_to_relax <- function(spec_den_data, type_name) {
	sigma_value <- spec_den_data$atom_pairs[, "sigma"]
	n_relax_rates <- sum(!is.na(sigma_value))
	make_spec_den_relax_data(
		atom_relax_data = list(
			atom_pairs = spec_den_data$atom_pairs,
			unit = FALSE,
			relax_data_list = list(
				sigma = list(
					value = unname(sigma_value[!is.na(sigma_value)]),
					spec_den_term_array = make_sigma_spec_den_term_array(
						n_pairs = n_relax_rates,
						proton_mhz = proton_mhz
					)
				)
			)
		),
		base_rate_mat = base_rate_mat,
		base_rates = base_rates,
		perm_rates = perm_rates_by_type[[type_name]]
	)
}

spec_den_data_list <- make_spec_den_data(ke_data, sigma_pairs, perm_internal = TRUE)

# calculate synthetic cross relaxation rates
sigma_synthetic <- coord_array_to_sigma(coord_synthetic_aperm, c(base_rates, kmethyl = kmethyl, karo = karo, kc = d_rot), spec_den_data_list, proton_mhz)

# insert synthetic cross relaxation rates into spec_den_data_list
for (i in seq_along(spec_den_data_list)) {
	spec_den_data_list[[i]]$atom_pairs <- data.frame(spec_den_data_list[[i]]$atom_pairs, sigma = NA_real_)
	spec_den_data_list[[i]]$atom_pairs[seq_along(sigma_synthetic[[i]]), "sigma"] <- sigma_synthetic[[i]]
}

spec_den_relax_data_list <- Map(spec_den_to_relax, spec_den_data_list, names(spec_den_data_list))
names(spec_den_relax_data_list) <- names(spec_den_data_list)

relax_energy <- coord_array_to_relax_energy(
	coord_array_aperm,
	rates,
	spec_den_relax_data_list,
	gradient = TRUE
)

relax_deriv_check <- deriv_check(
	coord_array_to_relax_energy,
	coord_array_aperm,
	dv = 1e-6,
	vdims = 1:3,
	gdims = 1:3,
	rates = rates,
	spec_den_relax_data_list = spec_den_relax_data_list
)

gradient <- attr(relax_energy, "gradient")
gradient_fd <- gradient - relax_deriv_check
overall_slope <- lsfit(as.vector(gradient_fd), as.vector(gradient), intercept = FALSE)$coef

message(
	"Combined max absolute finite-difference discrepancy:\n",
	max(abs(relax_deriv_check))
)
message(
	"Combined analytical vs finite-difference slope:\n",
	unname(overall_slope)
)

oldpar <- par(no.readonly = TRUE)
par(mfrow = c(1, 1))
plot(gradient_fd, gradient, xlab = "Finite Difference Gradient", ylab = "Analytical Gradient")
abline(0, 1, col = "red")

spec_den_data_ind_list <- make_spec_den_data(ke_data, sigma_pairs, perm_internal = FALSE)
spec_den_data_int_list <- make_spec_den_data(ke_data, matrix(character(), ncol = 2), perm_internal = TRUE)
names(spec_den_data_int_list) <- c("2 Internal", "3 Internal")
spec_den_data_ind_list <- c(spec_den_data_ind_list, spec_den_data_int_list)

# calculate synthetic cross relaxation rates
sigma_synthetic_ind <- coord_array_to_sigma(
	coord_synthetic_aperm,
	c(base_rates, kmethyl = kmethyl, karo = karo, kc = d_rot),
	spec_den_data_ind_list,
	proton_mhz
)

# insert synthetic cross relaxation rates
for (i in seq_along(spec_den_data_ind_list)) {
	spec_den_data_ind_list[[i]]$atom_pairs <- data.frame(spec_den_data_ind_list[[i]]$atom_pairs, sigma = NA_real_)
	spec_den_data_ind_list[[i]]$atom_pairs[seq_along(sigma_synthetic_ind[[i]]), "sigma"] <- sigma_synthetic_ind[[i]]
}

spec_den_relax_data_ind_list <- Map(spec_den_to_relax, spec_den_data_ind_list, names(spec_den_data_ind_list))
names(spec_den_relax_data_ind_list) <- names(spec_den_data_ind_list)

if (!"gradient_list" %in% ls()) {
	gradient_list <- lapply(seq_along(spec_den_relax_data_ind_list), function(i) {
		relax_energy_i <- coord_array_to_relax_energy(
			coord_array_aperm,
			rates,
			spec_den_relax_data_ind_list[i],
			gradient = TRUE
		)

		relax_deriv_check_i <- deriv_check(
			coord_array_to_relax_energy,
			coord_array_aperm,
			dv = 1e-6,
			vdims = 1:3,
			gdims = 1:3,
			rates = rates,
			spec_den_relax_data_list = spec_den_relax_data_ind_list[i]
		)

		gradient_i <- attr(relax_energy_i, "gradient")
		gradient_fd_i <- gradient_i - relax_deriv_check_i

		list(analytical = gradient_i, finite_difference = gradient_fd_i)
	})
	names(gradient_list) <- names(spec_den_relax_data_ind_list)
}

summary_df <- do.call(
	rbind,
	lapply(names(gradient_list), function(type) {
		gradient_fd_i <- gradient_list[[type]]$finite_difference
		gradient_i <- gradient_list[[type]]$analytical
		data.frame(
			type = type,
			max_abs_diff = max(abs(gradient_i - gradient_fd_i)),
			slope = unname(lsfit(as.vector(gradient_fd_i), as.vector(gradient_i), intercept = FALSE)$coef)
		)
	})
)
row.names(summary_df) <- NULL
message(
	"Per-group derivative summary:\n",
	paste(capture.output(print(summary_df)), collapse = "\n")
)

par(mfrow = c(2, 4))
for (type in names(gradient_list)) {
	gradient_fd_i <- gradient_list[[type]]$finite_difference
	gradient_i <- gradient_list[[type]]$analytical
	slope <- lsfit(gradient_fd_i, gradient_i, intercept = FALSE)$coef

	plot(
		gradient_fd_i,
		gradient_i,
		asp = 1,
		xlab = "Finite Difference Gradient",
		ylab = "Analytical Gradient",
		main = type
	)

	abline(0, slope, col = "red")
	abline(0, 1, col = "blue")

	legend(
		"topleft",
		legend = c("y = x", sprintf("y = %0.6fx", slope)),
		lwd = 1,
		col = c("blue", "red"),
		bty = "n"
	)
}

par(oldpar)

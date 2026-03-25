test_that("rank-2 rotation matches rotating xyz vectors before tensor conversion", {
	r_vectors <- rbind(
		c(1, 0, 0),
		c(0, 1, 0),
		c(0, 0, 1),
		c(1, 1, 0),
		c(1, -2, 3),
		c(-0.4, 0.7, 1.2)
	)

	angle_sets <- rbind(
		c(0, 0, 0),
		c(pi / 2, 0, 0),
		c(0, pi / 2, 0),
		c(0.3, 0.4, 0.5),
		c(-0.9, 1.1, 0.2)
	)

	for (i in seq_len(nrow(angle_sets))) {
		alpha <- angle_sets[i, 1]
		beta <- angle_sets[i, 2]
		gamma <- angle_sets[i, 3]
		rmat <- euler_zyz_matrix(alpha, beta, gamma)
		dmat <- d_real_rank2(alpha, beta, gamma)

		d_array <- r_array_to_d_array(r_vectors, dist = FALSE, unit = TRUE)
		d_rot_from_tensor <- t(dmat %*% t(d_array))
		r_rot <- t(rmat %*% t(r_vectors))
		d_rot_from_xyz <- r_array_to_d_array(r_rot, dist = FALSE, unit = TRUE)

		expect_equal(unname(d_rot_from_tensor), unname(d_rot_from_xyz), tolerance = 1e-12)
	}
})

test_that("rank-2 rotation commutes with averaging over vectors", {
	r_vectors <- rbind(
		c(1, 0, 0),
		c(0, 1, 1),
		c(1, -1, 2),
		c(-2, 0.5, 1),
		c(0.2, -0.7, 1.4)
	)

	angle_sets <- rbind(
		c(0.3, 0.4, 0.5),
		c(-1.1, 0.7, 1.3),
		c(2.0, 1.2, -0.8)
	)

	mean_tensor <- function(d_array) {
		colMeans(d_array)
	}

	for (i in seq_len(nrow(angle_sets))) {
		alpha <- angle_sets[i, 1]
		beta <- angle_sets[i, 2]
		gamma <- angle_sets[i, 3]
		rmat <- euler_zyz_matrix(alpha, beta, gamma)
		dmat <- d_real_rank2(alpha, beta, gamma)

		d_array <- r_array_to_d_array(r_vectors, dist = FALSE, unit = TRUE)
		avg_then_rotate <- drop(dmat %*% mean_tensor(d_array))

		r_rot <- t(rmat %*% t(r_vectors))
		d_rot_from_xyz <- r_array_to_d_array(r_rot, dist = FALSE, unit = TRUE)
		rotate_then_average <- mean_tensor(d_rot_from_xyz)

		d_rot_from_tensor <- t(dmat %*% t(d_array))
		rotate_tensor_then_average <- mean_tensor(d_rot_from_tensor)

		expect_equal(unname(avg_then_rotate), unname(rotate_then_average), tolerance = 1e-12)
		expect_equal(unname(avg_then_rotate), unname(rotate_tensor_then_average), tolerance = 1e-12)
	}
})

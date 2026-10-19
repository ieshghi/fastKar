test_that('distance_density matches the spline above 10 kb and is continuous at the first knot', {
<<<<<<< HEAD
    lmin = environment(fastKar::splineobj$distance_spline)$z$x[1]
    expect_equal(distance_density(exp(lmin) * (1 - 1e-9)), distance_density(exp(lmin)), tolerance = 1e-6)
    d = c(1e4, 1e5, 1e7, 2e8)
    expect_equal(distance_density(d), exp(fastKar::splineobj$distance_spline(log(d))))
=======
    lmin = environment(fastKarOptim::splineobj$distance_spline)$z$x[1]
    expect_equal(distance_density(exp(lmin) * (1 - 1e-9)), distance_density(exp(lmin)), tolerance = 1e-6)
    d = c(1e4, 1e5, 1e7, 2e8)
    expect_equal(distance_density(d), exp(fastKarOptim::splineobj$distance_spline(log(d))))
>>>>>>> optim
})

test_that('distance_density stays bounded and decreasing below 10 kb', {
    d = c(0.5, 2, 100, 1e3, 5e3, 1e4)
    dens = distance_density(d)
    expect_true(all(is.finite(dens)))
    expect_true(all(diff(dens) < 0))
    expect_lt(distance_density(2), 1e-3)
})

test_that('a walk that revisits a tiny node back to back does not blow up the simulated map', {
    walk.dt = data.table(orig.id = c(1, 2, 2, 2, 3), compartment = 'A', width = c(1e5, 2, 2, 2, 1e5), tile.id = 1:5)
    walk.dt[, end := cumsum(width)][, start := end - width + 1][, mid := (start + end) / 2]
    sim = simulate_map(walk.dt)
    expect_true(all(is.finite(sim$value)))
    expect_equal(max(sim$value), max(sim[i == j & i %in% c(1, 5)]$value))
})

test_that('diagonal (dist == 0) pixels still come from diag_spline only', {
    walk.dt = data.table(orig.id = 1:4, compartment = 'A', width = c(2, 1e3, 1e5, 1e6), tile.id = 1:4)
    walk.dt[, end := cumsum(width)][, start := end - width + 1][, mid := (start + end) / 2]
    sim = simulate_map(walk.dt)
<<<<<<< HEAD
    expect_equal(sim[dist == 0]$value, exp(fastKar::splineobj$diag_spline(log(walk.dt$width[sim[dist == 0]$i]^2))))
=======
    expect_equal(sim[dist == 0]$value, exp(fastKarOptim::splineobj$diag_spline(log(walk.dt$width[sim[dist == 0]$i]^2))))
>>>>>>> optim
})

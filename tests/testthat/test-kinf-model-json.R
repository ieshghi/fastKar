## inst/python/kinf/data/fastkar_model.json is a snapshot of the R forward model for the Python search.
## If this fails, the forward model changed: rerun scripts/export_kinf_model.R (and check the Python side).

test_that('the shipped kinf model json matches the current forward model', {
    f = system.file('python', 'kinf', 'data', 'fastkar_model.json', package = 'fastKar')
    expect_true(nzchar(f))
    m = jsonlite::read_json(f, simplifyVector = TRUE)
    chk = opt('.kinf_model_check')()
    expect_equal(m$check$dist, chk$dist)
    expect_equal(m$check$dens, chk$dens, tolerance = 1e-10)
    expect_equal(m$check$wp, chk$wp)
    expect_equal(m$check$diag, chk$diag, tolerance = 1e-10)
    expect_equal(m$interchrom_density, interchrom_density, tolerance = 1e-10)
    expect_equal(m$maxval_density, distance_density(1e8), tolerance = 1e-10)
    z = environment(splineobj$distance_spline)$z
    for (k in c('x', 'y', 'b', 'c', 'd')) expect_equal(m$distance_spline[[k]], z[[k]], tolerance = 1e-10)
})

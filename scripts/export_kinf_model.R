# export_kinf_model.R: dump the fastKar Hi-C forward-model constants to inst/python/kinf/data/fastkar_model.json
# (read by inst/python/kinf/hic_model.py). Rerun whenever the forward model (splineobj, interchrom_density,
# distance_density) changes; tests/testthat/test-kinf-model-json.R fails until you do.
# Usage (from the repo root): Rscript scripts/export_kinf_model.R
suppressMessages(devtools::load_all('.', quiet = TRUE))
spl = function(f){
  z = environment(f)$z
  list(method = z$method, n = z$n, x = z$x, y = z$y, b = z$b, c = z$c, d = z$d, class = class(f), names = names(z))
}
fd = splineobj$distance_spline
fg = splineobj$diag_spline
zd = environment(fd)$z
lmin = zd$x[1]
# diag_spline is a linear function of log(widthprod) (function(x_new){a + b*x_new}); recover a, b by evaluation
ga = unname(fg(0)); gb = unname(fg(1)) - ga
stopifnot(all(abs(unname(fg(c(3,7,12,20))) - (ga + gb*c(3,7,12,20))) < 1e-10))
chk = .kinf_model_check()
out = list(distance_spline = spl(fd), diag_spline = list(type = 'linear_in_log', a = ga, b = gb),
  lmin = lmin, slope = fd(lmin, deriv = 1),
  interchrom_density = interchrom_density, maxval = 1e8,
  maxval_density = distance_density(1e8),
  check = chk)
writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, digits = NA, na = 'null'), 'inst/python/kinf/data/fastkar_model.json')

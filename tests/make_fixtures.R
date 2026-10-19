## Builds the input fixtures for the fixture-based property tests.
## Inputs only (graphs, footprints, karyotypes); outputs are computed live by the tests.
##   Rscript tests/make_fixtures.R
suppressMessages({
    library(data.table); setDTthreads(1)
    pkgload::load_all('~/git/fastKar', quiet = TRUE)
})
out_dir = '~/git/fastKar/tests/testthat/fixtures'
nkar = 8

## karyotypes stored as plain lists so fixtures stay small and package-independent
as_plain = function(gw) list(snode.id = unname(gw$snode.id), circular = unname(gw$circular),
                             cn = if (is.null(gw$dt$cn)) NULL else gw$dt$cn)
fixture = function(name, gg, ft, kars, note) {
    list(name = name, gg = gg, ft = ft, kars = lapply(kars, as_plain), note = note)
}
samp = function(gg, n = nkar) sample.gwalks(gg, n, verbose = FALSE)

fx = list()
set.seed(1)

## toy double translocation (live_site event 1)
tiles = gr.tile(parse.gr('1:1-1e6,2:1-3e6,3:1-1e6'), 1e6)
gg = makewalk(tiles, list(c(2, 3, 4), c(1, 3, 5)))$graph
fx$toy_dtrans = fixture('toy_dtrans', gg, gg$footprint, samp(gg), 'double translocation, 5 nodes')

## three-way reciprocal translocation; first karyotype = the constructed one
gw = make_threeway_repdup(1e6)
fx$toy_threeway = fixture('toy_threeway', gw$graph, gw$graph$footprint, c(list(gw), samp(gw$graph)), 'three-way repdup')

## reciprocal duplication, cis1 geometry
gw = makerepdup(1e6, 'cis1')
fx$toy_repdup = fixture('toy_repdup', gw$graph, gw$graph$footprint, c(list(gw), samp(gw$graph)), 'r-dup cis1, revisits node B')

## BFB, high CN with foldbacks
gw = gimme_bfb(1e5, 10, cycles = 4)
fx$toy_bfb = fixture('toy_bfb', gw$graph, gw$graph$footprint, c(list(gw), samp(gw$graph)), 'BFB 10 x 100 kb, 4 cycles')

## circular walk (ecDNA) plus linear background with CN 2
tiles = gr.tile(parse.gr('1:1-2e6'), 2e5)
segs = gG(nodes = tiles)$gr[, c()]
gw = gW(grl = GRangesList(segs[3:6], segs[1:10]), circular = c(TRUE, FALSE))
gw$set(cn = c(1, 2))
gw$nodes$mark(cn = ifelse(gw$nodes$dt$node.id %in% 3:6, 3, 2))
fx$toy_circ = fixture('toy_circ', gw$graph, gw$graph$footprint, list(gw), 'circular 4-node walk + CN 2 linear background')

## tiny 2 bp node revisited back to back (spline blow-up case), negative strand step
tiles = c(parse.gr('1:1-100000'), parse.gr('1:100001-100002'), parse.gr('1:100003-200000'))
gw = makewalk(tiles, list(c(1, 2, -2, 2, 3)))
fx$toy_tiny = fixture('toy_tiny', gw$graph, gw$graph$footprint, list(gw), '2 bp node visited 3x back to back')

## real CCLE events: full graph + footprint (+100 kb), and the local graph cut to the footprint
ccle_pairs = readRDS('~/projects/karyotyping/db/ccle.pairs.rds')
pick = function(ty, skip = 0) {
    for (p in ccle_pairs$pair) {
        g = readRDS(ccle_pairs[p, complex])
        ev = g$meta$events[type == ty]
        if (nrow(ev)) { if (skip == 0) return(list(p = p, gg = g, ev = ev[1])); skip = skip - 1 }
    }
}
for (ty in c('pyrgo', 'dm', 'chromothripsis', 'tyfonas')) {
    x = pick(ty)
    ft = streduce(parse.gr(x$ev$footprint) + 1e5)
    full = x$gg
    fx[[paste0('ccle_', ty)]] = fixture(paste0('ccle_', ty), full, ft, samp(full),
                                        paste(x$p, ty, 'full graph'))
    dj = x$gg$copy$disjoin(ft); loc = loosefix(dj[dj$nodes$gr %^% ft])
    fx[[paste0('ccle_', ty, '_local')]] = fixture(paste0('ccle_', ty, '_local'), loc, ft, samp(loc),
                                                  paste(x$p, ty, 'local graph'))
}
## every fixture's karyotypes must rebuild on its saved graph
for (f in fx) for (k in f$kars) invisible(gW(graph = f$gg, snode.id = k$snode.id, circular = k$circular))
for (nm in names(fx)) saveRDS(fx[[nm]], file.path(out_dir, paste0(nm, '.rds')))
print(rbindlist(lapply(fx, function(f) data.table(name = f$name, nodes = length(f$gg$nodes), kars = length(f$kars),
    walks = length(f$kars[[1]]$snode.id), ft_mb = round(sum(width(f$ft)) / 1e6, 2), note = f$note))))

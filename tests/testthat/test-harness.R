## Smoke test for the equivalence harness itself: the oracle loads next to the package
## under test and forward_simulate agrees on every toy fixture.

test_that('oracle and package under test are separate namespaces', {
    skip_if_no_oracle()
    expect_false(identical(environment(ref('forward_simulate')), environment(opt('forward_simulate'))))
})

test_that('forward_simulate matches the oracle on toy fixtures', {
    skip_if_no_oracle()
    for (nm in fixture_names(toy_only = TRUE)) {
        fx = load_fixture(nm)
        gw = fx$gw[[1]]
        new = opt('forward_simulate')(gw, target_region = fx$ft, pix.size = 1e5, depth = 10)
        old = ref('forward_simulate')(gw, target_region = fx$ft, pix.size = 1e5, depth = 10)
        expect_gm_equiv(new, old, info = nm)
    }
})

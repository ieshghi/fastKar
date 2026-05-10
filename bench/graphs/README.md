# bench/graphs/

Drop real-data gGraph `.rds` files here. `bench_all()` (in `R/bench.R`) will
pick them up automatically and run the same A/B comparison as on synthetic
graphs.

Each `.rds` must deserialize to either a `gGraph` directly, or to an object
whose `$graph` slot is a `gGraph`.

The directory is git-tracked but its contents (other than this file) are
ignored — see `bench/.gitignore`.

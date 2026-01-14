.PHONY: test bench verify

test:
	~/codexhome/cjulia -e 'using Pkg; Pkg.test()'

bench:
	~/codexhome/cjulia scripts/bench_sliced.jl

verify:
	~/codexhome/cjulia scripts/verify_cached_vs_projection.jl

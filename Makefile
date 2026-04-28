.PHONY: test bench verify

CODEX_HOME ?= $(abspath ../..)
CJULIA ?= $(CODEX_HOME)/cjulia

test:
	$(CJULIA) -e 'using Pkg; Pkg.test()'

bench:
	$(CJULIA) scripts/bench_sliced.jl

verify:
	$(CJULIA) scripts/verify_cached_vs_projection.jl

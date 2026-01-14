.PHONY: test

test:
	~/codexhome/cjulia -e 'using Pkg; Pkg.test()'

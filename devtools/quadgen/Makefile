# Regenerate the chnk.quadgalerkin auxiliary Galerkin quadrature tables.
# See README.md.

CAP2D ?= $(HOME)/NIST/cap2Dsolver
FSRC  := $(CAP2D)/fortran
# plain `FC ?=` would lose to make's built-in FC=f77 default; a
# command-line `make FC=...` still overrides this
FC    := gfortran
FFLAGS := -std=legacy -O2 -w

TABLEDIR := ../../chunkie/+chnk/+quadgalerkin

# orders currently shipped by chnk.quadgalerkin.getlogquad_aux
ORDERS := 1 2 3 4 5 6 8 10 12 16 20

DEPS := $(FSRC)/logquad.f $(FSRC)/logquad_aux.f \
        $(FSRC)/legeexps.f $(FSRC)/prini.f

.PHONY: all verify tables clean

all: verify

gentab: gentab.f $(DEPS)
	$(FC) $(FFLAGS) -o $@ $^

# Regenerate every shipped table into a scratch directory and diff it
# against what is committed. This must report no differences: it is the
# check that this generator really is the one the tables came from.
verify: gentab
	@rm -rf regen && mkdir -p regen
	@printf '%s\n' $(ORDERS) | ./gentab | python3 mktab.py regen >/dev/null
	@fail=0; for f in regen/*.m; do \
	  if ! diff -q "$$f" "$(TABLEDIR)/$$(basename $$f)" >/dev/null 2>&1; then \
	    echo "DIFFERS  $$(basename $$f)"; fail=1; \
	  else echo "ok       $$(basename $$f)"; fi; \
	done; \
	if [ $$fail -ne 0 ]; then echo "VERIFY FAILED"; exit 1; fi; \
	echo "verify: all $(words $(ORDERS)) tables reproduce byte-identically"

# Overwrite the committed tables in place. Run `make verify` first.
tables: gentab
	@printf '%s\n' $(ORDERS) | ./gentab | python3 mktab.py $(TABLEDIR)

clean:
	rm -rf gentab gentab.dSYM regen

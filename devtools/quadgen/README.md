# quadgen — regenerating the Galerkin auxiliary quadrature tables

Generator for the `chnk.quadgalerkin` auxiliary tables
(`chunkie/+chnk/+quadgalerkin/auxlog_nnodeNNN_nauxMMM.m`). Those files
are data dumps from the `chunkmatc` Fortran distribution; this directory
holds the dump and formatting scripts so they do not have to be
reconstructed again.

Nothing here is on the MATLAB path and nothing here runs in CI — it is
developer tooling, used only when adding or changing an order.

## What the tables contain

For panel order `k` the Galerkin backend needs, on `[-1,1]`:

- `ts_aux`, `ws_aux` — `naux = 2*k` auxiliary target nodes/weights, from
  `legeexps_log_lr`. These are `logquad(k)` mirrored about 0, i.e. a
  two-panel log-GGQ set; each half is exact for degree `<= k-1`, which is
  the exactness the degree-`(k-1)` L2 projection in
  `chnk.quadgalerkin.getauxquad` relies on.
- `xs0{j}`, `ws0{j}` — for each aux target `j`, a GGQ rule for the space
  `phi(x) + psi(x)*log|x - ts_aux(j)|` with `phi,psi` of degree `<= 2k-1`,
  from `hqsuppauxquad`.

## Requirements

- `gfortran`, `python3`
- a checkout of `cap2Dsolver` for the Fortran sources
  (`fortran/{logquad,logquad_aux,legeexps,prini}.f`). Point `CAP2D` at it
  if it is not at `~/NIST/cap2Dsolver`.

## Usage

```sh
make verify              # regenerate all shipped orders, diff vs committed
make tables              # overwrite the committed tables in place
make CAP2D=/path/to/cap2Dsolver verify
```

`make verify` must report that every table reproduces byte-identically.
That is the check that this really is the generator the committed tables
came from — treat any difference as a bug in the generator, not in the
tables.

To dump a single order by hand:

```sh
make gentab && echo 8 | ./gentab | python3 mktab.py /tmp/out
```

## Adding an order

1. Confirm `hqsuppauxquad` in `cap2Dsolver/fortran/logquad_aux.f` has a
   `logquad2xN` block for it. Available upstream:
   `1 2 3 4 5 6 8 10 12 14 16 18 20 24`. There are no odd orders above 7,
   so the Galerkin backend can never cover the full range that
   `chnk.quadggq` does. `18` and `24` are available but not currently
   shipped.
2. Add it to `ORDERS` in the `Makefile`, run `make verify` (existing
   tables must still match), then `make tables`.
3. Add a `case` to `chnk.quadgalerkin.getlogquad_aux` and update its
   docstring.
4. Add it to the sweep in
   `devtools/test/chunkermat_galerkin_singlepanelTest.m`.

Orders below 2 are unreachable in practice: `@chunker/chunker.m` asserts
`k >= 2`. The `k = 1` table ships only for parity with
`chnk.quadggq.ggqself_nnode001_npoly002`, which is unreachable for the
same reason.

The shipped set is currently `1 2 3 4 5 6 8 10 12 16 20`, which is exactly
the set of orders for which `chnk.quadggq` has `hqsupp`/`hsupp` PV and
hypersingular tables. That is deliberate: `chnk.quadgalerkin.setup` falls
back to those GGQ tables for `pv`/`hs` kernels, so a Galerkin order
outside that set (14 was one) would have no fallback for the `S'` term a
CFIE needs.

## Conventions, recovered from the committed tables

Both are easy to get wrong and neither is documented in the Fortran:

- `legeexps_log_lr(itype=1, n=k, nquad=k, ...)`. With `itype=1` the
  routine returns *before* its trailing remap to `[0,1]`, so nodes come
  back on `[-1,1]` and weights sum to 2. That is the form the tables
  store.
- `hqsuppauxquad(k, inode, xs, ws, npts)` returns nodes on `[0,1]`. The
  tables store `2*x-1` and `2*w`. Node order is left exactly as returned:
  targets in the right half carry mirrored, hence *descending*, nodes.
  Do not sort.

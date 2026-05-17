function geo = build_exact_aux_geo(chnkr,fcurve,ab,naux)
%CHNK.QUADGALERKIN.BUILD_EXACT_AUX_GEO  exact boundary geometry at the
% Galerkin auxiliary target parameters, sampled directly from the user
% parametrization fcurve (no polynomial interpolation).
%
% In the Fortran chunkmatc_aux scheme the user routine chunkpnt is
% called for the boundary position at every auxiliary target ts_aux.
% chunkie normally stores r/d/d2/n only at the k disc nodes, so when
% the Galerkin path needs them at the 2k aux targets it must
% polynomial-interpolate from those k samples. For coarsely resolved
% boundaries that interpolation error compounds with the (imperfect)
% closed-form ipw projection on log-tuned aux nodes and stalls the
% interior-solve error around 1e-6 even when GGQ reaches 1e-15.
%
% This helper pre-samples exact r/d/d2/n at ts_aux per panel by
% calling fcurve at the per-panel global parameter
%   t_global = ab(1,i) + (ab(2,i)-ab(1,i))*(ts_aux+1)/2.
%
% Pass the returned struct to chunkermat via
%   opts.exact_aux_geo = chnk.quadgalerkin.build_exact_aux_geo(chnkr,fcurve,ab);
% with the (chnkr, ab) pair returned by chunkerfunc.
%
% input
%   chnkr   - chunker (only chnkr.k and chnkr.nch are used)
%   fcurve  - user parametrization handle from chunkerfunc, callable as
%             [r,d,d2] = fcurve(t)
%   ab      - 2 x nch parameter intervals (chunkerfunc's second output)
%   naux    - aux target count (default 2*chnkr.k)
%
% output
%   geo - struct with fields r, d, d2, n of shape (dim, naux, nch)

if nargin < 4 || isempty(naux)
    naux = 2*chnkr.k;
end

aux = chnk.quadgalerkin.getauxquad(chnkr.k);
ts_aux = aux.ts_aux(:);
assert(numel(ts_aux)==naux, ...
    'build_exact_aux_geo: aux node count mismatch (table naux ~= requested)');

dim = chnkr.dim;
nch = chnkr.nch;

geo.r  = zeros(dim,naux,nch);
geo.d  = zeros(dim,naux,nch);
geo.d2 = zeros(dim,naux,nch);
geo.n  = zeros(dim,naux,nch);

for i = 1:nch
    a = ab(1,i); b = ab(2,i);
    h = (b-a)/2;
    ts_global = (a+b)/2 + h*ts_aux;
    [r_i,d_i,d2_i] = fcurve(ts_global(:).');
    geo.r (:,:,i) = reshape(r_i, dim, naux);
    geo.d (:,:,i) = reshape(d_i, dim, naux) * h;       % chunkie d w.r.t. [-1,1] param
    geo.d2(:,:,i) = reshape(d2_i,dim, naux) * h^2;
    dn = sqrt(sum(geo.d(:,:,i).^2,1));
    geo.n(:,:,i) = [geo.d(2,:,i); -geo.d(1,:,i)] ./ dn;
end

end

function aux = getauxquad(k,naux)
%CHNK.QUADGALERKIN.GETAUXQUAD
%
% Build the oversampled auxiliary-node infrastructure used by the
% chunkmatc_aux Galerkin scheme: a set of naux >= k Gauss-Legendre nodes
% on [-1,1] together with the L2 projection matrix ipw that collapses
% naux target rows of a kernel block onto k discretization rows.
%
% Mathematical setup. Let {x_i, w_i^d}, i=1..k, be the discretization
% Legendre nodes/weights and {y_j, w_j^a}, j=1..naux, be the auxiliary
% nodes/weights, both on [-1,1]. With L_i the i-th cardinal Legendre
% polynomial associated with the disc nodes, the L2 projection of a
% function tabulated at the aux nodes onto the disc nodes is
%
%   f(x_i) ~ (1/w_i^d) * sum_j w_j^a * L_i(y_j) * f(y_j),
%
% which exactly recovers f when f is a polynomial of degree <= naux-1.
% Therefore
%
%   ipw(i,j) = (w_j^a / w_i^d) * L_i(y_j)
%            = (w_j^a / w_i^d) * sum_l vmatr_aux(j,l) * umatr(l,i)
%
% where umatr maps disc values -> Legendre coefficients and vmatr_aux
% maps Legendre coefficients -> values at the aux nodes.
%
% input
%   k    - order of discretization (Legendre disc nodes)
%   naux - number of auxiliary nodes (default 2*k)
%
% output
%   aux - struct with fields
%       aux.k, aux.naux
%       aux.ts_disc, aux.whts_disc - k Legendre nodes/weights on [-1,1]
%       aux.umatr, aux.vmatr        - values<->coeffs at disc nodes
%       aux.ts_aux, aux.whts_aux    - naux Legendre nodes/weights on [-1,1]
%       aux.vmatr_aux               - coeffs -> values at aux nodes (naux x k)
%       aux.ainterp_aux             - lege.matrin(k, ts_aux) (== vmatr_aux*umatr)
%       aux.ipw                     - L2 projection (k x naux)

if nargin < 2 || isempty(naux)
    naux = 2*k;
end

[ts_disc,whts_disc,umatr,vmatr] = lege.exps(k);
[ts_aux,whts_aux] = lege.exps(naux);

ainterp_aux = lege.matrin(k,ts_aux);   % naux x k, interp disc -> aux
vmatr_aux = ainterp_aux;               % equivalent: coeffs -> aux values
                                       % since ainterp_aux = (pols at aux).' * umatr

% ipw(i,j) = (whts_aux(j) / whts_disc(i)) * L_i(ts_aux(j))
%          where L_i(ts_aux(j)) = ainterp_aux(j,i)
%
% Build via outer-product form to keep numerical stability obvious.
Lij = ainterp_aux.';                   % k x naux, Lij(i,j) = L_i(ts_aux(j))
ipw = (1./whts_disc(:)) .* Lij .* (whts_aux(:).');

aux = [];
aux.k = k;
aux.naux = naux;
aux.ts_disc = ts_disc;
aux.whts_disc = whts_disc;
aux.umatr = umatr;
aux.vmatr = vmatr;
aux.ts_aux = ts_aux;
aux.whts_aux = whts_aux;
aux.vmatr_aux = vmatr_aux;
aux.ainterp_aux = ainterp_aux;
aux.ipw = ipw;

end

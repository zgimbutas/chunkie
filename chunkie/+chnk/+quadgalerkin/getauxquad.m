function aux = getauxquad(k,ts_aux,ws_aux)
%CHNK.QUADGALERKIN.GETAUXQUAD
%
% Build the L2 projection infrastructure for the chunkmatc_aux Galerkin
% scheme given an auxiliary target node set (ts_aux, ws_aux) of length
% naux on [-1,1] and a disc order k.
%
% Mathematical setup. Let {x_i, w_i^d}, i=1..k, be the discretization
% Legendre nodes/weights and {y_j, w_j^a}, j=1..naux, the auxiliary
% target nodes/weights. With L_i the i-th cardinal Legendre polynomial
% associated with the disc nodes, the L2 projection of a function
% tabulated at the aux nodes onto the disc nodes is
%
%   f(x_i) ~ (1/w_i^d) * sum_j w_j^a * L_i(y_j) * f(y_j),
%
% so
%
%   ipw(i,j) = (w_j^a / w_i^d) * L_i(y_j)
%            = (w_j^a / w_i^d) * sum_l vmatr_aux(j,l) * umatr(l,i),
%
% where umatr maps disc values -> Legendre coefficients and vmatr_aux
% maps Legendre coefficients -> values at the aux nodes.
%
% If (ts_aux, ws_aux) are omitted, the function loads the
% legeexps_log_lr nodes from the corresponding chunkmatc table via
% chnk.quadgalerkin.getlogquad_aux.

if nargin < 2
    [ts_aux,ws_aux] = chnk.quadgalerkin.getlogquad_aux(k);
end

naux = numel(ts_aux);
ts_aux = ts_aux(:);
ws_aux = ws_aux(:);

[ts_disc,whts_disc,umatr,vmatr] = lege.exps(k);
ainterp_aux = lege.matrin(k,ts_aux);   % naux x k, interp disc values -> aux values
vmatr_aux = ainterp_aux;               % coeffs -> aux values (== ainterp_aux*[disc->coeffs]^-1)

% ipw(i,j) = (ws_aux(j) / whts_disc(i)) * L_i(ts_aux(j))
% with L_i(ts_aux(j)) = ainterp_aux(j,i).
%
% This is the chunkmatc_form_ipipw weighted-L2 projection: a
% quadrature-weighted approximation to <L_i, f>_L2 / <L_i, L_i>_L2.
% For log-singular integrands as a function of target (the actual
% structure of M_aux for self-block layer-potential integrals) this
% is the right operator; substituting a polynomial least-squares fit
% (e.g., pinv(ainterp_aux)) recovers polynomials of degree <= k-1 to
% machine precision but produces a less accurate matrix because the
% true M_aux has log-singular dependence on target.
Lij = ainterp_aux.';                   % k x naux
ipw = (1./whts_disc(:)) .* Lij .* (ws_aux(:).');

aux = [];
aux.k = k;
aux.naux = naux;
aux.ts_disc = ts_disc;
aux.whts_disc = whts_disc;
aux.umatr = umatr;
aux.vmatr = vmatr;
aux.ts_aux = ts_aux;
aux.whts_aux = ws_aux;
aux.vmatr_aux = vmatr_aux;
aux.ainterp_aux = ainterp_aux;
aux.ipw = ipw;

end

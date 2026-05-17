function auxquad = setup(k,type,naux)
%CHNK.QUADGALERKIN.SETUP
%
% Assemble the auxquads struct for the chunkmatc_aux Galerkin
% (aux-node L2-projection) quadrature backend. The struct combines:
%
%   * standard per-target singular sub-rules used by the self block
%     (xs0/wts0/ainterps0) and the neighbor block (xs1/wts1/ainterp1),
%     mirroring chnk.quadggq.setup; and
%   * the auxiliary target oversampling data
%     (ts_aux/whts_aux/vmatr_aux/ipw) used by the off-diagonal
%     aux-projection variant (see chnk.quadgalerkin.smoothbuildmat).
%
% quadratures will integrate functions of the form
%
%     f_0(x) + log(x-x_j)*f_1(x) + 1/(x-x_j)*f_2(x) + 1/(x-x_j)^2*f_3(x)
%   type =    'log'                       'pv'                'hs'
%
% where the f_i are polynomials of order 2*k.
%
% input
%   k    - order of Legendre nodes
%   type - 'log', 'pv', 'hs', 'smooth', or 'removable'
%   naux - number of auxiliary target nodes (default 2*k)
%
% output: auxquad struct with the same fields as chnk.quadggq.setup plus
%         the auxiliary-projection fields from chnk.quadgalerkin.getauxquad.

if nargin < 3 || isempty(naux)
    naux = 2*k;
end

npolyfac = 2;
[xs1,wts1,xs0,wts0] = chnk.quadgalerkin.getlogquad_aux(k,npolyfac);
if strcmpi(type,'pv')
    [xs0,wts0] = chnk.quadggq.gethqsuppquad(k,1);
elseif strcmpi(type,'hs')
    [xs0,wts0] = chnk.quadggq.gethqsuppquad(k,2);
elseif strcmpi(type,'removable')
    [xs0,wts0] = chnk.quadggq.getremovablequad(k,1);
end

ainterp1 = lege.matrin(k,xs1);

ainterps0 = cell(k,1);
for j = 1:k
    ainterps0{j} = lege.matrin(k,xs0{j});
end

aux = chnk.quadgalerkin.getauxquad(k,naux);

auxquad = [];
auxquad.k = k;
auxquad.type = type;

% standard singular sub-rules
auxquad.xs1 = xs1;
auxquad.wts1 = wts1;
auxquad.xs0 = xs0;
auxquad.wts0 = wts0;
auxquad.ainterp1 = ainterp1;
auxquad.ainterps0 = ainterps0;

% auxiliary target oversampling + L2 projection
auxquad.naux = aux.naux;
auxquad.ts_aux = aux.ts_aux;
auxquad.whts_aux = aux.whts_aux;
auxquad.ts_disc = aux.ts_disc;
auxquad.whts_disc = aux.whts_disc;
auxquad.umatr = aux.umatr;
auxquad.vmatr = aux.vmatr;
auxquad.vmatr_aux = aux.vmatr_aux;
auxquad.ainterp_aux = aux.ainterp_aux;
auxquad.ipw = aux.ipw;

end

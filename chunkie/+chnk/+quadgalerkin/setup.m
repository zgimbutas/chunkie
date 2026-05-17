function auxquad = setup(k,type)
%CHNK.QUADGALERKIN.SETUP
%
% Assemble the auxquads struct for the chunkmatc_aux Galerkin
% (aux-target L2-projection) quadrature backend. The struct combines:
%
%   * The auxiliary target nodes (ts_aux/whts_aux), per-target singular
%     source sub-rules (xs0{inode}/wts0{inode}), and the corresponding
%     source-interpolation matrices ainterps0{inode}. Self-block rows
%     are assembled at the naux aux targets then projected to the k
%     disc-node rows via the L2 projection matrix ipw.
%
%   * A neighbor-panel intermediate rule (xs1/wts1/ainterp1) reused
%     from the existing GGQ infrastructure.
%
% Phase 2 supports only 'log' singularity (Laplace single/double layer
% and their analogs). 'pv', 'hs', 'smooth', 'removable' fall back to
% the GGQ tables on the self block; the aux-projection apparatus is
% still built so that the off-diagonal smoothbuildmat path remains
% usable.
%
% input
%   k    - order of Legendre disc nodes
%   type - 'log' (default), 'pv', 'hs', 'smooth', or 'removable'
%
% output
%   auxquad - struct with fields
%     .k, .naux, .type
%     .ts_disc, .whts_disc, .umatr, .vmatr
%     .ts_aux, .whts_aux, .vmatr_aux, .ainterp_aux, .ipw
%     .xs0{inode}, .wts0{inode}, .ainterps0{inode}  (Galerkin or GGQ
%                  depending on type)
%     .xs1, .wts1, .ainterp1                        (neighbor rule)

if nargin < 2 || isempty(type)
    type = 'log';
end

% standard neighbor (intermediate) rule from GGQ
npolyfac = 2;
[xs1,wts1] = chnk.quadggq.getlogquad(k,npolyfac);
ainterp1 = lege.matrin(k,xs1);

if strcmpi(type,'log') || strcmpi(type,'smooth')
    % The chunkmatc_aux auxiliary Galerkin rules were designed for
    % weakly singular kernels only (smooth + log*smooth), with ~30
    % digits of precision. PV / HS / removable variants are not in
    % scope; those fall back to GGQ rules with identity ipw so the
    % self block reproduces the GGQ block exactly.
    [ts_aux,ws_aux,xs0,wts0] = chnk.quadgalerkin.getlogquad_aux(k);
    aux = chnk.quadgalerkin.getauxquad(k,ts_aux,ws_aux);
else
    [~,~,xs0,wts0] = chnk.quadggq.getlogquad(k,npolyfac);
    if strcmpi(type,'pv')
        [xs0,wts0] = chnk.quadggq.gethqsuppquad(k,1);
    elseif strcmpi(type,'hs')
        [xs0,wts0] = chnk.quadggq.gethqsuppquad(k,2);
    elseif strcmpi(type,'removable')
        [xs0,wts0] = chnk.quadggq.getremovablequad(k,1);
    end
    [ts_disc,whts_disc] = lege.exps(k);
    aux = chnk.quadgalerkin.getauxquad(k,ts_disc,whts_disc);
end

ainterps0 = cell(numel(xs0),1);
for j = 1:numel(xs0)
    ainterps0{j} = lege.matrin(k,xs0{j});
end

auxquad = [];
auxquad.k = k;
auxquad.type = type;

% per-target singular sub-rules (used on the self block)
auxquad.xs0 = xs0;
auxquad.wts0 = wts0;
auxquad.ainterps0 = ainterps0;

% neighbor (intermediate) rule
auxquad.xs1 = xs1;
auxquad.wts1 = wts1;
auxquad.ainterp1 = ainterp1;

% aux target infrastructure + L2 projection
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

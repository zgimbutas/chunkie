function sysmat = buildmat(chnkr,kern,opdims,type,auxquads,ilist)
%CHNK.QUADGALERKIN.BUILDMAT build matrix for given kernel and chnkr
% description of boundary, using the chunkmatc_aux Galerkin
% (aux-node L2-projection) quadrature scheme. Rows of every block are
% evaluated at the naux auxiliary target nodes and collapsed to the k
% disc-node rows by the L2 projection ipw, matching the Fortran
% reference's projection-everywhere design.
%
% Off-diagonal source quadrature ports chunkmatc_aux_od faithfully:
% the reference integrates every off-diagonal block adaptively
% (eps 1e-13) with a recursion early-exit once the target is far
% relative to the subinterval radius, so separated blocks reduce to a
% fixed Gauss rule. Here, any block whose aux targets come within
% 2.4x the source-panel radius (the cap2Dsolver FMM near-field
% criterion; in particular all neighbor blocks and close-but-not-
% adjacent pairs across corners) is computed by per-target adaptive
% quadrature (chnk.adapgausswts, the cadachunk analog) via
% chnk.quadgalerkin.nearbuildmat; genuinely separated blocks use the
% oversampled smooth rule of smoothbuildmat, which agrees with the
% adaptive result to machine precision there. The self block uses the
% per-target singular tables.
%
% Input:
%   chnkr - chunker object describing boundary
%   kern  - kernel function handle, kern(srcinfo,targinfo)
%   opdims - [m,n] kernel operator dimensions; scalar kernels use [1,1]
%
% Optional input (defaults in brackets):
%   type ['log']      - singularity type: 'log', 'pv', 'hs', 'removable',
%                       or 'smooth'
%   auxquads [setup]  - struct from chnk.quadgalerkin.setup
%   ilist []          - cell array of integer arrays of panels to skip
%
% Output:
%   sysmat - dense system matrix

if nargin < 3
    error('not enough arguments in chnk.quadgalerkin.buildmat');
end

if nargin < 4 || isempty(type)
    type = 'log';
end
if nargin < 5 || isempty(auxquads)
    auxquads = chnk.quadgalerkin.setup(chnkr.k,type);
end
if nargin < 6
    ilist = [];
end

k = chnkr.k;
nch = chnkr.nch;
r = chnkr.r;
adj = chnkr.adj;
d = chnkr.d;
d2 = chnkr.d2;
n = chnkr.n;

data = [];
if chnkr.hasdata
    data = chnkr.data;
end

temp = eye(opdims(2));

xs0 = auxquads.xs0;
wts0 = auxquads.wts0;

ainterps0 = auxquads.ainterps0;
nrules = numel(ainterps0);
ainterps0kron = cell(nrules,1);
for j = 1:nrules
    ainterps0kron{j} = kron(ainterps0{j},temp);
end

% aux-target / source-oversample data for far-block aux-projection
ainterp_aux = auxquads.ainterp_aux;
ipw = auxquads.ipw;
ts_src = auxquads.ts_src;
whts_src = auxquads.whts_src;
ainterp_src = auxquads.ainterp_src;
ainterp_src_kron = kron(ainterp_src,temp);

% adaptive-quadrature precomputes (older cached auxquads structs may
% predate these fields)
if isfield(auxquads,'ct')
    ct = auxquads.ct; bw = auxquads.bw;
    tadap = auxquads.tadap; wadap = auxquads.wadap;
else
    ct = lege.exps(k);
    bw = lege.barywts(k,ct);
    [tadap,wadap] = lege.exps(2*k+1);
end

% near/far decision at the reference family's 2.4-radius boundary
% (cap2Dsolver FMM near-field criterion, chunkfmm2d0npairs)
nearf = chnk.quadgalerkin.nearflags(r,ainterp_aux,2.4);

sysmat = zeros(k*nch*opdims(1),k*nch*opdims(2));
for it = 1:nch
    imat = 1 + (it-1)*k*opdims(1);
    imatend = it*k*opdims(1);
    for js = 1:nch
        if js == it || js == adj(1,it) || js == adj(2,it)
            continue
        end
        if ~isempty(ilist) && ismember(it,ilist) && ismember(js,ilist)
            continue
        end
        jmat_o = 1 + (js-1)*k*opdims(2);
        jmatend_o = js*k*opdims(2);
        if ~nearf(it,js)
            sysmat(imat:imatend,jmat_o:jmatend_o) = ...
                chnk.quadgalerkin.smoothbuildmat(r,d,n,d2,data,it,js,...
                    kern,opdims,ainterp_aux,ipw,...
                    ts_src,whts_src,ainterp_src,ainterp_src_kron);
        else
            sysmat(imat:imatend,jmat_o:jmatend_o) = ...
                chnk.quadgalerkin.nearbuildmat(r,d,n,d2,data,it,js,...
                    kern,opdims,ainterp_aux,ipw,ct,bw,tadap,wadap);
        end
    end
end

% nbor and self
for j = 1:nch

    jmat = 1 + (j-1)*k*opdims(2);
    jmatend = j*k*opdims(2);

    ibefore = adj(1,j);
    iafter = adj(2,j);

    if ibefore > 0
        if ~isempty(ilist) && ismember(ibefore,ilist) && ismember(j,ilist)
            % skip if both chunks are in the bad list
        else
            submat = chnk.quadgalerkin.nearbuildmat(r,d,n,d2,data,ibefore,j, ...
                kern,opdims,ainterp_aux,ipw,ct,bw,tadap,wadap);
            imat = 1 + (ibefore-1)*k*opdims(1);
            imatend = ibefore*k*opdims(1);
            sysmat(imat:imatend,jmat:jmatend) = submat;
        end
    end

    if iafter > 0
        if ~isempty(ilist) && ismember(iafter,ilist) && ismember(j,ilist)
            % skip
        else
            submat = chnk.quadgalerkin.nearbuildmat(r,d,n,d2,data,iafter,j, ...
                kern,opdims,ainterp_aux,ipw,ct,bw,tadap,wadap);
            imat = 1 + (iafter-1)*k*opdims(1);
            imatend = iafter*k*opdims(1);
            sysmat(imat:imatend,jmat:jmatend) = submat;
        end
    end

    if ~isempty(ilist) && ismember(j,ilist)
        % skip self
    else
        submat = chnk.quadgalerkin.diagbuildmat(r,d,n,d2,data,j,kern,opdims,...
            xs0,wts0,ainterps0kron,ainterps0,ainterp_aux,ipw);
        imat = 1 + (j-1)*k*opdims(1);
        imatend = j*k*opdims(1);
        sysmat(imat:imatend,jmat:jmatend) = submat;
    end

end

end

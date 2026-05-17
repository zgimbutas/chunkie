function sysmat = buildmat(chnkr,kern,opdims,type,auxquads,ilist)
%CHNK.QUADGALERKIN.BUILDMAT build matrix for given kernel and chnkr
% description of boundary, using the chunkmatc_aux Galerkin
% (aux-node L2-projection) quadrature scheme for self and neighbor
% panels and smooth Gauss-Legendre for the rest.
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

xs1 = auxquads.xs1;
wts1 = auxquads.wts1;
xs0 = auxquads.xs0;
wts0 = auxquads.wts0;

ainterp1 = auxquads.ainterp1;
ainterp1kron = kron(ainterp1,temp);

ainterps0 = auxquads.ainterps0;
nrules = numel(ainterps0);
ainterps0kron = cell(nrules,1);
for j = 1:nrules
    ainterps0kron{j} = kron(ainterps0{j},temp);
end

% smooth-everywhere baseline (Gauss-Legendre at k disc nodes)
wts = chnkr.wstor;
sysmat = chnk.quadnative.buildmat(chnkr,kern,opdims,1:nch,1:nch,wts);

% overwrite nbor and self
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
                kern,opdims,xs1,wts1,ainterp1kron,ainterp1);
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
                kern,opdims,xs1,wts1,ainterp1kron,ainterp1);
            imat = 1 + (iafter-1)*k*opdims(1);
            imatend = iafter*k*opdims(1);
            sysmat(imat:imatend,jmat:jmatend) = submat;
        end
    end

    if ~isempty(ilist) && ismember(j,ilist)
        % skip self
    else
        submat = chnk.quadgalerkin.diagbuildmat(r,d,n,d2,data,j,kern,opdims,...
            xs0,wts0,ainterps0kron,ainterps0,...
            auxquads.ts_aux,auxquads.ainterp_aux,auxquads.ipw);
        imat = 1 + (j-1)*k*opdims(1);
        imatend = j*k*opdims(1);
        sysmat(imat:imatend,jmat:jmatend) = submat;
    end

end

end

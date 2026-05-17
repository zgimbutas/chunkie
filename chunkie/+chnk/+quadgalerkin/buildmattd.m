function spmat = buildmattd(chnkr,kern,opdims,type,auxquads,ilist,corrections)
%CHNK.QUADGALERKIN.BUILDMATTD sparse self/near block of the Galerkin
% (chunkmatc_aux aux-projection) backend, optionally as a correction to
% the native smooth (Gauss-Legendre) baseline.
%
% Returns a sparse matrix containing only the entries on the
% (target_panel == source_panel) diagonal block and on the two neighbor
% off-diagonal blocks. All other entries are zero. When corrections is
% true, the entries are differences (Galerkin - native smooth) suitable
% for adding on top of a pre-built smooth matrix.
%
% Mirrors chnk.quadggq.buildmattd.

if nargin < 3
    error('not enough arguments in chnk.quadgalerkin.buildmattd');
end
if nargin < 6
    ilist = [];
end
if nargin < 7
    corrections = false;
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

if nargin < 4 || isempty(type)
    type = 'log';
end
if nargin < 5 || isempty(auxquads)
    auxquads = chnk.quadgalerkin.setup(k,type);
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

if corrections
    wtss = chnkr.wts;
    wtss = repmat(wtss(:).',opdims(2),1);
    wtss = reshape(wtss,opdims(2)*k,nch);
    indd = kron(eye(k),true(opdims(1),opdims(2)));
    indd = indd(:) > 0;
else
    wtss = [];
    indd = [];
end

mmat = k*nch*opdims(1);
nmat = k*nch*opdims(2);

nnz_max = k*nch*opdims(1)*k*3*opdims(2);
nnz1 = k*opdims(1)*k*opdims(2);
v = zeros(nnz_max,1);
iind = zeros(nnz_max,1);
jind = zeros(nnz_max,1);

[jj1,ii1] = meshgrid(0:k*opdims(2)-1,0:k*opdims(1)-1);
ict = 1;

for j = 1:nch
    jmat = 1 + (j-1)*k*opdims(2);
    ibefore = adj(1,j);
    iafter = adj(2,j);

    if ibefore > 0
        if ~isempty(ilist) && ismember(ibefore,ilist) && ismember(j,ilist)
            % skip
        else
            submat = chnk.quadgalerkin.nearbuildmat(r,d,n,d2,data,ibefore,j, ...
                kern,opdims,xs1,wts1,ainterp1kron,ainterp1,corrections,wtss);
            imat = 1 + (ibefore-1)*k*opdims(1);
            induse = ict:ict+nnz1-1;
            iind(induse) = ii1(:)+imat;
            jind(induse) = jj1(:)+jmat;
            v(induse) = submat(:);
            ict = ict + nnz1;
        end
    end

    if iafter > 0
        if ~isempty(ilist) && ismember(iafter,ilist) && ismember(j,ilist)
            % skip
        else
            submat = chnk.quadgalerkin.nearbuildmat(r,d,n,d2,data,iafter,j, ...
                kern,opdims,xs1,wts1,ainterp1kron,ainterp1,corrections,wtss);
            imat = 1 + (iafter-1)*k*opdims(1);
            induse = ict:ict+nnz1-1;
            iind(induse) = ii1(:)+imat;
            jind(induse) = jj1(:)+jmat;
            v(induse) = submat(:);
            ict = ict + nnz1;
        end
    end

    if ~isempty(ilist) && ismember(j,ilist)
        % skip
    else
        if isfield(auxquads,'exact_aux_geo') && ~isempty(auxquads.exact_aux_geo)
            eag = auxquads.exact_aux_geo;
        else
            eag = [];
        end
        submat = chnk.quadgalerkin.diagbuildmat(r,d,n,d2,data,j,kern,opdims,...
            xs0,wts0,ainterps0kron,ainterps0,...
            auxquads.ts_aux,auxquads.ainterp_aux,auxquads.ipw,...
            corrections,wtss,indd,eag);
        imat = 1 + (j-1)*k*opdims(1);
        induse = ict:ict+nnz1-1;
        iind(induse) = ii1(:)+imat;
        jind(induse) = jj1(:)+jmat;
        v(induse) = submat(:);
        ict = ict + nnz1;
    end
end

nz = ict-1;
spmat = sparse(iind(1:nz),jind(1:nz),v(1:nz),mmat,nmat);

end

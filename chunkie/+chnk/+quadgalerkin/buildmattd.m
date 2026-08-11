function spmat = buildmattd(chnkr,kern,opdims,type,auxquads,ilist,corrections)
%CHNK.QUADGALERKIN.BUILDMATTD sparse self/near block of the Galerkin
% (chunkmatc_aux aux-projection) backend, optionally as a correction to
% the native smooth (Gauss-Legendre) baseline.
%
% Returns a sparse matrix containing the (target_panel == source_panel)
% diagonal blocks, the two neighbor off-diagonal blocks, and every
% close-but-not-adjacent block within the 2.4-radius near boundary
% (chnk.quadgalerkin.nearflags -- the cap2Dsolver FMM near-field
% criterion), all computed with the adaptive aux-projection quadrature.
% All other entries are zero. When corrections is true, the entries are
% differences (Galerkin - native smooth) suitable for adding on top of
% a pre-built smooth matrix: this is the FMM-compatible decomposition
% (smooth rule applied by an FMM in the far field, sparse Galerkin
% corrections for all pairs within 2.4 radii, exactly as in the
% cap2Dsolver FMM path).
%
% Mirrors chnk.quadggq.buildmattd, with the extra near-pair blocks.

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

xs0 = auxquads.xs0;
wts0 = auxquads.wts0;
ainterps0 = auxquads.ainterps0;
nrules = numel(ainterps0);
ainterps0kron = cell(nrules,1);
for j = 1:nrules
    ainterps0kron{j} = kron(ainterps0{j},temp);
end

% adaptive-quadrature precomputes for the near blocks
if isfield(auxquads,'ct')
    ct = auxquads.ct; bw = auxquads.bw;
    tadap = auxquads.tadap; wadap = auxquads.wadap;
else
    ct = lege.exps(k);
    bw = lege.barywts(k,ct);
    [tadap,wadap] = lege.exps(2*k+1);
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

% close-but-not-adjacent pairs inside the 2.4-radius near boundary
nearf = chnk.quadgalerkin.nearflags(r,auxquads.ainterp_aux,2.4);
for j = 1:nch
    if adj(1,j) > 0, nearf(adj(1,j),j) = false; end
    if adj(2,j) > 0, nearf(adj(2,j),j) = false; end
end
npair = nnz(nearf);

nnz_max = k*nch*opdims(1)*k*3*opdims(2) + k*opdims(1)*k*opdims(2)*npair;
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
                kern,opdims,auxquads.ainterp_aux,auxquads.ipw,...
                ct,bw,tadap,wadap,corrections,wtss);
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
                kern,opdims,auxquads.ainterp_aux,auxquads.ipw,...
                ct,bw,tadap,wadap,corrections,wtss);
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
        submat = chnk.quadgalerkin.diagbuildmat(r,d,n,d2,data,j,kern,opdims,...
            xs0,wts0,ainterps0kron,ainterps0,...
            auxquads.ainterp_aux,auxquads.ipw,corrections,wtss,indd);
        imat = 1 + (j-1)*k*opdims(1);
        induse = ict:ict+nnz1-1;
        iind(induse) = ii1(:)+imat;
        jind(induse) = jj1(:)+jmat;
        v(induse) = submat(:);
        ict = ict + nnz1;
    end
end

% near-but-not-adjacent pairs (2.4-radius criterion)
[ilistn,jlistn] = find(nearf);
for ip = 1:npair
    it = ilistn(ip); js = jlistn(ip);
    if ~isempty(ilist) && ismember(it,ilist) && ismember(js,ilist)
        continue
    end
    submat = chnk.quadgalerkin.nearbuildmat(r,d,n,d2,data,it,js, ...
        kern,opdims,auxquads.ainterp_aux,auxquads.ipw,...
        ct,bw,tadap,wadap,corrections,wtss);
    imat = 1 + (it-1)*k*opdims(1);
    jmat = 1 + (js-1)*k*opdims(2);
    induse = ict:ict+nnz1-1;
    iind(induse) = ii1(:)+imat;
    jind(induse) = jj1(:)+jmat;
    v(induse) = submat(:);
    ict = ict + nnz1;
end

nz = ict-1;
spmat = sparse(iind(1:nz),jind(1:nz),v(1:nz),mmat,nmat);

end

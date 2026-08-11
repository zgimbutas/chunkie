function submat = nearbuildmat(r,d,n,d2,data,i,j,fkern,opdims,...
                               xs1,whts1,ainterp1kron,ainterp1,...
                               ainterp_aux,ipw,corrections,wtss)
%CHNK.QUADGALERKIN.NEARBUILDMAT
%
% Assemble a neighbor-panel block for the Galerkin (aux-projection)
% backend. Rows are evaluated at the naux auxiliary target nodes on
% panel i (geometry interpolated via ainterp_aux) and collapsed to the
% k disc-node rows by the L2 projection ipw, exactly as on the self
% block. The source side uses the GGQ log intermediate rule xs1 on
% panel j (via ainterp1), which resolves the near-log singularity seen
% by aux targets clustered toward the shared endpoint; a smooth
% oversampled rule loses ~4 digits there. This reproduces the adaptive
% quadrature used by the Fortran chunkmatc_aux_od for neighbor blocks
% to machine precision.
%
% When corrections is true, the native smooth (Gauss-Legendre) block is
% subtracted so the result can be added on top of a pre-built smooth
% matrix.

rs = r(:,:,j); ds = d(:,:,j); d2s = d2(:,:,j); ns = n(:,:,j);
rt = r(:,:,i); dt = d(:,:,i); d2t = d2(:,:,i); nt = n(:,:,i);

if isempty(data)
    dd = [];
    dds = [];
else
    dd = data(:,:,i);
    dds = data(:,:,j);
end

op1 = opdims(1); op2 = opdims(2);

% source geometry at the intermediate-rule nodes
rfine = (ainterp1*(rs.')).';
dfine = (ainterp1*(ds.')).';
d2fine = (ainterp1*(d2s.')).';
dfinenrm = sqrt(sum(dfine.^2,1));
nfine = [dfine(2,:); -dfine(1,:)]./dfinenrm;

% target geometry at the naux aux nodes
rt_aux = (ainterp_aux*(rt.')).';
dt_aux = (ainterp_aux*(dt.')).';
d2t_aux = (ainterp_aux*(d2t.')).';
dt_aux_nrm = sqrt(sum(dt_aux.^2,1));
nt_aux = [dt_aux(2,:); -dt_aux(1,:)]./dt_aux_nrm;

srcinfo = [];
srcinfo.r = rfine; srcinfo.d = dfine; srcinfo.n = nfine;
srcinfo.d2 = d2fine;
targinfo = [];
targinfo.r = rt_aux; targinfo.d = dt_aux; targinfo.n = nt_aux;
targinfo.d2 = d2t_aux;
if ~isempty(data)
    srcinfo.data = (ainterp1*(dds.')).';
    targinfo.data = (ainterp_aux*(dd.')).';
else
    srcinfo.data = [];
    targinfo.data = [];
end

dsdt = dfinenrm(:).*whts1(:);
dsdtndim2 = repmat(dsdt(:).',op2,1);
dsdtndim2 = dsdtndim2(:);

% K shape: (naux*op1, nxs1*op2)
smatbig = fkern(srcinfo,targinfo);
smatbig = bsxfun(@times,smatbig,dsdtndim2.');
aux_block = smatbig*ainterp1kron;      % (naux*op1, k*op2)

% target-side ipw projection: naux*op1 -> k*op1
if op1 == 1
    submat = ipw*aux_block;
else
    submat = kron(ipw,eye(op1))*aux_block;
end

if nargin >= 16 && corrections
    srcinfo = [];
    srcinfo.r = rs;  srcinfo.d = ds;
    srcinfo.d2 = d2s; srcinfo.n = ns;
    srcinfo.data = dds;
    targinfo.r = rt; targinfo.d = dt; targinfo.n = nt;
    targinfo.d2 = d2t; targinfo.data = dd;
    wtsj = wtss(:,j);
    submat = submat - fkern(srcinfo,targinfo).*(wtsj(:).');
end

end

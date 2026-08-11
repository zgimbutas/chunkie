function submat = nearbuildmat(r,d,n,d2,data,i,j,fkern,opdims,...
                               ainterp_aux,ipw,ct,bw,tadap,wadap,...
                               corrections,wtss)
%CHNK.QUADGALERKIN.NEARBUILDMAT
%
% Assemble a near (neighbor or otherwise close) panel block for the
% Galerkin (aux-projection) backend, faithfully porting the off-
% diagonal treatment of the Fortran chunkmatc_aux_od: rows are
% evaluated at the naux auxiliary target nodes on panel i (geometry
% interpolated via ainterp_aux) by PER-TARGET ADAPTIVE Gaussian
% quadrature over the whole source panel j (chnk.adapgausswts, the
% chunkie analog of the reference's cadachunk), then collapsed to the
% k disc-node rows by the L2 projection ipw, exactly as on the self
% block.
%
% Adaptive integration is required here because fixed intermediate
% rules are only valid for log-type near fields: across a corner the
% D / S' kernels leave the poly+log class, and the aux targets cluster
% within ~5e-4 panel-lengths of the shared vertex, where a fixed rule
% loses ~7 digits (measured 3.6e-7 for the 'd' kernel at k=16).
%
% ct, bw   - order-k Legendre nodes and barycentric weights
% tadap, wadap - integration rule for the adaptive quadrature
%                (lege.exps(2k+1), as in chunkermat's adaptive path)
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

% target geometry at the naux aux nodes
rt_aux = (ainterp_aux*(rt.')).';
dt_aux = (ainterp_aux*(dt.')).';
d2t_aux = (ainterp_aux*(d2t.')).';
dt_aux_nrm = sqrt(sum(dt_aux.^2,1));
nt_aux = [dt_aux(2,:); -dt_aux(1,:)]./dt_aux_nrm;
if ~isempty(dd)
    dd_aux = (ainterp_aux*(dd.')).';
else
    dd_aux = [];
end

% per-aux-target adaptive integration over source panel j
aux_block = chnk.adapgausswts(r,d,n,d2,data,ct,bw,j,...
    rt_aux,dt_aux,nt_aux,d2t_aux,dd_aux,...
    fkern,opdims,tadap,wadap);

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
    targinfo = [];
    targinfo.r = rt; targinfo.d = dt; targinfo.n = nt;
    targinfo.d2 = d2t; targinfo.data = dd;
    wtsj = wtss(:,j);
    submat = submat - fkern(srcinfo,targinfo).*(wtsj(:).');
end

end

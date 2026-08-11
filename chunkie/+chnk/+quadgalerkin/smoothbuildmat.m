function submat = smoothbuildmat(r,d,n,d2,data,i,j,fkern,opdims,...
                                 ainterp_aux,ipw,...
                                 ts_src,whts_src,ainterp_src,ainterp_src_kron)
%CHNK.QUADGALERKIN.SMOOTHBUILDMAT
%
% Well-separated (target panel i, source panel j not adjacent to i)
% block assembly for the chunkmatc_aux Galerkin scheme, matching
% chunkmatc_aux_od in chunkmatc.f. For each of the naux auxiliary
% target nodes on panel i:
%
%   1. Interpolate target geometry to the aux nodes via ainterp_aux.
%   2. Interpolate source geometry on panel j to the oversampled source
%      rule nodes ts_src via ainterp_src.
%   3. Evaluate the smooth kernel between every (aux_target, src) pair
%      and integrate against the cardinal Legendre basis at the disc
%      source nodes (via ainterp_src_kron).
%
% This produces a (naux*op1) x (k*op2) auxiliary block. The aux-target
% L2 projection ipw (k x naux) then collapses the naux aux-target rows
% to k disc-target rows, yielding the standard k*op1 x k*op2 block.
%
% Unlike chnk.quadnative.buildmat, which evaluates kernel rows directly
% at the k disc targets, this routine evaluates at the naux aux targets
% and projects back, matching the Fortran reference's projection-
% everywhere design. For well-separated smooth blocks the result is
% numerically indistinguishable from native Gauss-Legendre. The smooth
% source rule is only valid for well-separated panels: neighbor blocks
% must use chnk.quadgalerkin.nearbuildmat, whose intermediate rule
% resolves the near-log singularity seen by aux targets clustered
% toward the shared endpoint.

rs = r(:,:,j); ds_src = d(:,:,j); d2s = d2(:,:,j);
rt = r(:,:,i); dt_targ = d(:,:,i); d2t = d2(:,:,i);

if isempty(data)
    dds = []; ddt = [];
else
    dds = data(:,:,j);
    ddt = data(:,:,i);
end

op1 = opdims(1); op2 = opdims(2);

% target geometry at naux aux nodes
rt_aux = (ainterp_aux*(rt.')).';
dt_aux = (ainterp_aux*(dt_targ.')).';
d2t_aux = (ainterp_aux*(d2t.')).';
dt_aux_nrm = sqrt(sum(dt_aux.^2,1));
nt_aux = [dt_aux(2,:); -dt_aux(1,:)]./dt_aux_nrm;

% source geometry at nsrc oversampled nodes
rs_aux = (ainterp_src*(rs.')).';
ds_aux = (ainterp_src*(ds_src.')).';
d2s_aux = (ainterp_src*(d2s.')).';
ds_aux_nrm = sqrt(sum(ds_aux.^2,1));
ns_aux = [ds_aux(2,:); -ds_aux(1,:)]./ds_aux_nrm;
dsdt_src = ds_aux_nrm(:).*whts_src(:);

srcinfo = []; srcinfo.r = rs_aux; srcinfo.d = ds_aux;
srcinfo.d2 = d2s_aux; srcinfo.n = ns_aux;
targinfo = []; targinfo.r = rt_aux; targinfo.d = dt_aux;
targinfo.d2 = d2t_aux; targinfo.n = nt_aux;
if ~isempty(dds)
    srcinfo.data = (ainterp_src*(dds.')).';
    targinfo.data = (ainterp_aux*(ddt.')).';
else
    srcinfo.data = [];
    targinfo.data = [];
end

% K shape: (naux*op1, nsrc*op2)
K = fkern(srcinfo,targinfo);

% scale source columns by ds/dt
dsdtndim2 = repmat(dsdt_src(:).',op2,1); dsdtndim2 = dsdtndim2(:);
K = bsxfun(@times,K,dsdtndim2.');

% source-side: interpolate cardinal basis at ts_src back to k disc sources
aux_block = K*ainterp_src_kron;        % (naux*op1, k*op2)

% target-side: ipw projection
if op1 == 1
    submat = ipw*aux_block;
else
    submat = kron(ipw,eye(op1))*aux_block;
end

end

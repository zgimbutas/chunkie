function submat = diagbuildmat(r,d,n,d2,data,i,fkern,opdims,...
                               xs0,whts0,ainterps0kron,ainterps0,...
                               ts_aux,ainterp_aux,ipw,corrections,wtss,indd,...
                               exact_aux_geo)
%CHNK.QUADGALERKIN.DIAGBUILDMAT
%
% Self-panel block assembly for the chunkmatc_aux Galerkin scheme.
%
% For each of the naux auxiliary target nodes ts_aux(inode), use the
% per-target singular source sub-rule (xs0{inode}, whts0{inode}) to
% integrate the layer-potential kernel against the k disc-node source
% basis functions. This builds a naux*opdims(1) x k*opdims(2) aux
% block. Apply the L2 projection ipw (k x naux) on the target side to
% collapse the aux rows to the k disc-node rows, producing a standard
% k*opdims(1) x k*opdims(2) self block.
%
% Target geometry at the aux nodes (r, d, d2, n) is interpolated from
% the k disc nodes via ainterp_aux. Source geometry is interpolated to
% xs0{inode} via ainterps0{inode}, as in the GGQ self block.

rs = r(:,:,i); ds = d(:,:,i); d2s = d2(:,:,i); ns = n(:,:,i);
if isempty(data)
    dd = [];
else
    dd = data(:,:,i);
end

[dim,k] = size(rs);
naux = size(ipw,2);
op1 = opdims(1); op2 = opdims(2);

% -- target geometry at the naux aux nodes (exact if provided)
if nargin >= 19 && ~isempty(exact_aux_geo)
    rt_aux  = exact_aux_geo.r (:,:,i);
    dt_aux  = exact_aux_geo.d (:,:,i);
    d2t_aux = exact_aux_geo.d2(:,:,i);
    nt_aux  = exact_aux_geo.n (:,:,i);
else
    rt_aux = (ainterp_aux*(rs.')).';        % dim x naux
    dt_aux = (ainterp_aux*(ds.')).';
    d2t_aux = (ainterp_aux*(d2s.')).';
    dt_aux_nrm = sqrt(sum(dt_aux.^2,1));
    nt_aux = [dt_aux(2,:); -dt_aux(1,:)]./dt_aux_nrm;
end
if ~isempty(dd)
    dd_aux = (ainterp_aux*(dd.')).';
else
    dd_aux = [];
end

aux_block = zeros(naux*op1, k*op2);

srcinfo = []; targinfo = [];

use_exact_src = nargin >= 19 && ~isempty(exact_aux_geo) && ...
                isfield(exact_aux_geo,'r_src');

for inode = 1:naux
    % -- per-target source rule (exact if provided, else interpolated)
    xj = xs0{inode};
    wj = whts0{inode};
    nptsj = numel(xj);
    ainterp_j = ainterps0{inode};            % nptsj x k

    if use_exact_src
        rs_j  = exact_aux_geo.r_src {i,inode};
        ds_j  = exact_aux_geo.d_src {i,inode};
        d2s_j = exact_aux_geo.d2_src{i,inode};
        dfinenrm = sqrt(sum(ds_j.^2,1));
        ns_j  = exact_aux_geo.n_src {i,inode};
    else
        rs_j = (ainterp_j*(rs.')).';             % dim x nptsj
        ds_j = (ainterp_j*(ds.')).';
        d2s_j = (ainterp_j*(d2s.')).';
        dfinenrm = sqrt(sum(ds_j.^2,1));
        ns_j = [ds_j(2,:); -ds_j(1,:)]./dfinenrm;
    end
    dsdt_j = dfinenrm(:).*wj(:);

    srcinfo.r = rs_j;  srcinfo.d = ds_j;
    srcinfo.d2 = d2s_j; srcinfo.n = ns_j;
    targinfo.r = rt_aux(:,inode);  targinfo.d = dt_aux(:,inode);
    targinfo.d2 = d2t_aux(:,inode); targinfo.n = nt_aux(:,inode);
    if isempty(dd)
        srcinfo.data = [];  targinfo.data = [];
    else
        srcinfo.data = (ainterp_j*(dd.')).';
        targinfo.data = dd_aux(:,inode);
    end

    % -- row at this aux target, of shape op1 x nptsj*op2
    row_aux = fkern(srcinfo,targinfo);
    dsdtndim2 = repmat(dsdt_j(:).',op2,1); dsdtndim2 = dsdtndim2(:);
    row_aux = bsxfun(@times,row_aux,dsdtndim2.');

    % -- source-side interpolation: nptsj*op2 -> k*op2
    ainterp_j_kron = ainterps0kron{inode};   % nptsj*op2 x k*op2

    aux_block(op1*(inode-1)+1:op1*inode, :) = row_aux*ainterp_j_kron;
end

% -- target-side L2 projection: naux*op1 -> k*op1 via kron(ipw, I_op1)
if op1 == 1
    submat = ipw*aux_block;
else
    submat = kron(ipw,eye(op1))*aux_block;
end

% -- optional: return only the correction to the smooth (native) baseline
%    submat_corr = submat - K(disc_targ, disc_src) .* wts_disc(src)
if nargin >= 16 && corrections
    srcinfo.r = rs;  srcinfo.d = ds;
    srcinfo.d2 = d2s; srcinfo.n = ns;
    targinfo.r = rs; targinfo.d = ds;
    targinfo.d2 = d2s; targinfo.n = ns;
    targinfo.data = dd;  srcinfo.data = dd;
    sc = fkern(srcinfo,targinfo);
    sc(indd) = 0;             % drop the singular diagonal of the smooth K
    wtsi = wtss(:,i);
    submat = submat - sc.*(wtsi(:).');
end

end

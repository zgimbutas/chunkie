function submat = diagbuildmat(r,d,n,d2,data,i,fkern,opdims,...
                               xs0,whts0,ainterps0kron,ainterps0)
%CHNK.QUADGALERKIN.DIAGBUILDMAT
%
% Assemble the self-panel block for the Galerkin (aux-projection)
% backend. Per-target singular sub-rules xs0{j}/whts0{j} are applied
% at each of the k disc target nodes on panel i; the source data
% (r,d,d2,n,...) is interpolated from the k Legendre disc nodes to
% the singular-rule nodes via ainterps0{j}.
%
% Structurally identical to chnk.quadggq.diagbuildmat. The Galerkin
% L2 projection collapses to the standard interpolation-matrix path
% on the self block because each target uses its own dedicated rule;
% the aux-target oversampling is reserved for the off-diagonal
% smoothbuildmat path (see chnk.quadgalerkin.smoothbuildmat).

rs = r(:,:,i); ds = d(:,:,i); d2s = d2(:,:,i);
ns = n(:,:,i);
if isempty(data)
    dd = [];
else
    dd = data(:,:,i);
end

[~,k] = size(rs);

rfine = cell(k,1);
dfine = cell(k,1);
nfine = cell(k,1);
d2fine = cell(k,1);
dsdt = cell(k,1);
ddfine = cell(k,1);

for j = 1:k
    rfine{j} = (ainterps0{j}*(rs.')).';
    dfine{j} = (ainterps0{j}*(ds.')).';
    d2fine{j} = (ainterps0{j}*(d2s.')).';
    dfinenrm = sqrt(sum(dfine{j}.^2,1));
    nfine{j} = [dfine{j}(2,:); -dfine{j}(1,:)]./dfinenrm;
    dsdt{j} = (dfinenrm(:)).*whts0{j};
    if ~isempty(data)
        ddfine{j} = ((ainterps0{j}*(dd.'))).';
    end
end

srcinfo = [];
targinfo = [];

submat = zeros(k*opdims(1), k*opdims(2));

for j = 1:k
    srcinfo.r = rfine{j};  srcinfo.d = dfine{j};
    srcinfo.d2 = d2fine{j}; srcinfo.n = nfine{j};
    targinfo.r = rs(:,j);  targinfo.d = ds(:,j);
    targinfo.d2 = d2s(:,j); targinfo.n = ns(:,j);
    if isempty(dd)
        targinfo.data = [];
        srcinfo.data = [];
    else
        srcinfo.data = ddfine{j};
        targinfo.data = dd(:,j);
    end

    smatbigi = fkern(srcinfo,targinfo);
    dsdtndim2 = repmat(dsdt{j}.',opdims(2),1);
    dsdtndim2 = dsdtndim2(:);
    smatbigi = bsxfun(@times,smatbigi,dsdtndim2.');
    submat(opdims(1)*(j-1)+1:opdims(1)*j,:) = ...
        smatbigi*ainterps0kron{j};
end

end

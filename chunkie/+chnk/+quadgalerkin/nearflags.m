function nearf = nearflags(r,ainterp_aux,fac)
%CHNK.QUADGALERKIN.NEARFLAGS flag panel pairs needing adaptive quadrature
%
% nearf(it,js) is true when any auxiliary target node of panel it lies
% within fac * radius(js) of the center of source panel js. This is the
% near/far decision boundary of the chunkmatc family: the cap2Dsolver
% FMM near-field lists use fac = 2.4 (chunkfmm2d0npairs, written
% 1.2d0*2), which certifies ~1e-13 for the smooth far rule at npols=10
% and far better for the 2*k-point oversampled far rule used here.
%
% The diagonal is returned false; adjacency exclusions are the
% caller's responsibility.

[~,~,nch] = size(r);
naux = size(ainterp_aux,1);

rta = zeros(2,naux,nch);
cen = zeros(2,nch);
rad = zeros(1,nch);
for j = 1:nch
    rta(:,:,j) = (ainterp_aux*(r(:,:,j).')).';
    cen(:,j) = mean(r(:,:,j),2);
    rad(j) = sqrt(max(sum((r(:,:,j)-cen(:,j)).^2,1)));
end

nearf = false(nch,nch);
rtaf = reshape(rta,2,[]);            % 2 x naux*nch
for j = 1:nch
    d2 = sum((rtaf-cen(:,j)).^2,1);
    dmin = sqrt(min(reshape(d2,naux,nch),[],1));
    nearf(:,j) = dmin(:) <= fac*rad(j);
end
nearf(1:nch+1:end) = false;

end

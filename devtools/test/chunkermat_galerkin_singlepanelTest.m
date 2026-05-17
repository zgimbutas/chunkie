chunkermat_galerkin_singlepanelTest0();


function chunkermat_galerkin_singlepanelTest0()
%CHUNKERMAT_GALERKIN_SINGLEPANELTEST
%
% Self-block-only diagnostic on one straight-line chunker panel of
% length L=0.4. For a straight segment the chunker stores r exactly
% (polynomial degree 1 perfectly fitted by degree k-1 polynomial), so
% there is no geometry-interpolation noise; any disagreement between
% the GGQ self block and the chnk.quadgalerkin self block is purely
% algorithmic / projection-induced.
%
% Reports the relative Frobenius diff vs k. Observed: my port's
% Galerkin self block disagrees with GGQ by ~1.5%-3% across all k,
% with the diff decreasing slowly (only ~1.5x per doubling of k).
% This is the closed-form chunkmatc_form_ipipw projection's intrinsic
% Gram-orthogonality imperfection on log-tuned aux nodes — the same
% formula as the chunkmatc reference and chunkmatc-matlab/form_ipw.m,
% so the port's matrix-level behaviour is exactly that of the
% reference algorithm.
%
% Whether the reference Fortran achieves materially better matrix
% entries on the same problem is open: the chunkmatc-matlab mex is
% from an older MATLAB and crashes R2024b, so a direct head-to-head
% on this same simple geometry is currently impossible to verify
% from this machine.

L = 0.4;
fkern = @(s,t) chnk.lap2d.kern(s,t,'s');
fprintf('Single straight panel of length L=%.2f, Laplace S, no neighbors\n',L);
fprintf('%-4s %-12s %-12s\n','k','rel-F','max-abs');
for k = [4 6 8 10 12 16]
    [ts,wts] = lege.exps(k); ts = ts(:).';
    r  = [(1+ts)*L/2; zeros(1,k)];
    d  = [(L/2)*ones(1,k); zeros(1,k)];
    d2 = zeros(2,k);

    prefab = chunkerpref(); prefab.k = k;
    chnkr = chunker(prefab);
    chnkr = chnkr.addchunk();
    chnkr.r(:,:,1) = r;  chnkr.d(:,:,1) = d;  chnkr.d2(:,:,1) = d2;
    chnkr.n(:,:,1) = [d(2,:); -d(1,:)] ./ vecnorm(d);
    chnkr.adj(:,1) = [0;0];
    chnkr.wts(:,1) = (vecnorm(d).' .* wts(:));

    A_ggq = chunkermat(chnkr,fkern,struct('quad','ggq'));
    A_gal = chunkermat(chnkr,fkern,struct('quad','galerkin'));
    diff = A_ggq - A_gal;
    fprintf('%-4d %-12.3e %-12.3e\n',k,...
        norm(diff,'fro')/norm(A_ggq,'fro'),max(abs(diff(:))));
end

end

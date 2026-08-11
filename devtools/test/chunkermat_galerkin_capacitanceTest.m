chunkermat_galerkin_capacitanceTest0();


function chunkermat_galerkin_capacitanceTest0()
%CHUNKERMAT_GALERKIN_CAPACITANCETEST
%
% Replicates chunkmatc-matlab/test1cap.m: capacitance of the unit
% square via the first-kind integral equation S sigma = 1 (Laplace
% interior Dirichlet, single-layer representation, constant boundary
% data). The chunkmatc-matlab reference value (their normalization,
% S = -log r) is cap = -6.03125250071049; chunkie's normalization
% (S = -log r / (2 pi)) gives cap = 2*pi * (-6.03125250071049).
%
% The test sweeps cparams.depth (corner dyadic refinement, controllable
% in chunkerpoly via cparams.depth) and reports the convergence rates
% of both backends. RCIP is OFF.
%
% Observed: both backends converge spectrally to the reference, with
% Galerkin ~2x ahead of GGQ at each depth (1.9e-11 vs 4.4e-11 at
% depth=20). An earlier version of the Galerkin port stalled at ~2e-4
% because neighbor blocks were assembled with a smooth oversampled
% source rule, which cannot resolve the near-log singularity seen by
% the aux targets clustered toward the shared panel endpoint; that
% error is scale-invariant under dyadic refinement and so never
% decayed. Neighbor blocks now use the GGQ log intermediate rule
% evaluated at the aux targets followed by the same ipw projection
% (chnk.quadgalerkin.nearbuildmat), which matches the reference
% Fortran's adaptive off-diagonal quadrature to machine precision.

cap_chunkmatc_ref = -6.03125250071049;
cap_ref = cap_chunkmatc_ref * (2*pi);   % chunkie convention

square = [-1 1 1 -1; -1 -1 1 1];
fkern_s = @(s,t) chnk.lap2d.kern(s,t,'s');

fprintf('square cap (chunkie convention) reference = %.14e\n',cap_ref);
fprintf('%-8s %-8s %-14s %-14s %-14s %-14s\n','depth','nch','ggq cap','ggq err','gal cap','gal err');
e_ggq_last = 1; e_gal_last = 1;
for depth = [0 4 8 12 20]
    pref = []; pref.k = 10;
    cparams = []; cparams.eps = 1e-9; cparams.nover = 1;
    cparams.rounded = false; cparams.depth = depth;
    chnkr = chunkerpoly(square,cparams,pref);
    npts = chnkr.k*chnkr.nch;
    rhs = ones(npts,1);

    S_ggq = chunkermat(chnkr,fkern_s,struct('quad','ggq','rcip',false));
    sigma_ggq = S_ggq\rhs;
    cap_ggq = sum(chnkr.wts(:).*sigma_ggq);
    e_ggq = abs(cap_ggq - cap_ref);

    S_gal = chunkermat(chnkr,fkern_s,struct('quad','galerkin','rcip',false));
    sigma_gal = S_gal\rhs;
    cap_gal = sum(chnkr.wts(:).*sigma_gal);
    e_gal = abs(cap_gal - cap_ref);

    fprintf('%-8d %-8d %-14.10f %-14.2e %-14.10f %-14.2e\n', ...
        depth, chnkr.nch, cap_ggq, e_ggq, cap_gal, e_gal);
    e_ggq_last = e_ggq; e_gal_last = e_gal;
end

% both backends should reach machine precision at the deepest refinement
assert(e_ggq_last < 1e-9, 'ggq did not converge to capacitance reference');
assert(e_gal_last < 1e-9, 'galerkin did not converge to capacitance reference');

end

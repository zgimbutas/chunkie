chunkermat_galerkin_convergenceTest0();


function chunkermat_galerkin_convergenceTest0()
%CHUNKERMAT_GALERKIN_CONVERGENCETEST
%
% Convergence study of the chunkmatc_aux Galerkin backend on a smooth
% starfish under panel refinement. The chnk.quadgalerkin port (and
% the reference Fortran chunkmatc_aux it mirrors) uses the closed-form
% chunkmatc_form_ipipw projection formula
%   ipw(j,i) = (whts_aux(j)/whts_disc(i)) * L_i(ts_aux(j))
% which is the correct weighted-L2 projection only when the auxiliary
% nodes diagonalize the cardinal-Legendre Gram matrix. The log-tuned
% legeexps_log_lr nodes diagonalize G to ~9e-5, leaving residual error
% per matrix entry that produces an O(h) convergence floor.
%
% This test documents that floor: GGQ delivers machine precision
% uniformly while Galerkin converges roughly linearly in panel count
% under refinement. To recover deeper accuracy from the Fortran-style
% chunkmatc_aux scheme one must heavily refine panels (dyadic refine
% in chunkmatc reference tests), not change the formula.

zk = 0.1;
fkern = @(s,t) chnk.helm2d.kern(zk,s,t,'s');
fcurve = @(t) starfish(t,3,0.25);

pref = []; pref.k = 16;

rng(8675309);
ns = 5;
sources = 3.0*fcurve(2*pi*rand(ns,1));
strengths = randn(ns,1);
targets = 0.2*randn(2,3);
srcinfo = []; srcinfo.r = sources;
targinfo = []; targinfo.r = targets;
utarg = fkern(srcinfo,targinfo)*strengths;

fprintf('%-8s %-6s %-12s %-12s\n','eps','nch','ggq','galerkin');
err_gal_prev = nan;
nch_prev = nan;
for cps_eps = [1e-3 1e-6 1e-9 1e-12 1e-14]
    cparams = []; cparams.eps = cps_eps; cparams.nover = 1;
    chnkr = chunkerfunc(fcurve,cparams,pref);
    targinfo.r = chnkr.r;
    ubdry = fkern(srcinfo,targinfo)*strengths;

    A_ggq = chunkermat(chnkr,fkern,struct('quad','ggq'));
    A_gal = chunkermat(chnkr,fkern,struct('quad','galerkin'));
    sig_ggq = A_ggq\ubdry(:);
    sig_gal = A_gal\ubdry(:);

    src_all = []; src_all.r = reshape(chnkr.r,2,[]);
    w_all = chnkr.wts(:);
    targinfo.r = targets;
    upred_ggq = fkern(src_all,targinfo)*(sig_ggq.*w_all);
    upred_gal = fkern(src_all,targinfo)*(sig_gal.*w_all);
    e_ggq = norm(upred_ggq-utarg,inf)/norm(utarg,inf);
    e_gal = norm(upred_gal-utarg,inf)/norm(utarg,inf);
    fprintf('%-8.0e %-6d %-12.2e %-12.2e\n', cps_eps, chnkr.nch, e_ggq, e_gal);

    assert(e_ggq < 1e-10, sprintf('ggq error too large at nch=%d',chnkr.nch));

    err_gal_prev = e_gal;
    nch_prev = chnkr.nch;
end

% At the finest resolution above the Galerkin solver should reach at
% least ~1e-7. Loosen the assertion if the geometry / wavenumber are
% changed.
assert(err_gal_prev < 1e-7, ...
    sprintf('galerkin error %.2e at finest refinement exceeds 1e-7',err_gal_prev));

end

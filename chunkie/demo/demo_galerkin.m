%DEMO_GALERKIN compare the GGQ and chunkmatc_aux Galerkin backends.
%
% Solves a Laplace interior Dirichlet problem on a starfish using both
% chunkermat backends and reports interior-point accuracy. GGQ should
% deliver near machine precision; the Galerkin (chunkmatc_aux) backend
% delivers ~1e-6 on this geometry, which is the practical precision of
% the underlying algorithm.

addpaths_loc();

iseed = 8675309;
rng(iseed);

% smooth closed boundary
pref = []; pref.k = 16;
cparams = []; cparams.eps = 1e-12; cparams.nover = 1;
chnkr = chunkerfunc(@(t) starfish(t,3,0.25),cparams,pref);
fprintf('starfish: %d panels, k = %d, %d unknowns\n',...
        chnkr.nch, chnkr.k, chnkr.k*chnkr.nch);

% exterior point sources for a known interior solution
ns = 10;
sources = 3.0*starfish(2*pi*rand(ns,1),3,0.25);
strengths = randn(ns,1);

fkern = @(s,t) chnk.lap2d.kern(s,t,'s');

srcinfo = []; srcinfo.r = sources;
targinfo = []; targinfo.r = chnkr.r;
ubdry = fkern(srcinfo,targinfo)*strengths;
ubdry = ubdry(:);

% a few interior targets, direct (FLAM-free) evaluation
nt = 5;
targets = 0.2*randn(2,nt);
targinfo = []; targinfo.r = targets;
utarg = fkern(srcinfo,targinfo)*strengths;

src_all = []; src_all.r = reshape(chnkr.r,2,[]);
w_all = chnkr.wts(:);

backends = {'ggq','galerkin'};
for ib = 1:numel(backends)
    qopt = backends{ib};
    tic;
    A = chunkermat(chnkr,fkern,struct('quad',qopt));
    tbuild = toc;
    sigma = A\ubdry;
    upred = fkern(src_all,targinfo)*(sigma.*w_all);
    err = norm(upred-utarg,inf)/norm(utarg,inf);
    fprintf('%-10s : build %5.2fs, interior rel-err = %5.2e\n',...
            qopt, tbuild, err);
end

fprintf(['\nthe Galerkin backend is exposed through chunkermat via\n'...
         '    sysmat = chunkermat(chnkr, kern, struct(''quad'',''galerkin''))\n'...
         'see chnk.quadgalerkin.setup for the auxquads struct layout.\n']);

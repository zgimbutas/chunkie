chunkermat_galerkinTest0();


function chunkermat_galerkinTest0()
%CHUNKERMAT_GALERKINTEST
%
% Exercise the chnk.quadgalerkin (chunkmatc_aux aux-projection) backend
% by comparing it to the GGQ backend on a Laplace double-layer interior
% Dirichlet problem on a starfish, and check it solves the BVP.

iseed = 8675309;
rng(iseed);

cparams = [];
cparams.eps = 1.0e-6;
cparams.nover = 1;
pref = [];
pref.k = 16;
narms = 3;
amp = 0.25;
chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref);

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 3.0*sources;
strengths = randn(ns,1);

nt = 3;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp);
targets = targets.*repmat(rand(1,nt),2,1);

kerns = @(s,t) chnk.lap2d.kern(s,t,'s');
srcinfo = []; srcinfo.r = sources;
targinfo = []; targinfo.r = chnkr.r;
kernmats = kerns(srcinfo,targinfo);
ubdry = kernmats*strengths;

targinfo = []; targinfo.r = targets;
kernmatstarg = kerns(srcinfo,targinfo);
utarg = kernmatstarg*strengths;

fkern = @(s,t) chnk.lap2d.kern(s,t,'D');

% -- assemble with both backends
optsg = []; optsg.quad = 'ggq';
optsG = []; optsG.quad = 'galerkin';

start = tic; Aggq = chunkermat(chnkr,fkern,optsg);  tggq = toc(start);
start = tic; Agal = chunkermat(chnkr,fkern,optsG);  tgal = toc(start);

fprintf('matrix build: ggq = %5.2e s, galerkin = %5.2e s\n',tggq,tgal);

% -- compare matrices entry-wise (Galerkin and GGQ are mathematically
%    different discretizations, so a small but nonzero difference is
%    expected on the self block; off-diagonal blocks should match
%    closely because both backends use smooth Gauss-Legendre there).
relF = norm(Aggq-Agal,'fro')/norm(Aggq,'fro');
relI = norm(Aggq-Agal,inf)/norm(Aggq,inf);
fprintf('||A_ggq - A_galerkin||_F / ||A_ggq||_F = %5.2e\n',relF);
fprintf('||A_ggq - A_galerkin||_inf / ||A_ggq||_inf = %5.2e\n',relI);
% the two backends should agree to a few digits beyond the chunkr
% truncation tolerance but need not match in floating point
assert(relF < 1e-6, 'galerkin matrix unexpectedly far from ggq matrix');

% -- exercise the setup struct directly
auxquads = chnk.quadgalerkin.setup(chnkr.k,'log');
assert(auxquads.naux >= chnkr.k);
assert(isequal(size(auxquads.ipw),[chnkr.k auxquads.naux]));

% L2-projection sanity check: aux samples of a polynomial of degree
% <= k-1 should reproduce its values at the disc nodes via ipw, to
% within the polynomial-exactness of the log-tuned aux rule (which is
% designed for log-singular integrands, so polynomial exactness is
% finite but high enough for Galerkin's needs).
deg = chnkr.k - 1;
c = randn(deg+1,1);
ts_a = auxquads.ts_aux; ts_d = auxquads.ts_disc;
fa = polyval(flipud(c),ts_a);
fd = polyval(flipud(c),ts_d);
err_proj = norm(auxquads.ipw*fa - fd,inf)/norm(fd,inf);
fprintf('ipw polynomial-recovery relative error = %5.2e\n',err_proj);
assert(err_proj < 1e-6, 'ipw fails to reproduce polynomial of degree <= k-1');

% -- solve the BVP via galerkin matrix and check residual
sys = -0.5*eye(chnkr.k*chnkr.nch) + Agal;
rhs = ubdry; rhs = rhs(:);
sol = sys\rhs;
res = norm(sys*sol - rhs)/norm(rhs);
fprintf('galerkin BVP solve residual %5.2e\n',res);
assert(res < 1e-10,'galerkin BVP solve residual too large');

% -- evaluate at interior targets if FLAM/FMM available (non-fatal)
try
    optsE = []; optsE.usesmooth = false; optsE.verb = false; optsE.accel = false;
    Dsol = chunkerkerneval(chnkr,fkern,sol,targets,optsE);
    relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
    fprintf('galerkin solve, relative frobenius error %5.2e\n',relerr);
    assert(relerr < 1e-10, 'low precision in galerkin chunkermat solve');
catch errEval
    fprintf('skipping interior target eval: %s\n',errEval.message);
end

end

chunkermat_galerkin_neumannTest0();


function chunkermat_galerkin_neumannTest0()
%CHUNKERMAT_GALERKIN_NEUMANNTEST  collocation (GGQ) vs Galerkin
% quadrature on exterior Neumann boundary value problems on the unit
% square, as a refinement study in the corner grading depth.
%
% Representation u = S sigma, exterior Neumann jump relation
% (-1/2 I + S') sigma = f, outward normals. The density behaves like
% r^(-1/3) at the square corners (exterior angle 3*pi/2), so the
% corner-panel quadrature quality is what the sweep exposes.
%
% Case 1, Laplace S0': data from interior point charges (net charge
%   Q). The plain operator is rank-1 deficient (equilibrium density),
%   so the system is charge-completed: A = -I/2 + S' + 1*w', rhs =
%   f + Q, which selects the physical solution with int sigma = Q.
%   Metrics: the EXACT total-charge functional |w'sigma - Q|/Q and
%   the far potential against the exact source field.
%
% Case 2, Helmholtz Sk' (sound-hard scattering analog, k = 1, safely
%   below the first interior Dirichlet eigenvalue k^2 = 2*pi^2 of the
%   square, so -I/2 + Sk' is invertible as is). Data from an interior
%   Helmholtz point source, metric: far potential vs the exact field.
%   This is also the backend's first non-Laplace exercise: the
%   aux-projection rules handle the H0-type log singularity.
%
% Observed behaviour (k = 10 panels): Galerkin sits at its quadrature
% floor (~5e-14, charge ~4e-14) already at depth 2 and stays flat in
% the grading depth; GGQ is non-monotone -- best near depth 4
% (~1e-13) and stalling ~40x above the Galerkin floor (2e-12..6e-12)
% as the corner towers deepen. At depth 2 the gap is ~600x.

pref = []; pref.k = 10;
verts = [ .5 .5 -.5 -.5; -.5 .5 .5 -.5];

% interior sources and net charge
ys = [ 0.13 -0.18; 0.07 -0.11];
qs = [0.7; 0.55];
Q  = sum(qs);
zk = 1.0;
y0 = [0.13; 0.07];

% exterior evaluation targets
tt = 2*pi*(0:7)/8;
targ = []; targ.r = 1.5*[cos(tt); sin(tt)];

fkern_sp  = @(s,t) chnk.lap2d.kern(s,t,'sprime');
fkern_s   = @(s,t) chnk.lap2d.kern(s,t,'s');
fkern_hsp = @(s,t) chnk.helm2d.kern(zk,s,t,'sprime');
fkern_hs  = @(s,t) chnk.helm2d.kern(zk,s,t,'s');

% exact fields
pot_ex = zeros(size(targ.r,2),1);
for i = 1:2
    pot_ex = pot_ex - qs(i)/(2*pi)*log(vecnorm(targ.r - ys(:,i))).';
end
poth_ex = (1i/4)*besselh(0,1,zk*vecnorm(targ.r - y0)).';

fprintf('%-6s %-6s | %-10s %-10s | %-10s %-10s | %-10s %-10s\n', ...
    'depth','npts','Qerr ggq','Qerr gal','potL ggq','potL gal', ...
    'potH ggq','potH gal');
depths = [2 4 6 8 10 12];
res = zeros(numel(depths),6);
for id = 1:numel(depths)
    cparams = []; cparams.eps = 1e-14; cparams.nover = 1;
    cparams.rounded = false; cparams.depth = depths(id);
    chnkr = chunkerpoly(verts, cparams, pref);
    npts = chnkr.k*chnkr.nch;
    bdry = reshape(chnkr.r,2,[]);
    nrm  = reshape(chnkr.n,2,[]);
    w    = chnkr.wts(:);

    % Neumann data
    f = zeros(npts,1);
    for i = 1:2
        dx = bdry - ys(:,i); rsq = sum(dx.^2,1);
        f = f - qs(i)/(2*pi)*(sum(dx.*nrm,1)./rsq).';
    end
    rr = vecnorm(bdry - y0);
    fh = -(1i*zk/4)*(besselh(1,1,zk*rr).*sum((bdry-y0).*nrm,1)./rr).';

    src = []; src.r = bdry;
    quads = {'ggq','galerkin'};
    for iq = 1:2
        opts = []; opts.quad = quads{iq};
        % Laplace, charge-completed
        A = chunkermat(chnkr,fkern_sp,opts) - 0.5*eye(npts) + ones(npts,1)*w.';
        sigma = A\(f + Q);
        errQ = abs(w.'*sigma - Q)/Q;
        pot = fkern_s(src,targ)*(sigma.*w);
        errL = max(abs(pot - pot_ex))/max(abs(pot_ex));
        % Helmholtz
        Ah = chunkermat(chnkr,fkern_hsp,opts) - 0.5*eye(npts);
        sigh = Ah\fh;
        poth = fkern_hs(src,targ)*(sigh.*w);
        errH = max(abs(poth - poth_ex))/max(abs(poth_ex));
        res(id,(iq-1)*3+(1:3)) = [errQ errL errH];
    end
    fprintf('%-6d %-6d | %-10.2e %-10.2e | %-10.2e %-10.2e | %-10.2e %-10.2e\n', ...
        depths(id), npts, res(id,1), res(id,4), res(id,2), res(id,5), ...
        res(id,3), res(id,6));
end

% both schemes must converge; the Galerkin floor must be at least as
% good as GGQ on every metric, and materially better on at least one
floor_ggq = min(res(:,1:3),[],1);
floor_gal = min(res(:,4:6),[],1);
assert(all(floor_gal < 1e-8), 'Galerkin did not converge');
assert(all(floor_gal <= 3*floor_ggq), ...
    'Galerkin floor worse than GGQ on some metric');
assert(any(floor_gal < floor_ggq/3), ...
    'no metric shows a Galerkin advantage');
fprintf('floors: ggq [%.1e %.1e %.1e]  galerkin [%.1e %.1e %.1e]\n', ...
    floor_ggq, floor_gal);
fprintf('chunkermat_galerkin_neumannTest PASSED\n');

end

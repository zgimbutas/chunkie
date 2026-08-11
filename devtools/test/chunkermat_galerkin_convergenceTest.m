chunkermat_galerkin_convergenceTest0();


function chunkermat_galerkin_convergenceTest0()
%CHUNKERMAT_GALERKIN_CONVERGENCETEST
%
% Corner-refinement convergence study comparing chnk.quadgalerkin
% (chunkmatc_aux aux-projection) and GGQ on the exterior Neumann
% Laplace problem (-1/2 I + S') sigma = du_inc/dn on a sharp 10-corner
% polygon, varying cparams.depth (dyadic-refinement depth into corners
% inside chunkerpoly). RCIP is disabled for both backends.
%
% With the adaptive near-block quadrature (chunkmatc_aux_od port) the
% Galerkin backend solves this corner problem to ~5e-14 already at
% depth 0 -- no corner refinement at all -- while GGQ needs depth 20
% to reach 2.7e-9: the chunkmatc Fortran reference's corner advantage,
% fully realized. The win ratio ranges from ~80x (depth 20) to ~6e8x
% (depth 0).

w = 1.0; h = 1.0; g = 0.05;
geo = -[ -w*g -w*g -w/2 -w/2 w/2 w/2 w*g w*g  w  w -w -w; ...
         -h   -.7*h -.7*h -h/3 -h/3 -.7*h -.7*h -h -h h  h -h];

source = [0.2; -0.1];
target = [1.0; 1.2];
charge = 1;
fkern_sp = @(s,t) chnk.lap2d.kern(s,t,'sp');
fkern_s  = @(s,t) chnk.lap2d.kern(s,t,'s');
pot_exact = (-charge/(2*pi)) * log(norm(target - source));

depths = [0 4 8 12 20];
fprintf('sharp 10-corner polygon, Sp Neumann, no RCIP\n');
fprintf('%-8s %-8s %-12s %-12s %-8s\n','depth','nch','ggq','galerkin','win');
e_ggq = zeros(1,numel(depths));
e_gal = zeros(1,numel(depths));
for id = 1:numel(depths)
    pref = []; pref.k = 10;
    cparams = []; cparams.eps = 1e-9; cparams.nover = 1;
    cparams.rounded = false;
    cparams.depth = depths(id);
    chnkr = chunkerpoly(geo,cparams,pref);
    npts = chnkr.k*chnkr.nch;

    boundary = reshape(chnkr.r,2,[]);
    normals  = reshape(chnkr.n,2,[]);
    dx = boundary - source;
    rsq = sum(dx.^2,1);
    rhs = (-charge/(2*pi)) * sum(dx .* normals, 1).' ./ rsq.';

    for iq = 1:2
        if iq == 1, qname = 'ggq'; else, qname = 'galerkin'; end
        A = chunkermat(chnkr,fkern_sp,struct('quad',qname,'rcip',false));
        A = A - 0.5*eye(npts);
        sigma = A\rhs;
        src_all = []; src_all.r = boundary;
        tgt = []; tgt.r = target;
        pot = fkern_s(src_all,tgt) * (sigma .* chnkr.wts(:));
        e = abs((pot - pot_exact)/pot_exact);
        if iq == 1, e_ggq(id) = e; else, e_gal(id) = e; end
    end
    win = e_ggq(id) / e_gal(id);
    fprintf('%-8d %-8d %-12.2e %-12.2e %-8.2f\n',...
        depths(id), chnkr.nch, e_ggq(id), e_gal(id), win);
end

% Across the realistic depths 0..20, Galerkin should be at least 2x
% better than GGQ on each refinement level (per chunkmatc reference)
ratios = e_ggq ./ e_gal;
assert(all(ratios > 2.0), ...
    'expected Galerkin to outperform GGQ by >2x on every depth 0..20');
end

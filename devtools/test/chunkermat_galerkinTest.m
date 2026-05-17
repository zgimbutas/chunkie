chunkermat_galerkinTest0();


function chunkermat_galerkinTest0()
%CHUNKERMAT_GALERKINTEST
%
% Side-by-side comparison of the chnk.quadgalerkin (chunkmatc_aux
% aux-projection Galerkin) backend and GGQ on the second-kind Laplace
% exterior-Neumann integral equation (-1/2 I + S') sigma = du_inc/dn,
% representation u = S sigma. The reference Fortran chunkmatc-matlab
% test2.m establishes the expected behaviour:
%
%   * Smooth boundaries: Galerkin and GGQ match in quality.
%   * Corner geometries: Galerkin slightly outperforms GGQ.

run_case('starfish (smooth)',  @smooth_starfish,    struct('eps',1e-12), 1e-9);
run_case('polygon 10 corners', @corner_poly,        struct('eps',1e-9),  1e-9);

end


function chnkr = corner_poly(cparams,pref)
    w = 1.0; h = 1.0; g = 0.05;
    geo = -[ -w*g -w*g -w/2 -w/2 w/2 w/2 w*g w*g  w  w -w -w; ...
             -h   -.7*h -.7*h -h/3 -h/3 -.7*h -.7*h -h -h h  h -h];
    chnkr = chunkerpoly(geo,cparams,pref);
end


function run_case(label, builder, cparam_opts, tol)
    pref = []; pref.k = 10;
    cparams = []; cparams.eps = cparam_opts.eps; cparams.nover = 1;
    chnkr = builder(cparams, pref);
    npts = chnkr.k*chnkr.nch;

    rng(8675309);
    source = [0.2; -0.1];        % interior to all test geometries
    target = [1.0; 1.2];         % exterior
    charge = 1;
    boundary = reshape(chnkr.r,2,[]);
    normals  = reshape(chnkr.n,2,[]);
    dx = boundary - source;
    rsq = sum(dx.^2,1);
    rhs = (-charge/(2*pi)) * sum(dx .* normals, 1).' ./ rsq.';

    fkern_sp = @(s,t) chnk.lap2d.kern(s,t,'sp');
    fkern_s  = @(s,t) chnk.lap2d.kern(s,t,'s');
    pot_exact = (-charge/(2*pi)) * log(norm(target - source));

    fprintf('%-22s (k=%d, nch=%d, npts=%d)\n',label,chnkr.k,chnkr.nch,npts);
    errs = struct();
    for q = {'ggq','galerkin'}
        A = chunkermat(chnkr,fkern_sp,struct('quad',q{1}));
        A = A - 0.5*eye(npts);
        sigma = A\rhs;
        src_all = []; src_all.r = boundary;
        tgt = []; tgt.r = target;
        pot = fkern_s(src_all,tgt) * (sigma .* chnkr.wts(:));
        err = abs((pot - pot_exact)/pot_exact);
        errs.(q{1}) = err;
        fprintf('  %-10s rel_err = %.2e\n',q{1},err);
        assert(err < tol, [label ' / ' q{1} ': error too large']);
    end
end


function chnkr = smooth_starfish(cparams,pref)
    chnkr = chunkerfunc(@(t) starfish(t,3,0.25),cparams,pref);
end



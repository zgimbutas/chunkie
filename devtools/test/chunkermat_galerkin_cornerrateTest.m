chunkermat_galerkin_cornerrateTest0();


function chunkermat_galerkin_cornerrateTest0()
%CHUNKERMAT_GALERKIN_CORNERRATETEST
%
% Corner LADDER RATE: the factor by which a functional's error drops per
% dyadic refinement level added at the corners. Singular-exponent theory
% predicts this factor exactly, so a measured rate can be attributed
% rather than just observed.
%
% For a Laplace corner the potential behaves like r^nu and the density
% like r^(nu-1), and a functional converges per added level at
%
%     2^nu       "unprotected"
%     2^(2 nu)   "protected"   -- a dual pairing squares the rate
%
% A 90-degree conductor corner has field arc 3*pi/2, so nu = 2/3, giving
% 1.5874 and 2.5198. Every geometry here is an axis-aligned polygon.
%
% Five cases compare collocation (GGQ) against the aux-projection
% Galerkin backend, all at k = 10 with RCIP off.
%
%   A  Square, first-kind S sigma = 1, charge functional.
%      Both backends 2.520. This pairing is self-dual, so collocation is
%      protected already; Galerkin wins a constant 2.5x, not a rate.
%
%   B  Square, combined field (1/2 I + D + S) sigma = 1, same charge
%      functional. Only the operator differs from A, and the branches
%      split: GGQ 1.585, Galerkin 2.520. B is also the well-conditioned
%      formulation (cond 7 -> 28 across these depths, against a
%      first-kind S that needs preconditioning), so under collocation
%      one trades conditioning against corner accuracy.
%
%   C  Coated square coax (dielectric): inner conductor at V=1, square
%      interface with eps 10 inside and air outside, outer conductor at
%      V=0. u = S sigma on all three curves; conductor rows set the
%      potential, interface rows impose continuity of eps du/dn.
%      Functional: charge on the inner conductor, LINEAR in sigma.
%      GGQ 1.70, Galerkin 2.65.
%
%   D  Same solve as C. Functional: pp1 = W1/(W1+W2), the fraction of
%      field energy stored in the coating, with W_r = (eps_r/2) oint
%      u du/dn by Green's identity. QUADRATIC in sigma and operator-
%      dependent -- it needs one-sided normal derivatives, so unlike the
%      charge it cannot be evaluated from sigma alone.
%      GGQ 1.68, Galerkin 2.73.
%
%   E  Slotted polygon, (-1/2 I + S') sigma = f, potential at an
%      exterior target. A robustness check on the GGQ rate: 1.588 again
%      on a different geometry, operator and functional, and at every
%      panel order (k = 4,6,8,10 all give 1.587-1.590) -- the rate
%      follows the corner exponent, not the quadrature. Galerkin has NO
%      measurable rate here: it starts at 5e-14 and sits on an unrelated
%      floor, so its r_gal column is that floor drifting, not error
%      growth. On a plain square this case does not discriminate at all
%      -- both backends floor out at once -- hence the slotted geometry.
%
% Reading across: A and B differ only in the operator, which isolates
% the split from geometry, functional and refinement. C then shows the
% linear charge functional -- protected under collocation in A -- losing
% that protection once the system carries S' interface rows. So what
% decides is the operator, not the functional; D confirms that a
% quadratic, operator-dependent functional fares no worse than a linear
% one. A through D all have both branches measurable; E stands apart as
% an independent check on the GGQ rate. The aux projection holds the
% protected branch throughout.
%
% Two practical notes on C and D. They assert BANDS rather than point
% values because two corner classes contribute: conductor corners with
% nu = 2/3 (1.587 / 2.520) and eps=10 coating corners with nu = 0.7317
% (1.661 / 2.759), and the measurements land near them. And since
% neither has a closed-form reference, the transmission solve and the
% region energies are first checked against the coated CIRCULAR coax,
% where both are exact; that check is part of the test, because on the
% square a sign error would pass as a plausible convergence rate.
%
% Not included: the plain electrostatic energy of conductors held at
% prescribed potential. W = (1/2) sum_k V_k Q_k is a combination of the
% charges, not an independent functional -- in case A the discrete form
% collapses onto Q to 5.6e-16, S sigma = 1 being the equation solved.
% Energy becomes independent only when decomposed by region, i.e. case D.
%
% ------------------------------------------------------------------
% REFERENCE OUTPUT (what this test prints when healthy)
%
%   case  depth   ggq err     gal err     r_ggq   r_gal
%   ----  -----   ---------   ---------   -----   -----
%   A     0       4.886e-03   1.984e-03       -       -
%   A     4       1.211e-04   4.919e-05   2.520   2.520
%   A     8       3.004e-06   1.220e-06   2.520   2.520
%   A     12      7.452e-08   3.026e-08   2.520   2.520
%   B     0       2.525e-02   1.393e-03       -       -
%   B     4       4.879e-03   3.448e-05   1.508   2.521
%   B     8       7.918e-04   8.557e-07   1.576   2.519
%   B     12      1.253e-04   2.293e-08   1.585   2.472
%   C     0       3.537e-03   2.266e-05       -       -
%   C     4       4.246e-04   4.624e-07   1.699   2.646
%   C     8       5.117e-05   9.387e-09   1.697   2.649
%   C     12      6.024e-06   1.294e-11   1.707   5.189 *
%   D     0       6.452e-04   2.280e-06       -       -
%   D     4       8.073e-05   4.082e-08   1.681   2.734
%   D     8       1.006e-05   6.391e-10   1.683   2.827
%   D     12      1.237e-06   7.369e-12   1.689   3.052 *
%   E     0       2.957e-05   5.031e-14       -       -
%   E     4       4.644e-06   5.960e-13   1.589   0.539 *
%   E     8       7.310e-07   1.875e-12   1.588   0.751 *
%   E     12      1.151e-07   5.068e-12   1.587   0.780 *
%
%   * r_gal at or past the Galerkin floor -- a floor drifting, not a
%     convergence rate. All of case E is in that regime; C and D enter
%     it at depth 12, which is why their rates are asserted on the
%     first two intervals only.
%
%   case B cond(A):  7.14  13.2  20.5  27.6   (identical for both
%                    backends -- same operator, different quadrature)
%
%   C/D circular-coax check vs closed form:  Q rel 4.4e-12,
%                    pp1 abs 5.3e-12, energy identity 1.0e-14
%
%   Runtime ~40 s. Numbers are k = 10 with RCIP off; they shift if
%   either changes.
% ------------------------------------------------------------------

nu = 2/3;
r_unprot = 2^nu;            % 1.5874
r_prot   = 2^(2*nu);        % 2.5198

fprintf(['nu = %.4f (90-deg corner)   unprotected 2^nu = %.4f   ' ...
         'protected 2^(2nu) = %.4f\n'], nu, r_unprot, r_prot);

depths = [0 4 8 12];

% ------------------------------------------------------------------
% Case A: first-kind S, charge functional -- self-dual, both protected
% ------------------------------------------------------------------
square = [-1 1 1 -1; -1 -1 1 1];
fkern_s = @(s,t) chnk.lap2d.kern(s,t,'s');
cap_ref = -6.03125250071049 * (2*pi);      % chunkie normalization

eA = zeros(2,numel(depths));
for id = 1:numel(depths)
    chnkr = build(square,depths(id));
    npts = chnkr.k*chnkr.nch;
    for iq = 1:2
        S = chunkermat(chnkr,fkern_s,struct('quad',qname(iq),'rcip',false));
        sigma = S\ones(npts,1);
        eA(iq,id) = abs(sum(chnkr.wts(:).*sigma) - cap_ref);
    end
end
rA = report('Case A: first-kind S, charge functional (self-dual)', ...
    depths, eA);

% ------------------------------------------------------------------
% Case B: combined field, same square/functional as A -- the two
% branches measured side by side on one problem. Use the
% +0.5I + D + S orientation: -0.5I - D + S converges to the same
% reference at the same rates but with cond(A) reaching ~1e11.
% ------------------------------------------------------------------
fkern_d = @(s,t) chnk.lap2d.kern(s,t,'d');
eB = zeros(2,numel(depths));
cB = zeros(2,numel(depths));
for id = 1:numel(depths)
    chnkr = build(square,depths(id));
    npts = chnkr.k*chnkr.nch;
    for iq = 1:2
        opts = struct('quad',qname(iq),'rcip',false);
        S = chunkermat(chnkr,fkern_s,opts);
        D = chunkermat(chnkr,fkern_d,opts);
        A = 0.5*eye(npts) + D + S;
        sigma = A\ones(npts,1);
        eB(iq,id) = abs(sum(chnkr.wts(:).*sigma) - cap_ref);
        cB(iq,id) = cond(A);
    end
end
rB = report('Case B: combined field (1/2 I + D + S), charge functional', ...
    depths, eB);
fprintf('cond(A)  ggq %s\n', sprintf('%-11.2e',cB(1,:)));
fprintf('cond(A)  gal %s\n', sprintf('%-11.2e',cB(2,:)));

% ------------------------------------------------------------------
% Cases C and D: coated coax (dielectric). First validate the
% transmission solve and the Green's-identity region energies against
% the closed-form COATED CIRCULAR coax, then run the cornered version.
% ------------------------------------------------------------------
ac = 1/3; cc = 0.6; bc = 1.4;        % bc ~= 1 avoids the log-capacity null
ep1 = 10.0; ep2 = 1.0;

Qx  = 2*pi/( log(cc/ac)/ep1 + log(bc/cc)/ep2 );
W1x = Qx^2*log(cc/ac)/(4*pi*ep1);
W2x = Qx^2*log(bc/cc)/(4*pi*ep2);
ppx = W1x/(W1x+W2x);

pref = []; pref.k = 16;
cp = []; cp.eps = 1e-12; cp.maxchunklen = 2*pi*bc/24;
circles = {@(t) circ(t,ac), @(t) circ(t,cc), @(t) circ(t,bc)};
chk = cell(1,3);
for j = 1:3, chk{j} = chunkerfunc(circles{j},cp,pref); end
[Qc,ppc,idc] = coax(chk{1},chk{2},chk{3},ep1,ep2,'galerkin');
fprintf(['\n=== Cases C/D setup check: coated CIRCULAR coax vs exact ===\n' ...
         'Q   %.12f  (exact %.12f, rel %.2e)\n' ...
         'pp1 %.12f  (exact %.12f, abs %.2e)\n' ...
         'energy identity |W1+W2 - Q/2| = %.2e\n'], ...
    Qc, Qx, abs(Qc-Qx)/abs(Qx), ppc, ppx, abs(ppc-ppx), idc);
assert(abs(Qc-Qx)/abs(Qx) < 1e-9, ...
    'coated-coax setup: charge off exact by %.2e', abs(Qc-Qx)/abs(Qx));
assert(abs(ppc-ppx) < 1e-9, ...
    'coated-coax setup: pp1 off exact by %.2e', abs(ppc-ppx));
assert(idc < 1e-11, ...
    'coated-coax setup: energy identity violated at %.2e', idc);

% cornered version. References are Richardson limits of the Galerkin
% sequence taken to depth 24; the circular check above is what
% establishes that the machinery producing them is right.
C_ref  = 7.563588120434;
pp_ref = 0.078144202559;
sq = @(hw) [-hw hw hw -hw; -hw -hw hw hw];

eC = zeros(2,numel(depths));
eD = zeros(2,numel(depths));
for id = 1:numel(depths)
    ck = cell(1,3); hw = [ac cc bc];
    for j = 1:3, ck{j} = build(sq(hw(j)),depths(id)); end
    for iq = 1:2
        [Cq,ppq] = coax(ck{1},ck{2},ck{3},ep1,ep2,qname(iq));
        eC(iq,id) = abs(Cq  - C_ref);
        eD(iq,id) = abs(ppq - pp_ref);
    end
end
rC = report('Case C: coated coax, charge functional (linear)', ...
    depths, eC);
rD = report('Case D: coated coax, region energy pp1 (quadratic)', ...
    depths, eD);

% ------------------------------------------------------------------
% Case E: Sp exterior Neumann on a slotted polygon -- collocation
% unprotected, Galerkin not corner-limited at all
% ------------------------------------------------------------------
w = 1.0; h = 1.0; g = 0.05;
poly = -[ -w*g -w*g -w/2 -w/2 w/2 w/2 w*g w*g  w  w -w -w; ...
          -h   -.7*h -.7*h -h/3 -h/3 -.7*h -.7*h -h -h h  h -h];
source = [0.2; -0.1];
target = [1.0; 1.2];
charge = 1;
fkern_sp = @(s,t) chnk.lap2d.kern(s,t,'sp');
pot_exact = (-charge/(2*pi)) * log(norm(target - source));

eE = zeros(2,numel(depths));
for id = 1:numel(depths)
    chnkr = build(poly,depths(id));
    npts = chnkr.k*chnkr.nch;
    boundary = reshape(chnkr.r,2,[]);
    normals  = reshape(chnkr.n,2,[]);
    dx = boundary - source;
    rhs = (-charge/(2*pi)) * sum(dx.*normals,1).' ./ sum(dx.^2,1).';
    for iq = 1:2
        A = chunkermat(chnkr,fkern_sp,struct('quad',qname(iq),'rcip',false));
        A = A - 0.5*eye(npts);
        sigma = A\rhs;
        src = []; src.r = boundary; tgt = []; tgt.r = target;
        pot = fkern_s(src,tgt)*(sigma.*chnkr.wts(:));
        eE(iq,id) = abs((pot - pot_exact)/pot_exact);
    end
end
rE = report('Case E: Sp exterior Neumann, potential functional', ...
    depths, eE);

% ------------------------------------------------------------------
% assertions
% ------------------------------------------------------------------

% A: both backends on the protected branch, and far from the
% unprotected one (the two branches differ by 59%, so 5% is ample
% separation while absorbing the drift near the reference floor)
for iq = 1:2
    assert(all(abs(rA(iq,:) - r_prot)/r_prot < 0.05), ...
        ['case A / ' qname(iq) ': rate %s not on the protected ' ...
         'branch 2^(2nu) = %.4f'], mat2str(round(rA(iq,:),3)), r_prot);
end

% B: the two branches split on one and the same problem. The first
% interval (depth 0->4) is still pre-asymptotic for GGQ (1.508), so the
% rates are asserted from depth 4 on.
asy = 2:size(rB,2);
assert(all(abs(rB(1,asy) - r_unprot)/r_unprot < 0.05), ...
    ['case B / ggq: rate %s not on the unprotected branch ' ...
     '2^nu = %.4f'], mat2str(round(rB(1,asy),3)), r_unprot);
assert(all(abs(rB(2,asy) - r_prot)/r_prot < 0.05), ...
    ['case B / galerkin: rate %s not on the protected branch ' ...
     '2^(2nu) = %.4f'], mat2str(round(rB(2,asy),3)), r_prot);

% B: and the combined field really is the well-conditioned one --
% otherwise the tradeoff this case demonstrates would not bite
assert(all(cB(:) < 1e3), ...
    'case B: combined field unexpectedly ill-conditioned, max %.2e', ...
    max(cB(:)));

% C/D: two corner classes contribute (conductor nu = 2/3 -> 1.587/2.520,
% eps=10 coating nu = 0.7317 -> 1.661/2.759), so assert disjoint BANDS
% spanning both predictions rather than a single value. Galerkin reaches
% its floor by depth 12, so its rates are read on the first two
% intervals only; GGQ is still converging throughout.
band_unprot = [1.50 1.85];
band_prot   = [2.40 2.95];
for cs = {{'C',rC},{'D',rD}}
    nm = cs{1}{1}; rr = cs{1}{2};
    assert(all(rr(1,:) > band_unprot(1) & rr(1,:) < band_unprot(2)), ...
        ['case ' nm ' / ggq: rate %s outside the unprotected band ' ...
         '[%.2f %.2f]'], mat2str(round(rr(1,:),3)), band_unprot);
    assert(all(rr(2,1:2) > band_prot(1) & rr(2,1:2) < band_prot(2)), ...
        ['case ' nm ' / galerkin: rate %s outside the protected band ' ...
         '[%.2f %.2f]'], mat2str(round(rr(2,1:2),3)), band_prot);
end

% E: collocation on the unprotected branch
assert(all(abs(rE(1,:) - r_unprot)/r_unprot < 0.05), ...
    ['case E / ggq: rate %s not on the unprotected branch ' ...
     '2^nu = %.4f'], mat2str(round(rE(1,:),3)), r_unprot);

% E: Galerkin not corner-limited -- orders of magnitude below GGQ at
% every depth, so no corner rate to measure
gap = eE(1,:)./eE(2,:);
assert(all(gap > 1e3), ...
    'case E: expected galerkin >1e3x below ggq at every depth, got %s', ...
    mat2str(round(gap)));

fprintf('\nchunkermat_galerkin_cornerrateTest PASSED\n');

end


function [Q,pp1,ident] = coax(ch1,ch2,ch3,eps1,eps2,quad)
%COAX  coated-coax transmission solve + per-region energies.
%
% u = S sigma on all three curves; conductor rows set the potential,
% interface rows impose continuity of eps du/dn. chunkie uses
% G = -log r/(2 pi) with outward normals, so du/dn|_pm = -+ sigma/2 +
% S'sigma (verified: S'1 = -1/2, exterior derivative -1, interior 0).
%
% Q    - charge on the inner conductor, -eps1 * oint du/dn
% pp1  - W1/(W1+W2), W_r = (eps_r/2) oint u du/dn with n out of region r
% ident- |W1+W2 - Q/2|, the energy identity (should be ~0)

chnkr = merge([ch1 ch2 ch3]);
n1 = ch1.k*ch1.nch; n2 = ch2.k*ch2.nch; n3 = ch3.k*ch3.nch;
N = n1+n2+n3;
i1 = 1:n1; i2 = n1+(1:n2); i3 = n1+n2+(1:n3);

opts = struct('quad',quad,'rcip',false);
S  = chunkermat(chnkr,@(s,t) chnk.lap2d.kern(s,t,'s') ,opts);
Sp = chunkermat(chnkr,@(s,t) chnk.lap2d.kern(s,t,'sp'),opts);
w  = chnkr.wts(:);

Ip1 = zeros(n1,N); Ip1(:,i1) = eye(n1);
Ip2 = zeros(n2,N); Ip2(:,i2) = eye(n2);
Ip3 = zeros(n3,N); Ip3(:,i3) = eye(n3);
dn1     = -0.5*Ip1 + Sp(i1,:);      % dielectric side of the inner cond.
dn2_out = -0.5*Ip2 + Sp(i2,:);      % region-2 side of the interface
dn2_in  =  0.5*Ip2 + Sp(i2,:);      % region-1 side of the interface
dn3     =  0.5*Ip3 + Sp(i3,:);      % dielectric side of the outer cond.

A = zeros(N,N); rhs = zeros(N,1);
A(i1,:) = S(i1,:);                    rhs(i1) = 1;
A(i3,:) = S(i3,:);                    rhs(i3) = 0;
A(i2,:) = eps1*dn2_in - eps2*dn2_out; rhs(i2) = 0;

sigma = A\rhs;
u = S*sigma;

Q = -eps1 * sum(w(i1).*(dn1*sigma));
W1 = 0.5*eps1*( -sum(w(i1).*(u(i1).*(dn1    *sigma))) ...
                +sum(w(i2).*(u(i2).*(dn2_in *sigma))) );
W2 = 0.5*eps2*( -sum(w(i2).*(u(i2).*(dn2_out*sigma))) ...
                +sum(w(i3).*(u(i3).*(dn3    *sigma))) );
pp1 = W1/(W1+W2);
ident = abs((W1+W2) - Q/2);
end


function [r,d,d2] = circ(t,R)
    t = t(:).';
    r  = R*[cos(t); sin(t)];
    d  = R*[-sin(t); cos(t)];
    d2 = -r;
end


function chnkr = build(verts,depth)
    pref = []; pref.k = 10;
    cparams = []; cparams.eps = 1e-9; cparams.nover = 1;
    cparams.rounded = false; cparams.depth = depth;
    chnkr = chunkerpoly(verts,cparams,pref);
end


function s = qname(iq)
    if iq == 1, s = 'ggq'; else, s = 'galerkin'; end
end


function rates = report(label,depths,e)
% per-level rate between successive depths, printed and returned
    fprintf('\n=== %s ===\n',label);
    fprintf('%-7s %-11s %-11s %-9s %-9s\n', ...
        'depth','ggq','galerkin','r_ggq','r_gal');
    rates = zeros(2,numel(depths)-1);
    fprintf('%-7d %-11.3e %-11.3e %-9s %-9s\n', ...
        depths(1), e(1,1), e(2,1), '-', '-');
    for id = 2:numel(depths)
        dl = depths(id) - depths(id-1);
        rates(:,id-1) = (e(:,id-1)./e(:,id)).^(1/dl);
        fprintf('%-7d %-11.3e %-11.3e %-9.3f %-9.3f\n', ...
            depths(id), e(1,id), e(2,id), rates(1,id-1), rates(2,id-1));
    end
end

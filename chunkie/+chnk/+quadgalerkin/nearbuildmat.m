function submat = nearbuildmat(r,d,n,d2,data,i,j,fkern,opdims,...
                               xs1,whts1,ainterp1kron,ainterp1,corrections,wtss)
%CHNK.QUADGALERKIN.NEARBUILDMAT
%
% Assemble a neighbor-panel block for the Galerkin (aux-projection)
% backend. Source data on panel j is interpolated to the intermediate
% rule nodes xs1 via ainterp1; the kernel is evaluated to the k disc
% targets on panel i. Structurally identical to chnk.quadggq.nearbuildmat.

rs = r(:,:,j); ds = d(:,:,j); d2s = d2(:,:,j); ns = n(:,:,j);
rt = r(:,:,i); dt = d(:,:,i); d2t = d2(:,:,i); nt = n(:,:,i);

if isempty(data)
    dd = [];
    dds = [];
else
    dd = data(:,:,i);
    dds = data(:,:,j);
end

rfine = (ainterp1*(rs.')).';
dfine = (ainterp1*(ds.')).';
d2fine = (ainterp1*(d2s.')).';
nfine = (ainterp1*(ns.')).';

if ~isempty(data)
    ddsfine = (ainterp1*(dds.')).';
end

srcinfo = [];
srcinfo.r = rfine; srcinfo.d = dfine; srcinfo.n = nfine;
srcinfo.d2 = d2fine;
if ~isempty(data)
    srcinfo.data = ddsfine;
end

targinfo = [];
targinfo.r = rt; targinfo.d = dt; targinfo.n = nt;
targinfo.d2 = d2t; targinfo.data = dd;

dfinenrm = sqrt(sum(dfine.^2,1));
dsdt = dfinenrm(:).*whts1(:);

dsdtndim2 = repmat(dsdt(:).',opdims(2),1);
dsdtndim2 = dsdtndim2(:);

smatbig = fkern(srcinfo,targinfo);
submat = smatbig*diag(dsdtndim2)*ainterp1kron;

if nargin >= 14 && corrections
    srcinfo = [];
    srcinfo.r = rs;  srcinfo.d = ds;
    srcinfo.d2 = d2s; srcinfo.n = ns;
    wtsj = wtss(:,j);
    submat = submat - fkern(srcinfo,targinfo).*(wtsj(:).');
end

end

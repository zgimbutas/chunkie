function [ts_aux,ws_aux,xs0,ws0] = getlogquad_aux(k)
%CHNK.QUADGALERKIN.GETLOGQUAD_AUX
%
% Load the Galerkin auxiliary log-singular quadrature tables produced
% by legeexps_log_lr / hqsuppauxquad in the chunkmatc Fortran code
% (file logquad_aux.f). For a given order k the table provides:
%
%   ts_aux,  ws_aux  - 2*k auxiliary target nodes/weights on [-1,1]
%   xs0{j},  ws0{j}  - per-target source rule integrating
%                       log|y - ts_aux(j)| * polynomial(y) on [-1,1]
%                       for j = 1..2*k. Length of each cell is the
%                       table-specific number of source nodes (~ 2*k).
%
% Tables ship for k in {4, 6, 8, 10, 12, 14, 16, 20}.

switch k
    case 4
        [ts_aux,ws_aux,xs0,ws0] = chnk.quadgalerkin.auxlog_nnode004_naux008();
    case 6
        [ts_aux,ws_aux,xs0,ws0] = chnk.quadgalerkin.auxlog_nnode006_naux012();
    case 8
        [ts_aux,ws_aux,xs0,ws0] = chnk.quadgalerkin.auxlog_nnode008_naux016();
    case 10
        [ts_aux,ws_aux,xs0,ws0] = chnk.quadgalerkin.auxlog_nnode010_naux020();
    case 12
        [ts_aux,ws_aux,xs0,ws0] = chnk.quadgalerkin.auxlog_nnode012_naux024();
    case 14
        [ts_aux,ws_aux,xs0,ws0] = chnk.quadgalerkin.auxlog_nnode014_naux028();
    case 16
        [ts_aux,ws_aux,xs0,ws0] = chnk.quadgalerkin.auxlog_nnode016_naux032();
    case 20
        [ts_aux,ws_aux,xs0,ws0] = chnk.quadgalerkin.auxlog_nnode020_naux040();
    otherwise
        error(['chnk.quadgalerkin.getlogquad_aux: no Galerkin auxiliary ' ...
               'table available for order k=%d'],k);
end

end

function [params_out] = adaptMaps(params, eigs_A)
    
    x = params.x;   % <-> r
    y = params.y;   % <-> q
    Lr = params.Lr_scaled;
    Lq = params.Lq_scaled;
    Np = params.Np;
    K = params.K;
    NODETOL = params.NODETOL;
    
    eigs_A_complete = repmat(eigs_A, 1, Np*K);
    eigs_A_complete = eigs_A_complete';
    eigs_A_complete = eigs_A_complete(:);
        
    vmapI_left   = find( eigs_A_complete > 0 && (abs(x + Lr/2) < NODETOL));
    vmapI_right  = find( eigs_A_complete < 0 && (abs(x - Lr/2) < NODETOL));
    params.vmapI = vmapI_left + vmapI_right;
    vmapO_left   = find( eigs_A_complete < 0 && (abs(x + Lr/2) < NODETOL));
    vmapO_right  = find( eigs_A_complete > 0 && (abs(x - Lr/2) < NODETOL));
    params.vmapO = vmapO_left + vmapO_right;
    
    params_out = params;
end
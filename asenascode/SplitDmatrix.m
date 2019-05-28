function [D_l_pos, D_l_neg, D_r_pos, D_r_neg] = SplitDmatrix(params, D, Phi)
    p = params.p_DFT;
    x = params.x;
    Lr = params.Lr_scaled;
    vmapD = params.vmapD;
    NODETOL = params.NODETOL;
    NpTimesK = params.Np * params.K;
    PhiT = Phi';    % within Phi' a single column stands for constant p
                    % within D a single column stands for fixed x and y
        
    map_l = Lr/2 + x(vmapD) < NODETOL;
    map_r = Lr/2 - x(vmapD) < NODETOL;
    map_pos = p >= 0;
    map_neg = p < 0;
    
    D_l = sparse(NpTimesK, NpTimesK);
    D_l_pos = sparse(NpTimesK, NpTimesK);
    D_l_neg = sparse(NpTimesK, NpTimesK);
    D_r_pos = sparse(NpTimesK, NpTimesK);
    D_r_neg = sparse(NpTimesK, NpTimesK);
    
    D_l(:,vmapD(map_l)) = D(:,vmapD(map_l)) ;
    D_l_pos(:, map_pos) = D_l * PhiT(:, map_pos);
    D_l_neg(:, map_neg) = D_l * PhiT(:, map_neg);
    D_r(:,vmapD(map_r)) = D(:,vmapD(map_r)) ;
    D_r_pos(:, map_pos) = D_r * PhiT(:, map_pos);
    D_r_neg(:, map_neg) = D_r * PhiT(:, map_neg);
end
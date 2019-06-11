function [D_l, D_r] = SplitDmatrix(params, D)
% function [D_l_pos, D_l_neg, D_r_pos, D_r_neg] = SplitDmatrix(params, D, R, eigs_A)

    Lr = params.Lr_scaled;
    NODETOL = params.NODETOL;
    Np = params.Np;
    Npy = params.Npy;   % number of nodes
    Ny = Npy - 1;
    K = params.K;
    x = params.x;
    systemsize = K*Np*Ny;
%     RT = R';
    
    map_l = find(Lr/2 + x < NODETOL);
    map_r = find(Lr/2 - x < NODETOL);
    
%     map_pos = find(eigs_A>0);
%     map_neg = find(eigs_A<0);
    
    D_l = sparse(systemsize, systemsize);
    D_r = sparse(systemsize, systemsize);
%     D_l_pos = sparse(systemsize, systemsize);
%     D_l_neg = sparse(systemsize, systemsize);
%     D_r_pos = sparse(systemsize, systemsize);
%     D_r_neg = sparse(systemsize, systemsize);
    
    D_l(:,map_l) = D(:,map_l) ;
%     D_l_pos(:, map_pos) = D_l * RT(:, map_pos);
%     D_l_neg(:, map_neg) = D_l;
    D_r(:,map_r) = D(:,map_r) ;
%     D_r_pos(:, map_pos) = D_r;
%     D_r_neg(:, map_neg) = D_r * RT(:, map_neg);
    
    if(params.testing)
        if (D ~= D_r + D_l)
            error("Something went wrong when splitting D.")
        end
        figure('Name', "Abs of D_l");
        spyc_grid(abs(D_l), 'cool', Np, Np*K);
        figure('Name', "Abs of D_r");
        spyc_grid(abs(D_r), 'cool', Np, Np*K);
%         figure('Name', "Abs of D_l_pos");
%         spyc_grid(abs(D_l_pos), 'cool', Np, Np*K);
%         figure('Name', "Abs of D_r_neg");
%         spyc_grid(abs(D_r_neg), 'cool', Np, Np*K);
    end
end
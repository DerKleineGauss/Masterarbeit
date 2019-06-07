function [D_l, D_r] = SplitDmatrix(params, D)

    Lr = params.Lr_scaled;
    NODETOL = params.NODETOL;
    Np = params.Np;
    Npy = params.Npy;   % number of nodes
    Ny = Npy - 1;
    K = params.K;
    x = params.x;
    systemsize = K*Np*Ny;
    
    map_l = find(Lr/2 + x < NODETOL);
    map_r = find(Lr/2 - x < NODETOL);
    
    D_l = sparse(systemsize, systemsize);
    D_r = sparse(systemsize, systemsize);
    
    D_l(:,map_l) = D(:,map_l) ;
    D_r(:,map_r) = D(:,map_r) ;
    
    if(params.testing)
        figure('Name', "Abs of D_l");
        spyc_grid(abs(D_l), 'cool', Np, Np*K);
        figure('Name', "Abs of D_r");
        spyc_grid(abs(D_r), 'cool', Np, Np*K);
    end
end
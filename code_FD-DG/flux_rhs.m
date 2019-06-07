function [uD] = flux_rhs(params, D)

    Lr = params.Lr_scaled;
    NODETOL = params.NODETOL;
    Np = params.Np;
    Npy = params.Npy;   % number of nodes
    Ny = Npy - 1;
    K = params.K;
    x = params.x;
    y = params.y;
    systemsize = K*Np*Ny;

    [a, b, c] = params.get_abc(params.gamma, params.epsilon);
    map_left_v = Lr/2 + x < NODETOL;
    map_right_v = Lr/2 - x < NODETOL;
    fh_left = @(k) a.*cos(y(map_left_v )*k).*log(1+exp(-b.*k^2+c));
    fh_right= @(k) a.*cos(y(map_right_v)*k).*log(1+exp(-b.*k^2+c));
    
    uD = zeros(systemsize, 1);
    uD(map_left_v) = integral(fh_left,-2*c,2*c,'ArrayValued', true);
    uD(map_right_v) = integral(fh_right,-2*c,2*c,'ArrayValued', true);
    
    uD = D*uD;
    
end
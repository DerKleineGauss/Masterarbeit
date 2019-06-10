function [rhs] = flux_rhs(params, D_l, D_r, R, eigs_A)

    Lr = params.Lr_scaled;
    NODETOL = params.NODETOL;
    Np = params.Np;
    Npy = params.Npy;   % number of nodes
    Ny = Npy - 1;
    K = params.K;
    x = params.x;
    y = params.y;
    systemsize = K*Np*Ny;
    RT = R';

    [a, b, c] = params.get_abc(params.gamma, params.epsilon);
    map_left_v = Lr/2 + x < NODETOL;
    map_right_v = Lr/2 - x < NODETOL;
    fh_left = @(k) a.*cos(y(map_left_v )*k).*log(1+exp(-b.*k^2+c));
    fh_right= @(k) a.*cos(y(map_right_v)*k).*log(1+exp(-b.*k^2+c));
    
    f = zeros(systemsize, 1);
    f(map_left_v) = integral(fh_left,-2*c,2*c,'ArrayValued', true);
    f(map_right_v) = integral(fh_right,-2*c,2*c,'ArrayValued', true);
    f = reshape(f, Np*K, Ny);
    
    map_pos = find(eigs_A>0);
    map_neg = find(eigs_A<0);
    
    RT_pos = 0*RT; RT_neg = 0*RT;
%     RT_pos(:, map_pos) = RT(:, map_pos);
%     RT_neg(:, map_neg) = RT(:, map_neg);
    RT_pos( map_pos, :) = RT( map_pos, :);
    RT_neg( map_neg, :) = RT( map_neg, :);
    rhs = D_l*strangeMatrixMultiplication(RT_pos , f) + D_r*strangeMatrixMultiplication(RT_neg , f);
    
end
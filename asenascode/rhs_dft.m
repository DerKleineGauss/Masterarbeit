function [uD] = rhs_dft(params, a,b,c, Phi, D_l_pos, D_r_neg, D_l_neg, D_r_pos)

    Lr = params.Lr_scaled;
    NODETOL = params.NODETOL;
    Np = params.Np;
    K = params.K;
    x = params.x;
    y = params.y;

    map_left_v = Lr/2 + x < NODETOL;
    map_right_v = Lr/2 - x < NODETOL;
    fh_left = @(k) a.*cos(y(map_left_v )*k).*log(1+exp(-b.*k^2+c));
    fh_right= @(k) a.*cos(y(map_right_v)*k).*log(1+exp(-b.*k^2+c));
    % uD_v caontains values in "volume speak" generated with x and y instead Fx
    % and Fy (as is the case for the later defined uD_complete)
    uD = zeros(Np*K, 1);
    uD(map_left_v) = integral(fh_left,-2*c,2*c,'ArrayValued', true);
    uD(map_right_v) = integral(fh_right,-2*c,2*c,'ArrayValued', true);
    % left and right sides must have same y values -> check here
    if (params.testing)
        assert ( norm(uD(map_left_v) - uD(map_right_v)) / norm(uD(map_left_v)) < 1e-16)
        figure(111);
        plot(y(map_left_v), uD(map_left_v), '.');
        % confirm we do not lose map control
        assert(norm(find(Phi*uD)-find(uD)) < NODETOL)
    end
    % apply DFT to rhs
    uD = (D_l_pos + D_r_neg + D_l_neg + D_r_pos)*uD;

return;

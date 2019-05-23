function [uD] = rhs_dft(params, a,b,c,Phi)

    Lr = params.Lr_scaled;
    NODETOL = params.NODETOL;
    Nfaces = params.Nfaces;
    Np = params.Np;
    Nfp = params.Nfp;
    K = params.K;
    x = params.x;
    y = params.y;
    Fx = params.Fx;
    P_DFT = params.P_DFT;
    mapD = params.mapD;
    Fmask = params.Fmask;


    map_left_v = Lr/2 + x < NODETOL;
    map_right_v = Lr/2 - x < NODETOL;
    fh_left = @(k) a.*cos(y(map_left_v )*k).*log(1+exp(-b.*k^2+c));
    fh_right= @(k) a.*cos(y(map_right_v)*k).*log(1+exp(-b.*k^2+c));
    % uD_v caontains values in "volume speak" generated with x and y instead Fx
    % and Fy (as is the case for the later defined uD_complete)
    uD_v = zeros(Np, K);
    uD_v(map_left_v) = integral(fh_left,-2*c,2*c,'ArrayValued', true);
    uD_v(map_right_v) = integral(fh_right,-2*c,2*c,'ArrayValued', true);
    % left and right sides must have same y values -> check here
    if (params.testing)
        assert ( norm(uD_v(map_left_v) - uD_v(map_right_v)) / norm(uD_v(map_left_v)) < 1e-16)
        figure(111);
        plot_sorted(y(map_left_v), uD_v(map_left_v));
        % confirm we do not lose map control
        assert(norm(find(Phi*uD_v(:))-find(uD_v)) < NODETOL)
    end
    % apply DFT to rhs
    uD_v = reshape(Phi*uD_v(:), Np, K);
    % go to needed definitions for each face as also done in StartUp2D for Fx
    uD_complete = uD_v(Fmask(:), :);
    uD = zeros(Nfp*Nfaces, K);
    bedleft = Lr/2+(Fx(mapD)) < NODETOL & P_DFT(mapD) >= 0;
    bedright = Lr/2-(Fx(mapD)) < NODETOL & P_DFT(mapD) <= 0;
    if (params.testing)
        figure(113)
        plot(Fx(:),P_DFT(:),'k.', Fx(mapD(bedleft)),P_DFT(mapD(bedleft)),'gx', Fx(mapD(bedright)),P_DFT(mapD(bedright)),'gx');
        legend('all face points', {['Dirichlet face points to'];['be set with fermi-dirac']}, 'Location','north')
    end
    uD(mapD(bedleft)) = uD_complete(mapD(bedleft));
    uD(mapD(bedright))= uD_complete(mapD(bedright));
return;

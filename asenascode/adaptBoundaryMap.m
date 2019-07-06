function [params_out] = adaptBoundaryMap(params)

  Lr = params.Lr_scaled;
  VX = params.VX;
  VY = params.VY;
  NODETOL = params.NODETOL;


    right_pos = find(Lr/2-VX < NODETOL & VY <= 0);  % negative y corresponds to positive p_DFT
    left_pos  = find(Lr/2+VX < NODETOL & VY <= 0);
    right_neg = find(Lr/2-VX < NODETOL & VY >= 0);
    left_neg  = find(Lr/2+VX < NODETOL & VY >= 0);
    [params] = CorrectBCTable_rectangular(params, left_neg, params.Neuman);
    [params] = CorrectBCTable_rectangular(params, right_pos, params.Neuman);
    [params] = BuildBCMaps2D(params);

    if (params.testing)
        P_DFT = params.P_DFT;
        mapD = params.mapD;
        mapN = params.mapN;
        Fx = params.Fx;
        figure(109)
        plot(VX,VY,'k.',VX(right_pos), VY(right_pos), 'ro', VX(left_pos), VY(left_pos), 'bo', ...
            VX(right_neg), VY(right_neg), 'yx', VX(left_neg), VY(left_neg), 'gx');
        figure(110)
        plot(Fx(:),P_DFT(:),'k.',Fx(mapN),P_DFT(mapN),'ro', Fx(mapD),P_DFT(mapD),'gx');
        legend('all face points','Neumann face points', 'Dirichlet face points', 'Location','north')
    end
params_out = params;
return;

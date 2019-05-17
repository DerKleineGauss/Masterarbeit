function adaptBoundaryMap(Lr)
    Globals2D;
    right_pos = find(Lr/2-VX < NODETOL & VY <= 0);  % negative y corresponds to positive p_DFT
    left_pos  = find(Lr/2+VX < NODETOL & VY <= 0);
    right_neg = find(Lr/2-VX < NODETOL & VY >= 0);
    left_neg  = find(Lr/2+VX < NODETOL & VY >= 0);
    CorrectBCTable(left_neg, Neuman);
    CorrectBCTable(right_pos, Neuman);
    BuildBCMaps2D;
    
    if (testing)
        figure(109)
        plot(VX,VY,'k.',VX(right_pos), VY(right_pos), 'ro', VX(left_pos), VY(left_pos), 'bo', ...
            VX(right_neg), VY(right_neg), 'yx', VX(left_neg), VY(left_neg), 'gx');
        figure(110)
        plot(Fx(:),P_DFT(:),'k.',Fx(mapN),P_DFT(mapN),'ro', Fx(mapD),P_DFT(mapD),'gx');
        legend('all face points','Neumann face points', 'Dirichlet face points', 'Location','north')
    end
return;
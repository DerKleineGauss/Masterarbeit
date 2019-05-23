function [out_params] = CorrectBCTable_rectangular(params, mapnodes,BCcode)

% function BCType = CorrectBCTable(EToV,BCType,mapnodes,BCcode);
% Purpose: Setup BCType for boundary conditions in 2D
% By Allan P. Engsig-Karup

VNUM = [1 2;2 3;3 4;4 1]; % face orientations

for k = 1:params.K
    % Test for each edge
    for f = 1:params.Nfaces
        m = params.EToV(k,VNUM(f,1)); n = params.EToV(k,VNUM(f,2));

        % if both points are on the boundary then it is a boundary point!
        ok=sum(ismember([m n],mapnodes));
        if ok==2
            params.BCType(k,f)=BCcode;
        end

    end
end
out_params = params;
return

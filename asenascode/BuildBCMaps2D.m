function [params_out] = BuildBCMaps2D(params)

  % function BuildMaps2DBC
  % Purpose: Build specialized nodal maps for various types of
  %          boundary conditions, specified in BCType.

  Nfp = params.Nfp;
  BCType = params.BCType;

  % create label of face nodes with boundary types from BCType
  bct    = BCType';
  bnodes = ones(Nfp, 1)*bct(:)';
  bnodes = bnodes(:);

  % find location of boundary nodes in face and volume node lists
  params.mapI = find(bnodes==params.In);           params.vmapI = params.vmapM(params.mapI);
  params.mapO = find(bnodes==params.Out);          params.vmapO = params.vmapM(params.mapO);
  params.mapW = find(bnodes==params.Wall);         params.vmapW = params.vmapM(params.mapW);
  params.mapF = find(bnodes==params.Far);          params.vmapF = params.vmapM(params.mapF);
  params.mapC = find(bnodes==params.Cyl);          params.vmapC = params.vmapM(params.mapC);
  params.mapD = find(bnodes==params.Dirichlet);    params.vmapD = params.vmapM(params.mapD);
  params.mapN = find(bnodes==params.Neuman);       params.vmapN = params.vmapM(params.mapN);
  params.mapS = find(bnodes==params.Slip);         params.vmapS = params.vmapM(params.mapS);
  params_out = params;
return;

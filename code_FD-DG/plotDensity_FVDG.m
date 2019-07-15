function [] = plotDensity_FVDG(v, params, R, n_linspace, fid, timelocal, color)

% function [] = plotDensity_FVDG(v, params, R, n_interpol)
% Purpose : Plot density function of finite volume / dg solution v which is
%           in eigenspace
%           R :         matrix of eigenvectors
%           n_linspace: size of linspace where v shall be evaluated
%           fid:        id of figure to get with fid =figure;

Lx = params.Lr_scaled;
Kx = params.K;
plot_logarithmic = params.plot_logarithmic;

Np = params.Np;
Ny = params.Ny;       % equals 'N_q' in Lukas Paper : the number of cells

v = reshape(v,params.Np*params.K, params.Ny);
u = strangeMatrixMultiplication(R , v);
u = reshape(u,params.Np*params.K, params.Ny);

% y = 0 is always the interface between cell Ny/2 and cell Ny/2+1
u_0 = (u(:,Ny/2)+u(:,Ny/2 + 1)) / 2;

xlin= linspace(-Lx/2,Lx/2,n_linspace);

rhoint= 0*xlin;

hx= params.VX(2)-params.VX(1);

cellx= int64((xlin+Lx/2-mod(xlin+Lx/2,hx))/hx+1);
cellx(end)= cellx(end)-1;

offx= 0;
for j=1:Kx
    indx= find(cellx==j);
    if ~isempty(indx)
        xloc= xlin(indx);
        r= 2*mod(xloc+Lx/2,hx)/hx-1;
        uc= u_0((offx+1):(offx+Np));
        Vsx= Vandermonde1D(params, r(:));
        invV = params.invV;
        rhoint(indx)= uc'*invV'*Vsx';
    end

    offx= offx+Np;
end

figure(fid)
% clf(fid)
yyaxis left
if (plot_logarithmic)
    semilogy(xlin,real(rhoint),'.', 'color', color)
    semilogy(xlin,real(rhoint), 'color', color)
else
    ylim([-1e23 6e23])
    plot(xlin,real(rhoint),'.', 'color', color)
    plot(xlin,real(rhoint), 'color', color)
end
hold on
yyaxis right
plot(xlin, Potential(xlin, timelocal, params))
% hold off

end

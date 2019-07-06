function [] = plotSolution(Kx,Ky,rho, param, VX, VY,invVx, invVy)

xlin= linspace(-param.Lx/2,param.Lx/2,Kx);
ylin= linspace(-param.Ly/2,param.Ly/2,Ky);
[X,Y]= meshgrid(xlin,ylin);
rhoint= zeros(Ky,Kx);

hx= abs(VX(2)-VX(1));
hy= abs(VY(2)-VY(1));



cellx= int64((xlin+param.Lx/2-mod(xlin+param.Lx/2,hx))/hx+1);
cellx(end)= cellx(end)-1;
celly= int64((ylin+param.Ly/2-mod(ylin+param.Ly/2,hy))/hy+1);
celly(end)= celly(end)-1;

offx= 0;
for j=1:param.Kx
    indx= find(cellx==j);
    offy=0;
    for k=1:param.Ky
        indy = find(celly==k);
        if ~isempty(indx) && ~isempty(indy)
            xloc= xlin(indx);
            yloc= ylin(indy);
            r= 2*mod(xloc+param.Lx/2,hx)/hx-1;
            s= 2*mod(yloc+param.Ly/2,hy)/hy-1;
            rhoc= rho((offy+1):(offy+param.Npy), (offx+1):(offx+param.Npx));
            Vsy= Vandermonde1D(param.Ny,s(:));
            Vsx= Vandermonde1D(param.Nx,r(:));
            rhoint(indy, indx)= Vsy*invVy*rhoc*invVx'*Vsx';
        end
        offy= offy+param.Ny;
    end
    offx= offx+param.Npx;
end

mesh(X,Y,real(rhoint))

end

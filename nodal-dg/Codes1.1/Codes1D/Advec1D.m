function [u] = Advec1D(u, FinalTime)

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;
time = 0;

% Runge-Kutta residual storage  
resu = zeros(Np,K); 

% advection speed
a = 2*pi;
a = a*10;

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.75; dt   = CFL/(a)*xmin; dt = .5*dt;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 

line_handle = plot(reshape(u,[K*Np, 1]));
axis([0 2 -1 1]);
% outer time step loop 
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        [rhsu] = AdvecRHS1D(u, timelocal, a);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u+rk4b(INTRK)*resu;
    end;
    % Increment time
    time = time+dt;
    set ( line_handle, 'Xdata', reshape (x(2:Np ,:), [(Np-1)*K,1]), 'Ydata', reshape (u(2:Np ,:), [(Np-1)*K,1]));
end;
return

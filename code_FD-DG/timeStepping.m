function [v] = timeStepping(params, v, moviename)


rk4a = params.rk4a; rk4b = params.rk4b; rk4c = params.rk4c;
time = 0;

% Runge-Kutta residual storage
resv = zeros(params.Np*params.K* params.Ny, 1);

% get time-independet parts of systemmatrix
[L_glob, rhs, R, maxSpeed] = LVN_systemmatrix_timeIndependant(params);

% compute time step size
CFL=0.75; dt   = CFL/(maxSpeed*params.hy)*params.hx; 
Nsteps = ceil(params.FinalTime/dt); dt = params.FinalTime/Nsteps;

fid = figure( 'name', "Dichte");
xlabel('$x / \xi$')
ylabel('$n(x)$')

if (params.makeMovie)
    % Set up the movie.
    writerObj = VideoWriter(moviename); % Name it.
    writerObj.FrameRate = 10; % How many frames per second.
    open(writerObj); 
end

% outer time step loop
gotLastG_glob = false;
for tstep=1:Nsteps
    if (mod(tstep,10)==0)
        display(['Timestep ',num2str(tstep),'/',num2str(Nsteps)])
    end
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        
        if (timelocal < params.rampTime || params.rampTime <= 0)    % only change G if potential changes
            % set potential
            B= functionB(params, timelocal);
            if(params.testing)
                figure(4)
                plot3(params.x_interface(:),params.y_interface(:), real(B(:)),'x')
                figure(5)
                plot3(params.x_interface(:),params.y_interface(:), imag(B(:)),'x')
            end
            [G_glob] = LVN_systemmatrix_timeDependant(params, B, R);
        elseif ( timelocal >= params.rampTime && ~gotLastG_glob)
            B= functionB(params, timelocal);
            G_glob = LVN_systemmatrix_timeDependant(params, B, R);
            gotLastG_glob = true;
        end
            
        rhsv = (-L_glob - G_glob)*v + rhs;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        v = v+rk4b(INTRK)*resv;
    end
    % Increment time
    time = time+dt;
    if mod(tstep,2)==0 && params.makeMovie
        n_linspace = 300;
        plotDensity_FVDG(v, params, R, n_linspace, fid, timelocal);
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
    end
end
if(params.makeMovie)
    close(writerObj);
end

end

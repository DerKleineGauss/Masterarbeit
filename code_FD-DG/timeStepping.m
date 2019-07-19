function [v] = timeStepping(params, v, v_final_stationary)


rk4a = params.rk4a; rk4b = params.rk4b; rk4c = params.rk4c;
time = 0;
timetick_counter = 1;

% Runge-Kutta residual storage
resv = zeros(params.Np*params.K* params.Ny, 1);

% get time-independet parts of systemmatrix
[L_glob, rhs, R, maxSpeed] = LVN_systemmatrix_timeIndependant(params);

% compute time step size
dt   = 1/(maxSpeed*params.hy)*params.hx / (2*params.N+1);  
dt = 2*dt;
if (strcmp(params.mode,'rk2ssp'))
    dt = dt*0.25;
end
Nsteps = ceil(params.FinalTime/dt); dt = params.FinalTime/Nsteps;

timeticks_in_fs = [0, 50, 200, 800, 16000];
timeticks_in_characteristic = timeticks_in_fs*1e-15 * params.epsilon;
timeticks_in_tsteps = ceil(timeticks_in_characteristic/dt)+1;


fid = figure( 'name', "Dichte");
xlabel('$x / \xi$')
ylabel('$n(x)$')

if (params.makeMovie)
    % Set up the movie.
    writerObj = VideoWriter(params.movieName); % Name it.
    writerObj.FrameRate = 10; % How many frames per second.
    open(writerObj);
end

% outer time step loop
gotLastG_glob = false;
for tstep=1:Nsteps
    if (mod(tstep,10)==0)
        display(['Timestep ',num2str(tstep),'/',num2str(Nsteps)])
    end
    if (strcmp(params.mode,'rk4'))
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
    elseif (strcmp(params.mode,'rk2ssp'))
        v_temp = v;

        for INTRK = 1:2
            timelocal = time + dt;

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

            rhsv = (-L_glob - G_glob)*v_temp + rhs;
            v_temp = v_temp + dt*rhsv;
        end
        v = 0.5*(v_temp + v);
    end
    % Increment time
    time = time+dt;
    if mod((tstep-1),1)==0 && params.makeMovie
        n_linspace = 300;
        clf(fid);
%         if (params.plot_logarithmic)
%             semilogy(1:10,1:10);
%         end

        plotDensity_FVDG(v, params, R, n_linspace, fid, timelocal, 'b');
        hold on
        plotDensity_FVDG(v_final_stationary, params, R, n_linspace, fid, timelocal, 'r');
        hold off
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
    end
    if (tstep == timeticks_in_tsteps(timetick_counter) && params.saveResults)
        save(['results/intermediate_t', num2str(timeticks_in_fs(timetick_counter)), 'fs'],'v');
        savefig(['results/intermediate_t', num2str(timeticks_in_fs(timetick_counter)), 'fs_plot'])
        timetick_counter = timetick_counter+1;
    end
end
if(params.makeMovie)
    close(writerObj);
end

end

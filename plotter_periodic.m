% Script 'plotter_periodic.m' v1.0.0, tested 28 September 2024
% Written by P.J. Shamberger (c) 2023

%% PLOTTER_PERIODIC  ------------------------------------------------------
% Script to generate plots illustrating transient thermal response.

%close all
%clear all

i = 1; %1e7 W/m2
j = 1;  

% --- FLAGS ---------------------------------------------------------------
% 1 = on, 0 = off
time_flag = 1;      % this is good to look at the output for a single set of conditions.
impedence_flag = 1; % this makes standard impedence plots


% #########################################################################
%% --- Time-dependent plots  ----------------------------------------------
% Note: if "wrapper_singlepulse" was run for a range of pulse lengths
% (t_on's), these plots will illustrate the results for the last condition
% run.

if time_flag
    
    % plot T_hot, q" vs time
    figure(1)
        [ax,h1,h2] = plotyy(time_out(:, 1, i, j), qpulse_out(:, 1, i, j),time_out(:, 1, i, j), T_hot_out(:, 1, i, j));
        hold on
        legend('q"/W*m^{-2}', 'T_{hot}/C')
        title('Transient thermal response during pulse')
        xlabel('time / s')
        ylabel('q", T_{hot}')
    
    % plot T_hot, T_cold, T_slab (T at 3 different defined positions) vs time
    figure(2)
        plot(time_out(:, 1, i, j), T_hot_out(:, 1, i, j),'r.');
        hold on
        plot(time_out(:, 1, i, j), T_cold_out(:, 1, i, j),'b.');
        plot(time_out(:, 1, i, j), T_slab_out(:, 1, i, j),'k.');
        title('Slab temperatures during heat pulse')
        legend('T_{hot}/C', 'T_{cold}/C', 'T_{slab}/C')
        xlabel('time / s')
        ylabel('T / C')
        
    % Track instantaneous heat asborption vs time (partitioned into
    % sensible, latent, and total heat)
    figure(3)
        plot(time1(4:length(time1)), qsense(4:length(time1)));
        hold on
        plot(time1(4:length(time1)), qlatent(4:length(time1)));
        plot(time1(4:length(time1)), qlatent(4:length(time1))+qsense(4:length(time1)));
        title('Track instantaneous heat asborption vs time')
        legend('q_{sense}', 'q_{latent}', 'q_{tot}')
        xlabel('time [s]')
        ylabel('q" [W/m2]')

    % plot a few example temperature profiles through the thickness of the PCM layer
    dimensions = size(T2);

    figure(4)
        plot(fliplr((1:dimensions(1)).*dely), T2(:,dimensions(2)),'Color',[1 0 0])
        hold on
        plot(fliplr((1:dimensions(1)).*dely), T2(:,dimensions(2)-floor(pts_per_cycle/10)+1),'Color',[0.9 0 0.1])
        plot(fliplr((1:dimensions(1)).*dely), T2(:,dimensions(2)-floor(pts_per_cycle/10*2)+1),'Color',[0.8 0 0.2])
        plot(fliplr((1:dimensions(1)).*dely), T2(:,dimensions(2)-floor(pts_per_cycle/10*3)+1),'Color',[0.7 0 0.3])
        plot(fliplr((1:dimensions(1)).*dely), T2(:,dimensions(2)-floor(pts_per_cycle/10*4)+1),'Color',[0.6 0 0.4])
        plot(fliplr((1:dimensions(1)).*dely), T2(:,dimensions(2)-floor(pts_per_cycle/10*5)+1),'Color',[0.5 0 0.5])
        plot(fliplr((1:dimensions(1)).*dely), T2(:,dimensions(2)-floor(pts_per_cycle/10*6)+1),'Color',[0.4 0 0.6])
        plot(fliplr((1:dimensions(1)).*dely), T2(:,dimensions(2)-floor(pts_per_cycle/10*7)+1),'Color',[0.3 0 0.7])
        plot(fliplr((1:dimensions(1)).*dely), T2(:,dimensions(2)-floor(pts_per_cycle/10*8)+1),'Color',[0.2 0 0.8])
        plot(fliplr((1:dimensions(1)).*dely), T2(:,dimensions(2)-floor(pts_per_cycle/10*9)+1),'Color',[0.1 0 0.9])
        %plot(fliplr((1:dimensions(1)).*dely), T2(:,dimensions(2)-floor(pts_per_cycle)+1),'Color',[0.0 0 1])
        title('Example internal temperature profiles')
        xlabel('distance / m')
        ylabel('Temperature / C')

end

% #########################################################################
%% --- Impedence plots  ---------------------------------------------------

if impedence_flag
    if ~BC_flag   % Periodic q BC
        for i = 1:length(df_vect) % duty factor index
            for j = 1:length(t1_vect) % pulse period index
               
                % Calculates the temperature rise of the junction over the cooling temperature
                % max_T_hot_mx(iii,jjj) = max_T_hot_out - boundary.T_cooling_hs;
                
            end
            
            r_index = (i-1)./(length(df_vect));

            if r_index == NaN
                r_index = 0;
            end

            % plots the temperature rise at different t_on / DF values
            figure(100)
                loglog(t1_vect, max_T_hot_mx(i,:),'x', 'color', [r_index 0 0])
                hold on

                title('Temperature rise')
                ylabel('max \DeltaT_j / C')
                xlabel('log_{10} time / s')
                %colorbar

            % Total power = q" x area = geometry.area x boundary.q
            power = geometry.area.*boundary.q;

            % plots the thermal impedence at different t_on / DF values
            figure(101)
                loglog(t1_vect, max_T_hot_mx(i,:)./power,'x', 'color', [r_index 0 0])
                hold on

                title('Thermal Impedance')
                ylabel('max /Delta T_j / Power / C/W')
                xlabel('log_{10} time / s')
                %colorbar
            
        end  
    end
        
end
    
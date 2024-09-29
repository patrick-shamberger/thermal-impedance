% Script 'plotter_impedance.m' v1.0.0, tested 28 September 2024
% Written by P.J. Shamberger (c) 2023

%% PLOTTER_IMPEDANCE  -----------------------------------------------------
% Script to generate plots illustrating transient thermal response.

% boundary.q = [1E5]; %[W/m2]
% geometry.area = 0.0001;     %[m2] = 1 cm2.  ONLY used for power normalization.

[dim1, dim2] = size(max_T_hot_mx);

% #########################################################################
%% --- Impedence plots  ---------------------------------------------------
for i_plot = 1:dim1 % duty factor index

    r_index = (i_plot-1)./(dim1);

    if r_index == NaN
        r_index = 0;
    end

    % plots the TEMPERATURE RISE at different t_on / DF values
    figure(100)
        loglog(t1_vect, max_T_hot_mx(i_plot,:),'x', 'color', [r_index 0 0])
        hold on
        loglog(t1_vect, max_T_hot_mx(i_plot,:),'-', 'color', [r_index 0 0])
        
        ylabel('max \DeltaT_j / C')
        xlabel('log_{10} time / s')
        %colorbar

    % Total power = q" x area = geometry.area x boundary.q
    power = geometry.area.*boundary.q;

    % plots the THERMAL IMPEDENCE at different t_on / DF values
    figure(101)
        loglog(t1_vect, max_T_hot_mx(i_plot,:)./power,'x', 'color', [r_index 0 0])
        hold on
        loglog(t1_vect, max_T_hot_mx(i_plot,:)./power,'-', 'color', [r_index 0 0])
        
        ylabel('max \DeltaT_j / Power / C/W')
        xlabel('log_{10} time / s')
        
        %xlim([1e-3 1e3])
        %ylim([3e-3 1e1])
        %colorbar

    % plots the UTILIZATION at different t_on / DF values
    figure(102)
        semilogx(t1_vect, utilization_mx(i_plot,:),'x', 'color', [r_index 0 0])
        hold on
        semilogx(t1_vect, utilization_mx(i_plot,:),'-', 'color', [r_index 0 0])
        
        ylabel('utilization')
        xlabel('log_{10} time / s')
        
        %xlim([1e-5 1e3])
        %ylim([0 1])
        %colorbar

    % plots the STORAGE FRACTION at different t_on / DF values
    figure(103)
        semilogx(t1_vect, storage_fraction_mx(i_plot,:),'x', 'color', [r_index 0 0])
        hold on
        semilogx(t1_vect, storage_fraction_mx(i_plot,:),'-', 'color', [r_index 0 0])
        
        ylabel('storage fraction')
        xlabel('log_{10} time / s')
        
        %xlim([1e-5 1e3])
        %ylim([0 1])
        %colorbar
end  
    
% #########################################################################

load('test_singlepulse.mat'); 

% #########################################################################
%% --- Impedence plots  ---------------------------------------------------
% Note: this section adds the DF = 0 (single pulse) line to the plots
% Note: this section will only work if 'wrapper_singlepulse.m' is run first

% plots the TEMPERATURE RISE at different t_on values / DF=0
figure(100)
    loglog(t1_vect, max_T_hot_mx(1,:),'b.')
    hold on
    loglog(t1_vect, max_T_hot_mx(1,:),'b-')

% plots the THERMAL IMPEDENCE at different t_on values / DF=0
figure(101)
    loglog(t1_vect, max_T_hot_mx(1,:)./power,'b.')
    hold on
    loglog(t1_vect, max_T_hot_mx(1,:)./power,'b-')

% plots the UTILIZATION at different t_on values / DF=0
figure(102)
    semilogx(t1_vect, utilization_mx(i_plot,:), 'bx')
    hold on
    semilogx(t1_vect, utilization_mx(i_plot,:), 'b-')
    
% plots the STORAGE FRACTION at different t_on values / DF=0
figure(103)
    semilogx(t1_vect, storage_fraction_mx(i_plot,:), 'bx')
    hold on
    semilogx(t1_vect, storage_fraction_mx(i_plot,:), 'b-')


% #########################################################################
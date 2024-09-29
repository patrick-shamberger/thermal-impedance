% Script 'wrapper_pulsetrain.m' v1.0.0, tested 28 September 2024
% Written by P.J. Shamberger (c) 2023

%% WRAPPER_PULSETRAIN  ----------------------------------------------------
% 'wrapper_pulsetrain.m' is a matlab script that calculates the thermal
% impedance of a slab of PCM in series configuration between a pulsed heat
% source and an (isothermal or convective) heat sink. This is run for a
% series of pulses, defined by a pulse length (t_on) and a duty factory
% (DF).
%
% Required functions / scripts / databases:
% 1) 'Material Properties.mat': a database that includes material
% properties for the PCM in question.
%
% 2) 'pulsed_systemdef.m': a script that defines the system geometry,
% boundary conditions, etc.
%
% 3a) 'single_pulse_1D.m': a script that calculates (by a finite difference
% computation scheme) the transient thermal response to a step change in
% the thermal boundary condition.
%
% 3b) 'pulsed_1D_time.m': a script that calculates (by a finite difference
% computation scheme) the transient thermal response to a pulsed boundary 
% condition.
%
% 4) 'plotter_periodic.m': a script to generate typical plots from output
% variables
%
% 

%% Initial Setup
clear all
clc
close all

warning('off','all')


%% step 1: load material properties as matrix
load(fullfile('Material Properties.mat')); 

%% step 2: load a description of the system to investigate
pulsed_systemdef;

%% step 3: Setup Tm and q" for this particular run
% done within pulsed_systemdef

%% step 4: Setup duty factor (DF) and pulse time (t1) values to loop over

% these are generally setup to loop through:
%   1) duty factor (stored in df_vect)
df_vect = [0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01];
%df_vect = [0.1, 0.01];
%df_vect = [0.5];

%   2) pulse time (stored in t1_vect)
%t1_vect = 1;
t1_vect = logspace(-3, 3, 7);   % evenly spaced from 0.001 to 1000 s

%% step 5: Run loops
max_value = length(df_vect).*length(t1_vect);
counter_wrapper = 0;
hhh = waitbar(counter_wrapper,'Running DF, t1 loops');

%% Tm = 18
materials.Tm = 18;

for iii = 1:length(df_vect) % duty factor index
    for jjj = 1:length(t1_vect) % pulse period index
        
        counter_wrapper = counter_wrapper +1;
        waitbar(counter_wrapper./max_value,hhh,'Running DF, t1 loops');

        fprintf('Duty Factor: %f \n', df_vect(iii));
        fprintf('Pulse Time: %f s\n', t1_vect(jjj));

        % step 3a: run to find an initial temperature distribution (given a power, and DF)
        % Note: this is helpful for short pulses.  The aim is to get the background
        % temperature distribution based on the *average* power to help the
        % convergence time.

        % runs a single pulse (DF = 0 limit!), at average power to attain quasisteady-state temperature distribution
        boundary.q = boundary.q.*df_vect(iii);   % average power
        boundary.df = 1;
        %boundary.t1 = 100;
        
        % #################################################################
        %   Original
        %boundary.t1 = 100;       
        
        %   Modified to work well for low t_on cases
        boundary.t1 = 1e9;
        timing.pts_per_cycle = 5000;     % # of sampling points/cycle
        % #################################################################
        
        tic
        single_pulse_1D   % Runs a *long* single pulse at the *average* power to get an initial temperature condition.
        toc
        
        dimensions = size(T2);  % Saves the output T distribution
        T_input = T2(:,dimensions(2));
        
%        plot(T_input)
%        stop
        
        % step 3b: run a pulse train experiment (define t_pulse, DF)
        pulsed_systemdef;       % reload system description
        materials.Tm = 18;

        % for each instance in the loop:
        %   1) define df and t1 for that instance
        boundary.df = df_vect(iii);   % duty factor.  df = t1 / (t1 + t2);
        boundary.t1 = t1_vect(jjj); % s.  Rectangular 'ON' pulse duration      
        
        %   2) call the script
        input_Temp_flag = 1;    % this uses the initial Temperature distribution from the steady-state pulse
        
        tic
        pulsed_1D_time;
        toc
        fprintf('\n');

        input_Temp_flag = 0;
        
        %   3) store whatever value you want to save
        max_T_hot_mx(iii,jjj) = max_T_hot_out - boundary.T_cooling_hs;
        cycle_latent_mx(iii,jjj) = cycle_latent_out;
        utilization_mx(iii,jjj) = utilization_out;
        storage_fraction_mx(iii,jjj) = storage_fraction_out;
    
        savefile = 'ongoing_save_Tm_18.mat';
        save(savefile, 'max_T_hot_mx','cycle_latent_mx','utilization_mx','storage_fraction_mx','boundary', 'geometry', 'df_vect', 't1_vect');
        
    end
end


%% step 6: Plot the results
plotter_periodic;

close(hhh)
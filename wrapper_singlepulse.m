% Script 'wrapper_singlepulse.m' v1.0.0, tested 28 September 2024
% Written by P.J. Shamberger (c) 2023

%% WRAPPER_SINGLEPULSE  ---------------------------------------------------
% 'wrapper_singlepulse.m' is a matlab script that calculates the thermal
% impedance of a slab of PCM in series configuration between a pulsed heat
% source and an (isothermal or convective) heat sink.  This is run for a
% single individual pulse (i.e., Duty Factor of 0 - an infinite period of
% time between pulses)
%
% Required functions / scripts / databases:
% 1) 'Material Properties.mat': a database that includes material
% properties for the PCM in question.
%
% 2) 'pulsed_systemdef.m': a script that defines the system geometry,
% boundary conditions, etc.
%
% 3) 'single_pulse_1D.m': a script that calculates (by a finite difference
% computation scheme) the transient thermal response to a step change in
% the thermal boundary condition.
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
savefile_singlepulse = 'test_singlepulse.mat';

%% step 1: load material properties as matrix
load(fullfile('Material Properties.mat')); 

%% step 2: load a description of the system to investigate
pulsed_systemdef;
 
%% step 3: Setup Tm and q" for this particular run
% done within pulsed_systemdef

%% step 4: run a single pulse experiment (define t_pulse, leave DF fixed)

% these are generally setup to loop through:
%   1) duty factor (stored in df_vect)
%   Note: for "single pulse" - this should be left as [1].
df_vect = [1];


%   2) pulse time (stored in t1_vect)
%t1_vect = [1e-3];    % evenly spaced from 0.001 to 1000 s
t1_vect = logspace(-3, 3, 7);

max_value = length(df_vect).*length(t1_vect);
counter_wrapper = 0;
hhh = waitbar(counter_wrapper,'Running DF, t1 loops');

for iii = 1:length(df_vect) % duty factor index
    for jjj = 1:length(t1_vect) % pulse period index
        counter_wrapper = counter_wrapper +1;
        waitbar(counter_wrapper./max_value,hhh,'Running DF, t1 loops');
        
        % for each instance in the loop:
        %   1) define df and t1 for that instance
        boundary.df = df_vect(iii);   % duty factor.  df = t1 / (t1 + t2);
        boundary.t1 = t1_vect(jjj); % s.  Rectangular 'ON' pulse duration      
        
        %   1) call the script
        input_Temp_flag = 0;    % this means that the initial Temperature distribution is "unknown"
        
        tic
        single_pulse_1D;
        toc
        
        %   2) store whatever value you want to save
        max_T_hot_mx(iii,jjj) = max_T_hot_out - boundary.T_cooling_hs;
        save(savefile_singlepulse, 'max_T_hot_mx');

    end
end

% 3) then plot the output
plotter_periodic;

close(hhh)

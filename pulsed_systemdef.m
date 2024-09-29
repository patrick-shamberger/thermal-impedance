% Script 'pulsed_systemdef.m' v1.0.0, tested 28 September 2024
% Written by P.J. Shamberger (c) 2023

%% PULSED_SYSTEM_DEF ------------------------------------------------------
% Input file that contains system description required to call:
% 'pulsed_1D_time' nad 'single_pulse_1D'.  This contains 5 differents 
% structs:
%
% boundary.  Contains parameters related to the thermal boundary conditions
% materials.  Contains parameters related to the material properties
% timing.  Contains parameters related to the periodic pulses
% geometry.  Contains parameters describing the geometry of the system
% booleans.  Contains control parameters which select options for output.


%% BOUNDARY CONDITIONS ....................................................
% ++++ Heating (hotside) BC's +++++++++++++++++++++++++++++++++++++++++++++
    boundary.BC_flag = 0; % [0 = heat flux BC; 1 = temperature BC].
    boundary.pulse_flag = 1; % [1 = pulse train; 0 = sine wave].


        % .... PULSE TRAIN ........................................
        % Guidance: 10^-5 to 10^3 s
        boundary.t1 = 0.001; % s.  Rectangular 'ON' pulse duration
        
        % Guidance: 0.5 to 0.01
        % note: t2 = boundary.t1/boundary.df - boundary.t1;    % s.  Rectangular 'OFF' pulse duration        
        boundary.df = 0.1;   % duty factor.  df = t1 / (t1 + t2);
        
        % Debug tool: show q" pulse function
        % boundary.df = 1;  % constant q"


        % .... SINE WAVE ..........................................
        % ---- frequency of BC's ----------------------------------
        %boundary.heat_frequency_matrix = logspace(-2, 5, 71);   % high res sampling
        boundary.heat_frequency = [1E-0]; %0.01


        % .... HEATING BC's .......................................
        % ---- periodic q" BC's -----------------------------------
        boundary.q = [0.433E4];         % [W/m2]  = 0.433 W/cm2
        boundary.q_offset_pulse = 0;    %"background" power (power dissipated during "off" pulse).
            
        % ---- periodic T BC's -----------------------------------
        boundary.T_heat_base = [1E1];       % [C], average T for periodic BC
        boundary.T_heat_amplitude = boundary.T_heat_base;     % [C], amplitude of periodic T for periodic BC


% ++++ Cooling (cool side) BC's +++++++++++++++++++++++++++++++++++
    % For Reference:
        % Free Convection - air, gases and dry vapors : 0.5 - 1000 (W/(m2K))
        % Free Convection - water and liquids: 50 - 3000 (W/(m2K))
        % Forced Convection - air, gases and dry vapors:  10 - 1000 (W/(m2K))
        % Forced Convection - water and liquids:  50 - 10000 (W/(m2K))
     
        % OR...
        % 1e3 (W/(m2K)) equivalent to 10 K/W heat sink (free air...)
        % 1e4 (W/(m2K)) equivalent to 1 K/W heat sink (forced air...)
        % 1e5 (W/(m2K)) equivalent to 0.1 K/W heat sink (forced water...)
    
    % Cooling on heat sink side -----------------------------------
    boundary.Cooling_hs_flag = 1;       % [0 = adiabatic BC; 1 = cooling BC].
    boundary.Cooling_hs_BC_flag = 1;    % [1 = Constant T BC; 2 = convection BC; 3 = constant q BC].
    
    %Note: if h_cooling is too large for a given timestep, goes crazy (overcools!)
    boundary.h_cooling_hs = [10000]; % h [W·m-2·K-1], convective cooling coeff
    % 10, 15, 25 C
    boundary.T_cooling_hs = [15]; % T [C], cooling BC temperature
    boundary.q_cooling_hs = [1]; % q" [W/m2], cooling BC heat flux

%% MATERIALS ..............................................................
    % indices, referring to a particular row in 'Materials Properties.mat'
    materials.PCM = 1;      % [1 = octadecane; rows 1-9 a].  PCM type.
    materials.Tm = 18;      % Tm of PCM (allows a user to over-ride the defined Tm)
    materials.HighK = 20;   % [20 = OHFC Cu/C10100]. high-K material.
    materials.Vfp = 0.967;    % Porosity of Cu foam
    
    materials.top_bond = 29;       % top interface: thermal epoxy/Al foil/thermal epoxy
    materials.bottom_bond = 30;    % bottom interface: thermal paste/Al foil/thermal epoxy

%% TIMING .................................................................
    %***NOTE: define pt & dt_step dynamically within each loop 
        % period = 1/frequency
        % dt = period./pts_per_cycle
        % pt = period.*n_cycles 
    % High res: (good for higher numerical resolution.  increase as needed)
    %timing.pts_per_cycle = 5000;      % # of sampling points/cycle
    timing.pts_per_cycle = 10000;     % # of sampling points/cycle
    
    % Test Low res: (good for testing)
    %timing.pts_per_cycle = 1000;      % # of sampling points/cycle
    
    % Variables related to how many cycles to run before a system achieves
    % quasi-steady state.
    timing.n_cycles = 1000;             % total (max) number of cycles
    %timing.n_cycles = 20;             % total (max) number of cycles
    timing.n_cycles_min = 10;          % minimum # of cycles
    timing.n_oscill_test = 5;          % Melt interface oscillations over last n cycle (note: n_cycles_min > n_oscill_test)
    %timing.convergence_delta = 1E-5;   % ratio of delta(T_max)/T_max for "convergence"
    timing.convergence_delta = 1E-4;   % ratio of delta(T_max)/T_max for "convergence"


%% CONTROL parameters .....................................................
    control.recordvid_flag = 0;     % Record a video with boolean?
    control.recorddata_flag = 0;    % Record data to excel sheets?
    control.endsimulation_flag = 0; % ends the simulation after all PCM has melted, experimental
    control.plotdata_flag = 0;      % Generate plots?
    control.convergence_flag = 1;   % Stop cycles after meeting "convergence"

    control.savefile = 'output.mat';
    
%% GEOMETRY ...............................................................

%           ___top_(heated)___ 
%          |                  |
%          | 1) thermal epoxy |  H_bond_top
%          |__________________|<<<  Heated at die top surface
%          |                  |
%          |                  |  
%  Adiabat | 2) PCM slab      |  H_slab      Adiabat
%          |                  |  
%          |__________________|
%          |                  |
%          | 3) thermal paste |  H_bond_bot
%          |__________________|
%             Bottom (cooled)

    geometry.m = 1001;    % # nodes in y direction.  Note: convergence testing.  within 5% of analytical solution!
    %geometry.m = 200;    % # nodes in y direction.  Note: Low RESOLUTION. ok for testing.  Not OK for publication figs.
    geometry.n = 1;       % # nodes in x direction (1D model!)
        
    geometry.H_bond_top = 0.00025;   %[m] top interface: thermal epoxy/Al foil/thermal epoxy
    geometry.H_slab = 0.0107;       %[m] slab thickness
    geometry.H_bond_bot = 0.00025;   %[m] bottom interface: thermal paste/Al foil/thermal epoxy
    
    % Total height
    geometry.H = geometry.H_bond_top + geometry.H_slab + geometry.H_bond_bot;

    geometry.L = 1;      %[m] width.  set to 1 m (1D model!)
    geometry.w2= 0;      %[m] thickness of metal base layer.

    % Area of die
    geometry.area = 0.0006002;     %[m2] = 1 cm2.  ONLY used for power normalization.
    % Total power = q" x area = geometry.area x boundary.q
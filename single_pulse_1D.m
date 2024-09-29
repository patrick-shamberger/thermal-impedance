% Script 'single_pulse_1D.m' v1.0.0, tested 28 September 2024
% Written by P.J. Shamberger (c) 2023

%% SINGLE_PULSE_1D ------------------------------------------------------
% Finite Difference Scheme: 1D model, periodic boundary conditions
%
% single_pulse_1D.m is a finite difference scheme that performs finite 
% difference analysis on a 1D PCM slab with either a periodic heat flux or 
% temperature boundary condition. The 2 vertical sides are adiabatic, and 
% the top surface has a constant BC (constant q, T, or convect BC).


%% Initialize Variables

% MODEL VARIABLES (start) -------------------------------------------------
    % BOUNDARY CONDITIONS .................................................
    % ++++ Heating (hotside) BC's +++++++++++++++++++++++++++++++++++++
        BC_flag = boundary.BC_flag; % [0 = heat flux BC; 1 = temperature BC].   
        pulse_flag = boundary.pulse_flag; % [1 = pulse train; 0 = sine wave].
        
            % .... SINE WAVE ..........................................
            % ---- frequency of BC's ----------------------------------
            % heat_frequency_matrix = boundary.heat_frequency_matrix;
            heat_frequency = boundary.heat_frequency;
            % ---- periodic q" BC's -----------------------------------
            q_amplitude = boundary.q;
            offset_sin = boundary.q;
            % ---- periodic T BC's -----------------------------------
            T_heat_base = boundary.T_heat_base;           % [C], average T for periodic BC
            T_heat_amplitude = boundary.T_heat_amplitude; % [C], amplitude of periodic T for periodic BC

            % .... PULSE TRAIN ..........................................
            t1 = boundary.t1;     % s.  Rectangular 'ON' pulse duration.  Set up as a matrix to sweep through multiple values of 'ON' pulse
            df = boundary.df;   % duty factor.  df = t1 / (t1 + t2);
            
    % ++++ Cooling (heat sink side) BC's ++++++++++++++++++++++++++++++++++
        Cooling_hs_flag = boundary.Cooling_hs_flag;       % [0 = adiabatic BC; 1 = cooling BC].
        Cooling_hs_BC_flag = boundary.Cooling_hs_BC_flag; % [1 = Constant T BC; 2 = convection BC; 3 = constant q BC].

        %Note: if h_cooling is too large for a given timestep, goes crazy (overcools!)
        h_cool = boundary.h_cooling_hs;  % h [W·m-2·C-1], convective cooling coeff
        T_cool = boundary.T_cooling_hs;  % T [C], cooling BC temperature
        q_cool = boundary.q_cooling_hs; % q" [W·m-2], cooling BC heat flux

    % MATERIALS ...........................................................
        PCM = materials.PCM; %[4 = LNH].  PCM type.  See list below.
        HighK = materials.HighK; %[10 = Cu]. high-K material.  See list below.
        Vfp = materials.Vfp; %Volume fraction of PCM

        index_top_bond = materials.top_bond;        % top interface: thermal epoxy/Al foil/thermal epoxy
        index_bottom_bond = materials.bottom_bond;  % bottom interface: thermal paste/Al foil/thermal epoxy
        
    % TIMING ..............................................................
        %***NOTE: define pt & dt_step dynamically within each loop 
            % period = 1/frequency
            % dt = period./pts_per_cycle
            % pt = period.*n_cycles 
        pts_per_cycle = timing.pts_per_cycle;    % # of sampling points/cycle
        n_cycles = 1;       % total number of cycles
        n_cycles_min = 1;   % minimum # of cycles
        n_oscill_test = 1;  % Melt interface oscillations over last n cycle (note: n_cycles_min > n_oscill_test)
        convergence_delta = timing.convergence_delta;   % ratio of delta(T_max)/T_max for "convergence"

    % GEOMETRY ............................................................
        m = geometry.m;     % # nodes in y direction.        
        n = geometry.n;     % # nodes in x direction (1D model!)
        H = geometry.H;     %[m] height. (direction perpendicular to wall)
        L = geometry.L;     %[m] width.  set to 1 m (1D model!)
        w2= geometry.w2;    %[m] thickness of metal base layer.
        
        H_bond_top = geometry.H_bond_top;   %[m] top interface: thermal epoxy/Al foil/thermal epoxy
        H_slab = geometry.H_slab;           %[m] slab thickness
        H_bond_bot = geometry.H_bond_bot;   %[m] bottom interface: thermal paste/Al foil/thermal epoxy
    
    % CONTROL parameters ..................................................
        recordvid_flag = control.recordvid_flag;            % Record a video with boolean?
        recorddata_flag = control.recorddata_flag;          % Record data to excel sheets?
        endsimulation_flag = control.endsimulation_flag;    % ends the simulation after all PCM has melted, experimental
        plotdata_flag = control.plotdata_flag;              % Generate plots?
        convergence_flag = control.convergence_flag;        % Stop cycles after meeting "convergence"

        savefile = control.savefile;
        
% MODEL VARIABLES (end) ---------------------------------------------------            

% DEFINE OUTPUT MATRICES---------------------------------------------------
mmm = pts_per_cycle*n_cycles - 1;

time_out = zeros(mmm, 1);
Ymelt_out = zeros(mmm, 1);
T_cold_out = zeros(mmm, 1);
T_slab_out = zeros(mmm, 1);
T_hot_out = zeros(mmm, 1);

if ~BC_flag      % selects amplitude of periodic q" BC
    qsin_out = zeros(mmm, 1);
    qpulse_out  = zeros(mmm, 1); 
end

if BC_flag
    qbottom_out = zeros(mmm, 1);
end  

% DEFINE OUTPUT MATRICES-(End)---------------------------------------------
tic


%% Setup each particular iteration ..........................%%
%clear variables
clear qlatent
clear qtot2
clear qsense
clear qtot1
clear qfelts
clear qlog
clear qanaly
clear qnorm
clear q
clear DelT

% setup BC parameters for this particular iteration ............................................................

t2 = t1/df - t1;    % s.  Rectangular 'OFF' pulse duration
t_tot = t1 + t2;    % s.  Total duration of period


% dynamic calculation of dt, pt, etc.
% % % % pts_per_cycle = # of sampling points/cycle
% % % % n_cycles = total number of cycles

if pulse_flag % this happens if q" = pulse train
    % Define time vector.......................................
    period = t_tot;                     % [s], period of periodic BC
    dt_step = period./pts_per_cycle;    % [s], simulation timestep
    pt = period.*n_cycles;              % [s], total time of simulation

    time1 = dt_step:dt_step:pt;         % [s], time vector
    dt = dt_step.*ones(1, n_cycles.*pts_per_cycle);   % [s], dt vector   

    % Define q" vector.........................................
    offset_pulse = boundary.q_offset_pulse;

    dummy_vec = find(mod(time1,period)<t1); % find all points during "on" pulse
    qpulse = zeros(1,length(time1));
    qpulse(dummy_vec) = 1;

    qpulse = q_amplitude.*qpulse + offset_pulse;

    % Debug tool: show q" pulse function
    % figure(1)
    %     plot(time1, qpulse,'r.')

else % this happens if q" = sinewave
    % Define time vector.......................................
    period = 1./heat_frequency;         % [s], period of periodic BC
    dt_step = period./pts_per_cycle;    % [s], simulation timestep
    pt = period.*n_cycles;              % [s], total time of simulation

    time1 = dt_step:dt_step:pt;         % [s], time vector
    dt = dt_step.*ones(1, n_cycles.*pts_per_cycle);   % [s], dt vector

    % Define q" vector.........................................
    if ~BC_flag % this happens if variable q" BC
        qsin = q_amplitude*cos(2.*pi().*heat_frequency.*time1) + offset_sin;
    else  % this happens if variable T BC
        T_heating = T_heat_base + T_heat_amplitude*cos(2.*pi().*heat_frequency*time1);
    end
end

% -------------------------------------------------------------------------
%Derived Material Property Variables
% 1 = Top Bond
% 2 = Composite PCM (Cu foam + PCM)
% 3 = High-K Scaffold (Empty Cu foam)
% 4 = Bottom Bond

% Note: PCM is treated as a laminar composite using the following relations:
%   PCM = materials.PCM; %[1 = octadecane].  PCM type.
%   HighK = materials.HighK; %[20 = OHFC Cu/C10100]. high-K material.
%   Vfp = materials.Vfp; %Volume fraction of PCM

% PROPERTIES of composite PCM
rho_comp = Mat(PCM,2).*Vfp + Mat(HighK,2).*(1-Vfp); % Density (of solid); kg/m^3
% Note: thermal conductivity can be calculated for composite materials using a rule of mixtures or based on a specific measured thermal conductivity.  Here, we do the latter.
% kl_comp = Mat(PCM,1).*Vfp + Mat(HighK,8).*(1-Vfp); % thermal conductivity (of liquid); W/m/K
% ks_comp = Mat(PCM,8).*Vfp + Mat(HighK,8).*(1-Vfp); % thermal conductivity (of solid); W/m/K
kl_comp = 4.8; % thermal conductivity (of liquid); W/m/K
ks_comp = 4.8; % thermal conductivity (of solid); W/m/K
cpl_comp = Mat(PCM,3).*Vfp.*Mat(PCM,2)./rho_comp + Mat(HighK,9).*(1-Vfp).*Mat(HighK,2)./rho_comp; % volumetric heat capacity (of liquid); kj/m3/K 
cps_comp = Mat(PCM,9).*Vfp.*Mat(PCM,2)./rho_comp + Mat(HighK,9).*(1-Vfp).*Mat(HighK,2)./rho_comp; % volumetric heat capacity (of solid); kj/m3/K 
Lw_comp = Mat(PCM,4).*Vfp.*Mat(PCM,2)./rho_comp;    % Gravimetric Latent Heat of Fusion; j/kg 
Lv_comp = Mat(PCM,5).*Vfp;  % Volumetric Latent Heat of Fusion; j/m^3

% PROPERTIES of *just conductive* copper scaffold (units as above)
rho_scaff = Mat(HighK,2).*(1-Vfp);
%kl_scaff = Mat(HighK,8).*(1-Vfp);
%ks_scaff = Mat(HighK,8).*(1-Vfp);
kl_scaff = 4.8;
ks_scaff = 4.8;
cpl_scaff = Mat(HighK,9).*(1-Vfp).*Mat(HighK,2)./rho_comp;
cps_scaff = Mat(HighK,9).*(1-Vfp).*Mat(HighK,2)./rho_comp;
Lw_scaff = 1;
Lv_scaff = 1;

% Modified to account for "porosity" of bond layer
k_top_bond = 5; % W/m/K
k_bottom_bond = 0.1; % W/m/K

% Vectors of material properties, corresponding to:
% 1 = Top Bond
% 2 = Composite PCM (Cu foam + PCM)
% 3 = High-K Scaffold (Empty Cu foam)
% 4 = Bottom Bond
kl = [k_top_bond, kl_comp, kl_scaff, k_bottom_bond];    % liquid thermal conductivity in W/mK [PCM, metal]
ks = [k_top_bond, ks_comp, ks_scaff, k_bottom_bond];    % solid thermal conductivity in W/mK [PCM, metal]
rho = [Mat(index_top_bond,2), rho_comp, rho_scaff, Mat(index_bottom_bond,2)];  % density in kg/m^3 [PCM, metal]
cpl = [Mat(index_top_bond,3), cpl_comp, cpl_scaff, Mat(index_bottom_bond,3)];  % specific heat of liquid in J/kgK [PCM, metal]
cps = [Mat(index_top_bond,9), cps_comp, cps_scaff, Mat(index_bottom_bond,9)];  % specific heat of solid in J/kgK [PCM, metal]
Lw = [Mat(index_top_bond,4), Lw_comp, Lw_scaff, Mat(index_bottom_bond,4)];    % Specific Latent Heat of Fusion in J/kg [PCM,metal]

%Tm = Mat(PCM,7); %Melting temperature in [C] of the PCM
Tm = materials.Tm;
%alpha_l = kl./rho./cpl; %Thermal diffusivity in m^2/s
%FOM = sqrt(kl(1)*Lv_comp); %figure of merit for PCM materials for constant T BC


%% Setup initial state/geometry of nodes

%Tinf = Tm; %Initial temperature of PCM [C]
Tinf = T_cool; %Initial temperature of PCM [C]

% Be smart about initial temperature...
% if ~BC_flag  % [0 = heat flux BC].
%     if Cooling_hs_flag == 1 & Cooling_hs_BC_flag == 2;       % [2 = convection BC].
%         %Tinf = q./h_cool; %equilibrium temperature of PCM [C]
%         Tinf = T_cool;
%     elseif Cooling_hs_flag == 1 & Cooling_hs_BC_flag == 1;       % [1 = constant T].
%         %Tinf = q./h_cool; %equilibrium temperature of PCM [C]
%         Tinf = T_cool;
%     end
% end

%Set up 1D nodal network (1D: along y-axis)
Y=zeros(m,1);
dely = H/(m-1); % Vertical distance between nodes [m]

% Define materials
P = ones(m,1).*7; % Initializes all nodes to be heat sink alloy to start with, P(i,j) determines if a material is PCM or high K

dtn = length(dt); %number of iterations performed in analysis

%Declare a bunch of zeros to be populated during analysis
walltemp = zeros(dtn,1);
wTanaly = zeros(dtn,1);
Ymelt = zeros(dtn,1);
Tgraph_analy = zeros(m,1);
wTnorm = zeros(dtn,1);
qtot2 = zeros(dtn,1);
Melt_Difference = zeros(dtn,1);

%% Determine material of nodes

%           ___top_(heated)___ 
%          |                  |
%          | 1) thermal epoxy |  H_bond_top
%          |__________________|<<<  Heated at die top surface
%          |                  |
%          |                  |  
%  Adiabat | 2 or 3) PCM slab |  H_slab      Adiabat
%          |                  |  
%          |__________________|
%          |                  |
%          | 4) thermal paste |  H_bond_bot
%          |__________________|
%             Bottom (cooled)


%Derived Material Property Variables
% 1 = Top Bond
% 2 = Composite PCM (Cu foam + PCM)
% 3 = High-K Scaffold (Empty Cu foam)
% 4 = Bottom Bond

% Determine material for each element
for i = 1:1:m %step in y
    Y(i) = (i-1)*dely;

    %          origin starts at top
    %          1           | (0)
    %          |           |
    %          |           |
    %          (index i)   Y(i)
    %          |           |
    %          |           |
    %          V           V (H)
    
    if  Y(i) < H_bond_top
        P(i) = 1; % 1 = Top Bond (thermal epoxy/Al foil/thermal epoxy)
    elseif  Y(i) < (H_bond_top + H_slab)
        P(i) = 2; % 2 = composite PCM; 3 = empty Cu Foam (switch)
        % NOTE: to look at PCM-filled foams, switch this to 2!!!
    else
        P(i) = 4; % 4 = Bottom Bond (thermal epoxy/Al foil/thermal paste)
    end
end

% Find the index of the slab closest to the heater
index_slab = find(Y > H_bond_top,1);

YboundaryIndex = find(P(:,end) == 1,1,'last');%determines y index of node that borders high-k material (address of last node with PCM)

if isempty(YboundaryIndex)
    YboundaryIndex = m;
end

%% Perform Finite Element Analysis
% Create array for distance of each node from the origin (top element)

N = m; %total number of nodes in the network
A = zeros(N);
B = zeros(N,1);
hft = zeros(N,1);
sl = zeros(N,1);
T = ones(N,1)*Tinf;
T2 = zeros(m,1); %Temperature of each node at each time
sl2 = zeros(m,1); %amount of energy pumped into a node that goes to melting that node expressed as a percent #percentmelted
K(1:m) = ks(P(:));
CP(1:m) = cps(P(:));

% %  (inherited from 2D)
% %                 vol1 = [0.5*delx*0.5*dely; ones(m-2,1)*dely*0.5*delx; 0.5*delx*0.5*dely]; %creates matrix of volumes for nodes of first and/or last column
% %                 vol2 = delx*0.5*dely*ones(1,n-2); %creates matrix of volumes for top and bottom nodes excluding corners
% %                 vol3 = delx*dely*ones(m-2,n-2); %creates matrix of volumes for interior nodes
% %                 vol4 = [vol2; vol3; vol2]; %vertically concates vol2 and 3 to create volume matrix of nodes excluding side walls
% %                 voltot = [vol1 vol4 vol1]; %horizontally concates vol1 and vol4 to create the complete volume matrix of the nodes

voltot = [0.5*dely; dely*ones(m-2,1); 0.5*dely;];

cycle_index = 0;    % keep track of how many "cycles" have gone through



for ii = 1:1:(length(time1)-1) %step in time

    %% TRANSPORT heat via conduction (Finite Difference Part)
    co = 1; % count index.  Same as ii?
    comov = 1;

    for i = 1:1:m %step in y
        if i > 1 && i < m % Internal Node 
            %All Fourier values use material properties of the node it is in, with
            %the exception of conductivity, which is taken as the average of
            %properties between nodes.

            % MODIFIED
            K3_eff = 2.*K(i).*K(i-1)./(K(i) + K(i-1));
            K4_eff = 2.*K(i).*K(i+1)./(K(i) + K(i+1));

            Fo3 = dt(ii)*K3_eff/(rho(P(i))*CP(P(i))*dely^2);
            Fo4 = dt(ii)*K4_eff/(rho(P(i))*CP(P(i))*dely^2);

            A(co,co) = -(Fo3 + Fo4 + 1); %equivalent of 1-4Fo in normal finite difference -(Fo1+Fo2+Fo3+Fo4 + 1); %equivalent of 1-4Fo in normal finite difference
            A(co,co-n) = Fo3; %co-n represents a negative step in y, quirk of A matrix
            A(co,co+n) = Fo4; %co+n represents the node above, or a positive step in y

            B(co) = -T(co);
        end

        % -----------------------------------------------------
        % TOP WALL (HEATED surface): Adiabatic OR set to heated BC (T, q" BC's, as specified)
        if i == 1   
            % MODIFIED
            K4_eff = 2.*K(i).*K(i+1)./(K(i) + K(i+1));

            Fo4 = dt(ii)*K4_eff/(rho(P(i))*CP(P(i))*dely^2);

            A(co,co) = -(2*Fo4 + 1);
            A(co,co+n) = 2*Fo4;

            B(co) = -T(co);

            if ~BC_flag % this happens if variable q" BC
                % (defined above for whole time vector)
                % qsin(ii) = Amplitude*cos(2.*pi().*heat_frequency*time1(ii)) + Offset;

                if pulse_flag % this happens if q" = pulse train
                    B(co) = -T(co) - 2*(qpulse(ii))*dt(ii)/(rho(P(i))*CP(P(i))*dely);
                else % this happens if q" = sine wave
                    % Note: ORIG
                    % B(co) = -T(co) -2*(qsin(ii))*dt(ii)/(rho(P(i))*CP(P(i))*dely);
                    B(co) = -T(co) - 2*(qsin(ii))*dt(ii)/(rho(P(i))*CP(P(i))*dely);
                end

            else  % this happens if variable T BC Note: DOUBLE CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                A(co,co) = 1;

                % (defined above for whole time vector)
                % T_heating(ii) = T_heat_base + T_heat_amplitude*cos(2.*pi().*heat_frequency*time1(ii));

                B(co) = T_heating(ii);

            end
        end

        % -----------------------------------------------------
        % BOTTOM WALL (HS surface): Adiabatic OR set to constant cooling (convec, T, q" BC's, as specified)
        if i == m
            % MODIFIED
            K3_eff = 2.*K(i).*K(i-1)./(K(i) + K(i-1));

            Fo3 = dt(ii)*K3_eff/(rho(P(i))*CP(P(i))*dely^2);

            A(co,co) = -(2*Fo3 + 1);
            A(co,co-n) = 2*Fo3; %co-n represents a negative step in y, quirk of A matrix

            B(co) = -T(co);

            if Cooling_hs_flag
                % [1 = Constant T BC; 2 = convection BC; 3 = constant q BC].
                if Cooling_hs_BC_flag == 1 % Constant T BC          GOOD!!!!
                    A(co,co) = 1;
                    A(co,co-n) = 0;

                    B(co) = T_cool;

                elseif Cooling_hs_BC_flag == 2 % Convection BC      GOOD!!!!

                    % Use A values as defined above.
                    B(co) = -T(co) -2*h_cool*dt(ii)/(rho(P(i))*CP(P(i))*dely)*(T_cool - T(co));

                elseif Cooling_hs_BC_flag == 3 % Constant q" BC.    GOOD!!!!
                    % Use A values as defined above.

                    % Note: for now, setup for q_cool to be average of q_heat
                    B(co) = -T(co) +2*q_cool*dt(ii)/(rho(P(i))*CP(P(i))*dely);

                else
                    disp('Cooling_hs_BC_flag out of bounds.')
                    stop                                
                end
            end
        end

        co = co + 1;
    end

    A = sparse(A);
    B = sparse(B);

    Tn1 = full(A\B);
    Tn2 = Tn1;

    %% PHASE TRANSITION (Enthalpy Method Part)
    for i = 1:1:length(Tn1)

        if P(i) == 2 %if the material is PCM

            if Tn1(i) > Tm && hft(i) < Lw(2) %is the temp higher than melting point?

                if hft(i) + (Tn1(i)-Tm)*CP(2) > Lw(2) %is energy enough to fully melt element?
                    Tn1(i) = Tm + (hft(i) + (Tn1(i)-Tm)*CP(2) - Lw(2))/CP(2); %determines the temp of element after melting
                    hft(i) = Lw(2);
                else
                    hft(i) = hft(i) + (Tn1(i)-Tm)*CP(2); %determine amount of element melted
                    Tn1(i) = Tm;
                end
            end

            if Tn1(i) < Tm && hft(i) > 0 %this section just in case something somehow solidifies
                if hft(i) + (Tn1(i)-Tm)*CP(2) < 0   % Note: there was a typo on this line originally!
                    Tn1(i) = Tm + (hft(i) + (Tn1(i)-Tm)*CP(2) - 0)/CP(2);
                    hft(i) = 0;
                else
                    hft(i) = hft(i) + (Tn1(i)-Tm)*CP(2);
                    Tn1(i) = Tm;
                    %                                     Tn1(i) = (Tn1(i-1)+Tn1(i+1))/2;
                end
            end

            sl(i) = hft(i)/Lw(2);

        end



    end

    co2 = 1;

    for i = 1:1:m %step in y

        T2(i,ii) = full(Tn1(co2)); %full() changes matlab clasifying T2 from a sparse matrix to a full matrix, allowing 3 indicies
        sl2(i,ii) = sl(co2);
        K(i) = kl(P(i))*sl2(i,ii)+ks(P(i))*(1-sl2(i,ii));

        %T2 is the temperature of each node in space and time. ii is time
        %sl2 is the fraction melted of each node in space and time. 0
        %is not melted at all, 1 is fully melted
        co2 = co2 + 1;

    end

    if ~isempty(find(sl2(:,ii) == 1,1,'first'))
        Ymelt_Index = find(sl2(:,ii) == 1,1,'first');
        Ymelt(ii) = dely *(m - Ymelt_Index);
        %                         T2(Ymelt_Index-1,ii) = (T2(Ymelt_Index,ii)+T2(Ymelt_Index-2,ii))/2; %%averages artifact of first unmelted node
    else
        Ymelt(ii) = H_bond_top;
    end

    %% Collect Results
    % following code calculates the total energy in each node,
    % in the sensable and latent forms, and subtracts that by the total
    % energy in the previous time step.
    % basically takes the time derivative of the energy in the system to
    % determine heat flux

    T_cold(ii,1) = T2(m,ii);
    T_slab_top(ii,1) = T2(index_slab,ii);
    T_hot(ii,1) = T2(1,ii);

    % Add up change in heat across all nodes.  (Latent + sensible contributions)
    if ii == 1
        qlatent(1) = sum(sl2(:,ii).*Lw(P(:))'.*rho(P(:))'.*voltot)/(dt(ii));
        qsense(1) = sum((T2(:,ii)-Tinf*ones(m,1)).*CP(P(:))'.*rho(P(:))'.*voltot)/(dt(ii));
    else
        qlatent(ii) = sum((sl2(:,ii)-sl2(:,ii-1)).*Lw(P(:))'.*rho(P(:))'.*voltot)/(dt(ii));
        qsense(ii) = sum((T2(:,ii)-T2(:,ii-1)).*CP(P(:))'.*rho(P(:))'.*voltot)/(dt(ii));
    end

    qtot1(ii)=qsense(ii) + qlatent(ii); %total heat flux, given for each node, in that time step
    qtot2(ii) = sum(qtot1(ii)); %total heat flux in given time step


    if ~BC_flag  % track "wall temperature" for q" BC case
        walltemp(ii)=mean(T2(end,ii)); %wall temperature for the timestep
        Th1 = walltemp(ii);
    else    % track heat flux for T BC case
        buffer = 1;     % avoid edge effects!
        qbottom(ii)= -(sum((T2(m-1-buffer,ii) - T2(m-buffer,ii))*kl(1)/dely))/L; %Dr. Johnny Felt's method of heat flux calculation

    end

    T = Tn1;

    if endsimulation_flag==1
        if sl2(1,ii) == 1
            %kills the current FEA loop and moves to the next one when the
            %melt front has reached the far side of the control volume.
            %Only designed for pure PCM experiments.
            break;
        end
    end

    % ###############################################################################################################################################################################
    if convergence_flag == 1
        if mod(ii,pts_per_cycle)==0
            cycle_index = cycle_index + 1;  % increment "cycles" index


            if cycle_index > n_cycles_min   % only check for convergence after n_cycles_min cycles
                max1 = max(T_slab_top((ii-pts_per_cycle):ii,1));
                max2 = max(T_slab_top((ii-n_oscill_test*pts_per_cycle):(ii-(n_oscill_test-1)*pts_per_cycle),1));

                test_value = abs((max1-max2)/max1);
                if test_value < convergence_delta

                    time1(ii+1:end) = [];
                    walltemp(ii+1:end) = [];
                    wTanaly(ii+1:end) = [];
                    Ymelt(ii+1:end) = [];
                    Tgraph_analy(ii+1:end) = [];
                    wTnorm(ii+1:end) = [];
                    T_cold(ii+1:end) = [];
                    T_slab_top(ii+1:end) = [];
                    T_hot(ii+1:end) = [];
                    qtot2(ii+1:end) = [];

                    break            
                end
            end
        end
    end
    % ###############################################################################################################################################################################
end



time1(ii+1:end) = [];
walltemp(ii+1:end) = [];
wTanaly(ii+1:end) = [];
Ymelt(ii+1:end) = [];
Tgraph_analy(ii+1:end) = [];
wTnorm(ii+1:end) = [];
T_cold(ii+1:end) = [];
T_slab_top(ii+1:end) = [];
T_hot(ii+1:end) = [];
qtot2(ii+1:end) = [];
%kills zeros left at the end of declared matricies when endsimulation
%activates. Shortens each matrix to size ii, reduces file bloat.

% *************************************************************
%% Principal outputs
% Time series to aggregate for each loop.
endpoint = length(time1);

time_out(1:endpoint) = time1;
Ymelt_out(1:endpoint) = Ymelt;

if ~BC_flag      % selects amplitude of periodic q" BC
    if pulse_flag % this happens if q" = pulse train
        qpulse_out(1:endpoint) = qpulse(1:endpoint);
    else % this happens if q" = sine wave
        qsin_out(1:endpoint) = qsin(1:endpoint);
    end
end

T_cold_out(1:endpoint) = T_cold;
T_slab_out(1:endpoint) = T_slab_top;
T_hot_out(1:endpoint) = T_hot;

if BC_flag
    qbottom_out(1:endpoint) = qbottom;
end

% Melt interface oscillations over last n cycle


for i=1:n_oscill_test
    max_Y(i) = max(Ymelt);
    min_Y(i) = min(Ymelt);
    
    del_Y(i) = max_Y(i) - min_Y(i);
end

avg_del_Y_out = mean(del_Y);
% diff_max_Y_out = max_Y(1)-max_Y(2);
max_Y_out = mean(max_Y);
min_Y_out = mean(min_Y);

%Lv = Volumetric Latent Heat of Fusion of PCM [j/m^3]
%Lv * del_Y = energy/area going to melting!
E_melt_out = Lv_comp.*mean(del_Y); % [J/m2]

%q"_avg * period = energy/area heating the system
if ~BC_flag      % selects amplitude of periodic q" BC
    E_input_out = offset_sin.*period; % [J/m2]
    E_ratio_out = (Lv_comp.*mean(del_Y))./(offset_sin.*period);
end

% CALCULATE PHASE LAG (Disabled for pulse train)
if ~BC_flag      % selects amplitude of periodic q" BC
    for i=1:n_oscill_test
        max_T_hot(i) = max(T_hot);
        %index_maxT(i) = min(find(T_j((endpoint - pts_per_cycle*i):(endpoint - pts_per_cycle*(i-1))) == max_Tj(i)));
        %phase_lag(i) = 2.*pi().*(index_maxT(i)-1)./pts_per_cycle;   % [radians]
        %phase_lag(find(phase_lag > 6 & phase_lag < 7)) = phase_lag(find(phase_lag > 6 & phase_lag < 7)) - 2*pi();
    end
    max_T_hot_out = mean(max_T_hot);
    %phase_lag_out = mean(phase_lag);
else
    for i=1:n_oscill_test
        max_qbase(i) = max(qbottom((endpoint - pts_per_cycle*i):(endpoint - pts_per_cycle*(i-1))));
        %index_maxq(i) = min(find(qbottom((endpoint - pts_per_cycle*i):(endpoint - pts_per_cycle*(i-1))) == max_qbase(i)));
        %phase_lag(i) = 2.*pi().*(index_maxq(i)-1)./pts_per_cycle;   % [radians]
    end
    max_q_out = mean(max_qbase);
    %phase_lag_out = mean(phase_lag);               
end


% *************************************************************       

if recorddata_flag

    if ~BC_flag
        save(fullfile(cd,'Data','Variable Heat','Raw Data Files',vidname),'Ymelt','qsin','Y','T2','Tinf','Th1','ii','Tm','H','qtot2','T2','time1')
    else
        save(fullfile(cd,'Data','1D','Raw Data Files',vidname),'Ymelt','Y','Tinf','Th1','Tm','H','T_heating','q','ii','qtot2','T2','time1')
    end

    %saves data in a matlab file with filename=vidname and location the
    %same as the file location of this matlab program. To access data, load
    %the file in matlab, similar to how the material properties are loaded
    %at the start of this program. See save documentation for more details.
end

if recordvid_flag
    close(v)
end

save(savefile)
    

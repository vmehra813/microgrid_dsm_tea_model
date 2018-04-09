%% TESTING kibam model
%04/09/16
%takes one consumer, one 75 W panel, 75 Ah battery, and conducts SOC w/
%KiBAM over 1 year
clc;
clearvars;
addpath('/Users/vmehra813/Dropbox (MIT)/[uLink]/Simulation/MatLab/uLink/Modeling Analysis/TPS');
addpath('/Users/vmehra813/Dropbox (MIT)/[uLink]/Simulation/MatLab/uLink/Modeling Analysis/Generation + Reliability Model');

addpath('/Users/ramatya/Dropbox (MIT)/[uLink]/Simulation/MatLab/uLink/Modeling Analysis/TPS');
addpath('/Users/ramatya/Dropbox (MIT)/[uLink]/Simulation/MatLab/uLink/Modeling Analysis/Generation + Reliability Model');


num.consumers = 1;
num.A = 1;
num.B = 0;
num.graph = num2str(num.consumers);
num.years = 5;

config.pv_batt = 'PVBatteryConfigurations.xlsx';
config.sheet_name = 'Data';
config.PanelSize = xlsread(config.pv_batt, config.sheet_name, 'PanelSize');
%config.PVArea = xlsread(config.pv_batt, config.sheet_name,'PVArea');
config.BattCap = xlsread(config.pv_batt, config.sheet_name, 'BattCap');
%config.PVEff = xlsread(config.pv_batt, config.sheet_name, 'PVEff');
%config.system = 8; % 7 is 75 watts // 10 is 250 watts, 180 ah 8 is 125
%Watts, 180 Ah battery [old version]

config.PVEff = 0.15; % 15 %
config.option = 1;
config.PV_Size = 75;

%include logic for battery $/kWh
config.MinDist = 10;    % meters
config.MaxDist = 200;   % meters

config.PV_area = pv_area(config);


solar.solaroutputs=xlsread('pvwatts_hourly_patamda_jharkhand.xlsx','DATA','A20:K8779');
%solar.irradiance_beam = solar.solaroutputs(:,4);       % plane of array irradiance (W/m^2) (beam = 4)
solar.irradiance_plane = solar.solaroutputs(:,8);       % plane of array irradiance (W/m^2) (beam = 4)
solar.tAmbient=solar.solaroutputs(:,6);                 % ambient temperature in C
solar.hourofday=solar.solaroutputs(:,3);
solar.month=solar.solaroutputs(:,1);
%solar.temp_coeff = 0.06;   % degradation per degree C above STC // Tata Power Solar
solar.temp_coeff = 0.0041;  % literature on mono-si
solar.temp_STC = 25;        % 25 degrees C reference temperature

%visual plot of distances / nodes

% Call Function
[a] = singlehome_generation_VM(config, solar,num);

demand.inputs = 'demandInputs_uLink.xlsx';
demand.type.residential = 'Residential';
demand.type.anchor = 'Anchor';

demand.distr_loss = 0.9;

%incorporating losses
demand.min_distance = config.MinDist;   % meters
demand.max_distance = config.MaxDist;   % meters
demand.wire_resistivity = 0.00659;      %Ohm / foot
% NEED WIRE SIZE
demand.distances_b = demand.wire_resistivity*((demand.max_distance - demand.min_distance)*rand(num.consumers,1)+demand.min_distance); %vector of distances for b1-bn in a radial format from A box
demand.distr_voltage = 24;
demand.CNSE = 0.01;

%Call Function
[loadprof, load] = CreateDemandProfiles_VM(num,demand, solar);

%% Battery Parameters
%NOT DIVIDING BY 1000 FOR KWH CURRENTLY

config.Batt_Size = 75;                              % Battery capacity (Ah)
battery.charge_eff = 0.95;                          % Battery charging efficiency
battery.discharge_eff = 0.95;                       % Battery discharging efficieny
battery.voltage = 12;                               % Battery operating voltage [V]
battery.size = config.Batt_Size*battery.voltage;    % (7) is 75 Ah; (10) is 180 Ah % capacity (Wh)
%battery.initial_capacity = battery.size / 1;       % initial capacity (100%)
battery.MDOD = .6;                                  % maximum depth of discharge (example = 60% i.e. SOC can be upto 40% before cut off)
batt = config.Batt_Size;                            % Wh
num_batt = 1;

deltaT = 1;     % Time step (1 hour)
alpha_c = 1;    % max charge rate in A/Ah i.e same as saying 1 Ah capacity battery charging at 1 hour rate with 1 A input current
Imax = 12;      % 12 amps max current

%define these for HBL battery! have inputs in main script under battery
%data structure
% c = battery.c;
% k = battery.k;

c = 0.305;% capacity ratio
k = 2.12; % rate constant

%round trip efficiency sqrt is what HOMER uses for charge/discharge eff


Qmax_batt = batt/1;                 % currently have this in Ah (maximum battery capacity)
Qmax_battbank = Qmax_batt*num_batt; % Qmax_battbank is max capacity of the battery bank in Ah

e_raised = exp(-k*deltaT);

Q_1 = zeros(size(loadprof.total)); % some have an added a row bc of t+1 iteration later on Available charge (Ah)
Q_2 = zeros(size(loadprof.total)); % bound charge (Ah)
%battEcum = zeros(size(loadprof.total));
battEoutMax = zeros(size(loadprof.total));  % energy out max from the battery during discharge (Wh)
battEinMax = zeros(size(loadprof.total));   % energy in max into the battery during charge (Wh)
Ebatt_max_in_lim1 = zeros(size(loadprof.total));
Ebatt_max_in_lim2 = zeros(size(loadprof.total));
Ebatt_max_in_lim3 = zeros(size(loadprof.total));
Ebatt_max_in_limSOCmax = zeros(size(loadprof.total));
EoutMax_SOCminlimit = zeros(size(loadprof.total));
EoutMax_kineticmodellimit = zeros(size(loadprof.total));
battSOC = zeros(size(loadprof.total));      % battery state of charge (as %)
Q_total = zeros(size(loadprof.total));      % total charge (available + bound) (Ah)


battery.initial_capacity = 0.8; % initial capacity as % (initial available and bound charge is 80% of max capacity * c) <-- this is the initial SOC

%% Defining Battery Charge / Discharge
%these are matrices
%power.in_out_batt_noncritical = a.available_solar - loadprof.total;
%power.in_out_batt_critical = a.available_solar - loadprof.final.critical;

power.in_out_batt_noncritical = loadprof.total- a.available_solar ; % (unit: W) sign convention: positive: discharge, negetive: charge 
power.in_out_batt_critical = loadprof.final.critical - a.available_solar;


current.in_out_batt_noncritical = power.in_out_batt_noncritical / battery.voltage; % (unit: A) sign convention: positive - discharge, negetive - charge 
current.in_out_batt_critical = power.in_out_batt_critical / battery.voltage;


%% Initial Condition 
Q_1(1) = c*battery.initial_capacity*Qmax_battbank;      % these are the starting values of Q_1 and Q_2 for the first iteration - this is a property of the battery bank (Available charge Ah)
Q_2(1) = (1-c)*battery.initial_capacity*Qmax_battbank;  % bound charge (Ah)

%battEcum(1) = 0;
Q_total(1) = Q_1(1) + Q_2(1);               % unit Ah
battSOC(1) = Q_total(1)/(Qmax_battbank);    % should be equal to initial capacity (0.8 in this case)

%removed -1*
EoutMax_kineticmodellimit(1) = battery.discharge_eff*deltaT*((-k*c*Qmax_battbank + k*Q_1(1)*e_raised + Q_total(1)*k*c*(1-e_raised))/...
    (1-e_raised+c*(k*deltaT-1+e_raised)));

%EoutMax_SOCminlimit = (battSOC(1)-battery.MDOD)*Qmax_battbank;
EoutMax_SOCminlimit(1) = (1-battery.MDOD)*Qmax_battbank;

battEoutMax(1)= max(EoutMax_kineticmodellimit(1),EoutMax_SOCminlimit(1));


Ebatt_max_in_lim1(1) = deltaT*k*Q_1(1)*e_raised+Q_total(1)*k*c*(1-e_raised)/(1-e_raised+c*(k*deltaT-1+e_raised));
Ebatt_max_in_lim2(1) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-Q_total(1))/deltaT;
Ebatt_max_in_lim3(1) = deltaT*num_batt*Imax*battery.voltage/1;
%Ebatt_max_in_limSOCmax = (BESS.SOCmax-battSOC(1))*Qmax_battbank;
%Ebatt_max_in_limSOCmax = (1-battSOC(1))*Qmax_battbank;
Ebatt_max_in_limSOCmax(1) = Qmax_battbank;

% energy into the battery during charge (Wh)
battEinMax(1) = min([Ebatt_max_in_lim1(1), Ebatt_max_in_lim2(1), Ebatt_max_in_lim3(1), Ebatt_max_in_limSOCmax(1)])/battery.charge_eff;

%NonCritical
%CHANGED THIS TO BE FOR NON-CRITICAL LOADS
%initialize vectors
network.total.reliability = zeros(num.consumers,1);
network.total.amountserved.cumulative = zeros(length(loadprof.total),num.consumers);
network.total.nonserved = zeros(length(loadprof.total),num.consumers);
%network.total.amountserved.cumulative(1) = 20;
network.total.solar_load_coincidence = zeros(length(loadprof.total),num.consumers);


SOC.total = zeros(length(loadprof.total), num.consumers);
SOC.total_spillage = zeros(length(loadprof.total),num.consumers);
SOC.total_cycles = zeros(num.consumers,1);

%% Noncritical

for i=1
    
    if current.in_out_batt_noncritical(i) < 0     %charging
        
        % need other 1st condition to include solar/load for hour 1
        % Q_1(i) = Q_1(i)*e_raised+(Q_total(i)*k*c+current.in_out_batt_noncritical(i))*(1-e_raised)/k+power.in_out_batt_noncritical(i)*c*(k*deltaT-1+e_raised)/k;
        % Q_2(i) = Q_2(i)*e_raised+(Q_total(i)*(1-c))*(1-e_raised)+current.in_out_batt_noncritical(i)*(1-c)*(k*deltaT-1+e_raised)/k;
        % Q_total(i) = Q_1(i) + Q_2(i);
        % battSOC(i) = Q_total(i)/Qmax_battbank;
         
        Ebatt_max_in_lim1(i) = deltaT*(k*Q_1(i)*e_raised+Q_total(i)*k*c*(1-e_raised))/(1-e_raised+c*(k*deltaT-1+e_raised));
        Ebatt_max_in_lim2(i) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-Q_total(i))/deltaT;
        Ebatt_max_in_lim3(i) = deltaT*num_batt*Imax*battery.voltage/1;
        Ebatt_max_in_limSOCmax(i) = (1-battSOC(i))*Qmax_battbank;
        % energy into the battery during charge (Wh)
        battEinMax(i) = min([Ebatt_max_in_lim1(i), Ebatt_max_in_lim2(i), Ebatt_max_in_lim3(i), Ebatt_max_in_limSOCmax(i)])/battery.charge_eff;
        
        % Total spillage (Wh)
        SOC.total_spillage(i) = abs(a.available_solar(i)-loadprof.final.non_critical(i)-battEinMax(i))*ge(a.available_solar(i)-loadprof.final.non_critical(i),battEinMax(i));
      
        
    elseif current.in_out_batt_noncritical(i) > 0 %discharging
        
        % need other 1st condition to include solar/load for hour 1
        
        % Q_1(i) = Q_1(i)*e_raised+(Q_total(i)*k*c+current.in_out_batt_noncritical(i))*(1-e_raised)/k+power.in_out_batt_noncritical(i)*c*(k*deltaT-1+e_raised)/k;
        % Q_2(i) = Q_2(i)*e_raised+(Q_total(i)*(1-c))*(1-e_raised)+current.in_out_batt_noncritical(i)*(1-c)*(k*deltaT-1+e_raised)/k;
        % Q_total(i) = Q_1(i) + Q_2(i);
        % battSOC(i) = Q_total(i)/Qmax_battbank;
        
       
        EoutMax_kineticmodellimit(i) = battery.discharge_eff*deltaT*((-k*c*Qmax_battbank + k*Q_1(i)*e_raised + Q_total(i)*k*c*(1-e_raised))/...
            (1-e_raised+c*(k*deltaT-1+e_raised)));
        EoutMax_SOCminlimit(i) = (battSOC(i)-battery.MDOD)*Qmax_battbank;
        % energy out from the battery during discharge (Wh)
        battEoutMax(i)= min([EoutMax_kineticmodellimit(i),EoutMax_SOCminlimit(i)]);
        
        % Total non-served energy (Wh)
        network.total.nonserved(i) = abs(loadprof.total(i) + battEoutMax(i)-a.available_solar(i))*ge(loadprof.total(i)-a.available_solar(i),abs(battEoutMax(i)));
        
    end
    
end

for i=2:length(loadprof.final.non_critical)
    if current.in_out_batt_noncritical(i) < 0       % charging 
        
        Q_1(i) = Q_1(i-1)*e_raised+(Q_total(i-1)*k*c-current.in_out_batt_noncritical(i))*(1-e_raised)/k-current.in_out_batt_noncritical(i)*c*(k*deltaT-1+e_raised)/k;
        Q_2(i) = Q_2(i-1)*e_raised+(Q_total(i-1)*(1-c))*(1-e_raised)-current.in_out_batt_noncritical(i)*(1-c)*(k*deltaT-1+e_raised)/k;
        Q_total(i) = Q_1(i) + Q_2(i);
        battSOC(i) = abs(Q_total(i)/Qmax_battbank);      % + battSOC(i-1);
        
 
        if battSOC(i)+battSOC(i-1) < 1
            battSOC(i) = battSOC(i)+battSOC(i-1);
        else
            %battSOC(i) = battSOC(i-1);
            battSOC(i) = 1;
        end
       
        %spillage: if Q_1 is larger than batt_e_in_max
        %set to previous condition
        %limits
        Ebatt_max_in_lim1(i) = deltaT*(k*Q_1(i)*e_raised+Q_total(i)*k*c*(1-e_raised))/(1-e_raised+c*(k*deltaT-1+e_raised));
        Ebatt_max_in_lim2(i) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-Q_total(i))/deltaT;
        Ebatt_max_in_lim3(i) = deltaT*num_batt*Imax*battery.voltage/1;
        Ebatt_max_in_limSOCmax(i) = (1-battSOC(i))*Qmax_battbank;
        battEinMax(i) = min([Ebatt_max_in_lim1(i), Ebatt_max_in_lim2(i), Ebatt_max_in_lim3(i), Ebatt_max_in_limSOCmax(i)])/battery.charge_eff;
        
        SOC.total_spillage(i) = abs(a.available_solar(i)-loadprof.final.non_critical(i)-battEinMax(i))*ge(a.available_solar(i)-loadprof.final.non_critical(i),battEinMax(i));
    
    elseif current.in_out_batt_noncritical(i) == 0
        
        %Q_1(i) = Q_1(i-1);
        %Q_2(i) = Q_2(i-1);
        %Q_total(i) = Q_total(i-1);
        battSOC(i) = battSOC(i-1);
        
    elseif current.in_out_batt_noncritical(i) > 0       % discharging 
        
        Q_1(i) = Q_1(i-1)*e_raised+(Q_total(i-1)*k*c - current.in_out_batt_noncritical(i))*(1-e_raised)/k - current.in_out_batt_noncritical(i)*c*(k*deltaT-1+e_raised)/k;
        Q_2(i) = Q_2(i-1)*e_raised+(Q_total(i-1)*(1-c))*(1-e_raised) - current.in_out_batt_noncritical(i)*(1-c)*(k*deltaT-1+e_raised)/k;
        Q_total(i) = Q_1(i) + Q_2(i);
        battSOC(i) = abs(Q_total(i)/Qmax_battbank);          % + battSOC(i-1);
        
        if battSOC(i-1) - battSOC(i) > (1- battery.MDOD)
            battSOC(i) = battSOC(i-1) - battSOC(i);
            
        else
            %battSOC(i) = battSOC(i-1);
            battSOC(i) = 1 - battery.MDOD;
            % is this right 
        end
        
        %limits
        EoutMax_kineticmodellimit(i) = battery.discharge_eff*deltaT*((-k*c*Qmax_battbank + k*Q_1(i)*e_raised + Q_total(i)*k*c*(1-e_raised))/...
            (1-e_raised+c*(k*deltaT-1+e_raised)));
        EoutMax_SOCminlimit(i) = (battSOC(i)-battery.MDOD)*Qmax_battbank;
        battEoutMax(i)= min([EoutMax_kineticmodellimit(i),EoutMax_SOCminlimit(i)]);
        
        %check conditions on limits
        %in current (A)
        %battEoutMax is negative
        network.total.nonserved(i) = abs(loadprof.total(i) + battEoutMax(i)-a.available_solar(i))*ge(loadprof.total(i)-a.available_solar(i),abs(battEoutMax(i)));
        
    end
    
end

reliability = 1 - sum(network.total.nonserved)/sum(loadprof.total);
disp(reliability);


figure(3);
plot(battSOC);
title('KiBam Battery SOC: 75Ah, 1 Customer');
ylabel('Percent SOC (%)');
xlabel('Hour');
axis([0 8800 0 1.05]);

figure(4);
plot(network.total.nonserved);
title('NonServed Energy, Hour by Hour');

figure(5);
plot(SOC.total_spillage);
title('Spilt Solar Energy, Hour by Hour');


% hold 
% plot(current.in_out_batt_noncritical/max(current.in_out_batt_noncritical))

%plot((solar.irradiance_plane)/(max(solar.irradiance_plane)))






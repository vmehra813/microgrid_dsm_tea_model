%% Main Script for Running Generation Asset / Reliability Model: N<10
% Updated 1/7/16, using for thesis.

% Now all in one script, not separating larger n

%% Initialize
clc
clearvars
fclose all;
close all;

tic;

%% Specify Number of Users

% uLink unit: 5 users; 1 generator (A+B box); 4 customers (B box)
num.consumers = 30;
% Note: this corresponds directly to the state-space of available
% solar/battery combinations 

%do not change num.A > 1 for now
num.A = 1;
num.B = num.consumers-num.A;
num.graph = num2str(num.consumers);
num.years = 5;
num.max_consumers = 15;

addpath('/Users/vmehra813/Dropbox (MIT)/[uLink]/Simulation/MatLab/uLink/Modeling Analysis/TPS');
addpath('/Users/ramatya/Dropbox (MIT)/[uLink]/Simulation/MatLab/uLink/Modeling Analysis/TPS');
addpath('D:\Claudio\Dropbox (MIT)\[uLink]\Simulation\matlab\uLink\Modeling Analysis\TPS');


%% Cost Inputs

% config.pv_batt = 'PVBatteryConfigurations.xlsx';
% config.sheet_name = 'Data';
% config.PanelSize = xlsread(config.pv_batt, config.sheet_name, 'PanelSize');
% %config.PVArea = xlsread(config.pv_batt, config.sheet_name,'PVArea');
% config.BattCap = xlsread(config.pv_batt, config.sheet_name, 'BattCap');
% %config.PVEff = xlsread(config.pv_batt, config.sheet_name, 'PVEff');
% %config.system = 8; % 7 is 75 watts // 10 is 250 watts, 180 ah 8 is 125
% %Watts, 180 Ah battery [old version]

config.PVEff = 0.15; % 15 %

% Capital Cost Data
% config.PVCosts = xlsread(config.pv_batt,config.sheet_name,'PVCosts'); %7 is 75W, 10 is 250
% config.BattCosts = xlsread(config.pv_batt,config.sheet_name,'BattCosts'); % 7 is 75Ah, 10 is 180


config.Rs_to_Dollar = 66; %exchange rate

% inputs taken from Jharkhand Trial Summer 2015
%config.WireCost = .28; % $ / meter
config.WireCost = .06; % $ / meter - Source: Wardah
%config.WireCost = 3600/(config.Rs_to_Dollar*90); %cost from Anisha, $/meter, 2.5 mmsq
%config.EthernetCost = 0.42; % $ / meter
config.EthernetCost = 0.03; % $ / meter - Source: Wardah
%config.EthernetCost = 5800/(config.Rs_to_Dollar*305); %cost from Anisha, $/meter cat 6
config.DistrVolt = 24; % volts of network (can be 48V)
config.WireRes = 0.0083; % AWG 14 0.00658*3.28; %ohm / meter (copper)


%% Distances
config.MinDist = 20; % meters
config.Avg_Dist = 40; %40 meters distance per house
config.MaxDist = (num.consumers*config.Avg_Dist-config.MinDist)+config.MinDist; % meters, between min and num.B*avg.dist
% should be num.B


%% Loads and Other Inputs
config.BattVolt = 12; % volts
config.LED_Cost = 1.67; % $ / unit - Source: General price for 5 W LED ~Rs 100.
%config.Fan_Cost = 24.2; % $/unit [original is $24.2 from June 2015 DC Fan]
config.Fan_Cost = 10; % $/unit - Source: Alibaba 500 piece order (12V, 15W)
config.LED_Num = 2; % number of units
config.Fan_Num = 1; % number of units
config.Install_Cost = 2.67; % $ / labor/household for installation (All the way to household wiring)
config.ABox_Cost = 40; % $ / unit
config.BBox_Cost = 15; % $ / unit
config.Pole_Cost = 3.33; % $ / pole - Source: Wardah
config.discount_rate = .12; %world bank reference
config.discount_rate_month = (1+config.discount_rate)^(1/12)-1;
config.payback = 5; % years
%config.Charger = 3; % 200 Rs per charger
config.Charger = 100/config.Rs_to_Dollar; % Source: Alibaba 100 Rs per charger

%have not included in cost function yet!!
config.cnse_led = 24.73 / 1000; % $/wh calculated for LED's
%config.cnse_fan = ;
%config.cnse_cell = ;

%% Timescale / Lifetime of Equipment / Financial Inputs 
config.payback = 5; % years for revenue payback

config.lifetime = 20; %years, lifetime of system
%assume this is always greater than any of equipment lifetimes!

%config.maint_rate = .11; %maintenance cost in percent of initial capital cost (equipment replacment + labor for maintenance)
config.maint_rate = .025; %maintenance cost in percent of initial capital cost (labor for maintenance)

config.Battery_Lifetime = 3; %years (can do by cycle?)
%connect_ulink_customer_VM variables: SOC.total_cycles &
%SOC.critical.cycles
%config.Battery_Cycles = 1500; %number of cycles

config.Poles_Lifetime = 5; %years
config.PV_Lifetime = 15; %years
config.Load_Lifetime = 2; %years
config.Box_Lifetime = 5; %years

config.connection = 0; % 0 connection fee ($)


%% Load Inputs
% NEED TO BE CONSISTENT W/ CREATE DEMAND PROFILES
config.hours_per_day = 4;
config.LED_Watt = 5; % Watts
config.Fan_Watt = 14; % Watts
config.Batt_Days = 2; %days w/ additional storage

config.Charger_Watt = 2.5; % 2.5 watts

config.Batt_Eff = 0.85; %roundtrip efficiency


config.payment = 2; % $ 2 / month payment
config.connection = 5; % $60 / connection fee (JUSCO)
%config.connection = .2; % percentage of capital costs per household
%config.MDOD = 0.6;

% 
% %% Generation Asset Sizing + Area of Panel
% % size panel for charging + lighting + fans
% %size battery for charging + lighting
% config.option = 1; % CHOOSE BETWEEN 1, 2 or 3
% % % 1 is PV Size based on all loads
% % % 2 is PV Size based on evening loads (should be daytime?)
% % % 3 is PV Size as 250 watts
% %[config.PV_Size, config.Batt_Size, config.PV_per_watt] = gen_asset_sizing(config,num);
% %
% % %include logic for battery $/kWh
% %
% % config.PV_area = pv_area(config);
% %
% % if config.Batt_Size > 90
% %     config.Batt_per_kwh = 110; %
% % else
% %     config.Batt_per_kwh = 150;
% % end



%% Solar Panel Generation

solar.solaroutputs=xlsread('pvwatts_hourly_patamda_jharkhand.xlsx','DATA','A20:K8779');
%solar.solaroutputs=xlsread('pvwatts_hourly_karnataka.xlsx','DATA','A20:K8779');

%solar.irradiance_beam = solar.solaroutputs(:,4);    % plane of array irradiance (W/m^2) (beam = 4)
solar.irradiance_plane = solar.solaroutputs(:,8);    % plane of array irradiance (W/m^2) (beam = 4)
solar.tAmbient=solar.solaroutputs(:,6);             % ambient temperature in C
solar.hourofday=solar.solaroutputs(:,3);
solar.month=solar.solaroutputs(:,1);
%solar.temp_coeff = 0.06; % % degradation per degree C above STC // Tata Power Solar
solar.temp_coeff = 0.0041; % literature on mono-si
solar.temp_STC = 25; % 25 degrees C reference temperature

% Charge Controller

solar.cc_data=xlsread('PowerElectronics_Eff','ChargeController','A2:B17');
solar.cc_power=solar.cc_data(:,1);
solar.cc_eff=solar.cc_data(:,2);
solar.cc_fcn=polyfit(solar.cc_power,solar.cc_eff,4);

%visual plot of distances / nodes

% Call Function
%[a] = singlehome_generation_VM(config, solar,num);
%a.available_solar is output

% figure;
% plot(solar.irradiance_plane);
% xlabel('Hour of Year');
% ylabel('Irradiance (W/m^2)');
% set(gca, 'fontsize',14);



%% Creating Demand Profiles

%demand.inputs = 'demandInputs_uLink_gen_reliability.xlsx';
%demand.inputs = 'demandInputs_uLink_gen_reliability_thesis.xlsx';
demand.inputs = 'demandInputs_uLink_gen_reliability_MGP.xlsx';


demand.type.residential = 'Residential';
demand.type.anchor = 'Anchor';
%two LED's, one fan, one cell phone charger

demand.distr_loss = 0.9;

%incorporating losses
demand.min_distance = config.MinDist; % meters
demand.max_distance = config.MaxDist; % meters
demand.wire_resistivity = config.WireRes; %Ohm / foot * foot/meter
% NEED WIRE SIZE
%demand.distances_b = demand.wire_resistivity*((demand.max_distance - demand.min_distance)*rand(num.consumers,1)+demand.min_distance); %vector of distances for b1-bn in a radial format from A box
demand.distances = randfixedsum(num.consumers,1,demand.max_distance,demand.min_distance,demand.max_distance);

%demand.distances = config.MaxDist;
demand.distances_b = demand.wire_resistivity*demand.distances;

demand.distr_voltage = config.DistrVolt;
demand.CNSE = 0.01;

% %Plotting for thesis (one household)
% figure;
% subplot(2,1,1),plot(loadprof.final.total(:,1),'c');
% xlabel('Hour of Year');
% ylabel('Watts');
% set(gca,'fontsize',14);
% legend('Total Load (1 House)');
% subplot(2,1,2),plot(loadprof.final.critical(:,1),'r');
% xlabel('Hour of Year');
% ylabel('Watts');
% set(gca,'fontsize',14);
% legend('Critical Load (1 House)');

% %Introduce Converter Losses
%can do effects of polynomial fitting / loss function / regularization

%Load Converter
demand.load_data =xlsread('PowerElectronics_Eff','LoadConverter','A2:B11');
demand.load_power =demand.load_data(:,1);
demand.load_eff=demand.load_data(:,2)/100;
demand.load_fcn = polyfit(demand.load_power,demand.load_eff,8);



%Link Converter
demand.link_data=xlsread('PowerElectronics_Eff','LinkConverter','A2:B15');
demand.link_power=demand.link_data(:,1);
demand.link_eff=demand.link_data(:,2)/100;
demand.link_fcn=polyfit(demand.link_power,demand.link_eff,4);
%loss calculations are in create_demand_profiles

%Call Function
[loadprof, load] = CreateDemandProfiles_1016_VM(num,demand, solar);

num.total_hour_threshold = 0.9;
num.critical_hour_threshold = 0.99;

%% Battery Inputs
global battery;

battery.charge_eff = .95;          % Battery charging efficiency
battery.discharge_eff = .95;       % Battery discharging efficieny
battery.voltage = 12;                                                         % Battery operating voltage [V]
%battery.size = config.Batt_Size*battery.voltage; % (7) is 75 Ah; (10) is 180 Ah
%battery.initial_capacity = battery.size / 1; %initial capacity


%battery.MDOD = .6; %max depth of discharge, 60% (can change for SOC inaccuracy)
battery.MLOC = 1; %max level of charge

%battery.MDOD = .3;
%battery.MDOD = .4;
%battery.MDOD = .5;
battery.MDOD = .6;
%battery.MDOD = .7;
%battery.MDOD = .8;

% if battery.MDOD == .3
%     battery.num_cycles = 1500; %lead acid cycle life
% elseif battery.MDOD == .4
%     battery.num_cycles = 1200; %lead acid cycle life
% elseif battery.MDOD == .5
%     
% elseif battery.MDOD == .6
%     battery.num_cycles = 1200; %lead acid cycle life    
% elseif battery.MDOD == .7
%     battery.num_cycles = 1000; %lead acid cycle life
% elseif battery.MDOD == .8
%     battery.num_cycles = 800; %lead acid cycle life
% end
% 
% 
% battery.deg_per_cycle = 0.01; %degradation % per cycle
voltage = battery.voltage;

%% Exponential Equation for Cycle Life versus MDOD 

battery.num_cycles = (5891*exp(-2.382*battery.MDOD));
battery.num_cycles = round(battery.num_cycles,-2);

config.Battery_Cycles = battery.num_cycles;


%% KIBAM Inputs

%Battery selection: VRLA Panasonic Example
%Trojan example: 5 hour: 185Ah, 10hr: 207 Ah; 20 hr: 225 Ah
battery.t = 1;
battery.t1 = 10;
battery.t2 = 20; %hours

battery.t_cap = 4.5; %Ah
%battery.t_cap = 185;
battery.t1_cap = 6.8; %Ah
%battery.t1_cap = 207;
battery.t2_cap = 7.2; %Ah
%battery.t2_cap = 225;

% battery.F1 = battery.t_cap / battery.t1_cap;
% battery.F2 = battery.t_cap / battery.t2_cap;

battery.deltaT = 1;
battery.alpha_c = 1;
battery.Imax = 12;

% default battery c and k values 

battery.c = .305;
battery.k = 2.12;


% disp(battery.c_1)
% disp(battery.k_1)
% disp(battery.c)
% disp(battery.k);



%% Iterate on Generation Sizing
%read in spreadsheet to get vector of sizes / ranges

%  [50,60,80,100,130,150,160,210,240,250,260,280,290,300];
%max = 5kW
%combine with integer values
%from tata power solar
config.PV_modules = [50; 60; 80; 100; 130; 150; 160; 210; 240; 250; 260; 280; 290; 300]; %atts

config.PV_modules = [100; 130; 150; 210; 250; 280; 300];

config.PV_modules = [50; 75; 80; 100; 120; 125; 150; 200; 250; 300];
%SU KAM FROM ANISHA

config.PV_modules = [150; 200; 250; 300; 400; 600; 800; 1000; 1200;];

%config.PV_modules = [200; 400; 600; 800; 1000; 1200; 1400; 1600; 1800; 2000;];


%config.PV_modules = [200];
%back of the envelope calculation
%config.PV_max = roundn(config.LED_Num*config.LED_Watt*num.max_consumers ...
%+ config.Charger_Watt*num.max_consumers,1) ;

%from HBL batteries india
%config.Batt_sizes = [60; 80; 100; 135; 150; 180]; %amp-hours, 20 hr rate
% %config.Batt_sizes_t2 = [10; 20; 60; 80; 100; 135; 150; 180]; %amp-hours 20
% %hr rate
% config.Batt_sizes_t2 = [80; 100; 135; 150; 180]; %20 hr rate
% % EXIDE
% %config.Batt_sizes_t2 = [40; 60; 75; 100; 120; 150; 200; 300]; %10hr rate
%
% config.Batt_sizes_t1 = config.Batt_sizes_t2 / 1.05; % 10 hr rate
%
% config.Batt_sizes = config.Batt_sizes_t2 / 1.5; %1 hour rate
%config.Batt_sizes = [150];

%back of the envelope for max
%config.Batt_max  = roundn((config.Batt_Days*config.hours_per_day*(config.LED_Num*config.LED_Watt*num.max_consumers ...
%    + config.Charger_Watt*num.max_consumers) / ...
%   (config.MDOD*config.Batt_Eff))/config.BattVolt,1);

% Amara Raja 
config.Batt_sizes = [60; 75; 100; 120; 150; 165; 180; 200];

config.Batt_sizes = [120; 150; 165; 180; 200; 300; 400; 500; 600; 700;];
%config.Batt_sizes = [200; 300; 400; 500; 600; 700; 800; 900; 1000; 1100; 1200; 1300; 1400; 1500;];


increment_pv = 200; %pv increment
increment_batt = 100; %ah increment

%thesis example: 250 and 100 ah is optimal for 5 hh 
%config.PV_modules = 250;
%config.Batt_sizes = 100;




%% CHOICE OF PV AND BATTERY RANGES

config.PV_Options = config.PV_modules;
config.Battery_Options = config.Batt_sizes;

%% Initialize Matrices for Storing Values 

reliability_critical_matrix = zeros(length(config.Battery_Options),...
    length(config.PV_Options));

reliability_total_matrix = zeros(length(config.Battery_Options),...
    length(config.PV_Options));

cost_pv_batt_matrix = zeros(length(config.Battery_Options),...
    length(config.PV_Options));

critical_amount_served = zeros(length(config.Battery_Options),...
    length(config.PV_Options));

total_amount_served_mat = zeros(length(config.Battery_Options),...
    length(config.PV_Options));

critical_amt_hours = zeros(length(config.Battery_Options),...
    length(config.PV_Options));

total_amt_hours = zeros(length(config.Battery_Options),...
    length(config.PV_Options));

cycles_critical = zeros(length(config.Battery_Options),...
    length(config.PV_Options));

cycles_total = zeros(length(config.Battery_Options),...
    length(config.PV_Options));


ratio_matrix_pv = zeros(length(config.Battery_Options),...
    length(config.PV_Options));


coin_fan = zeros(length(config.Battery_Options),...
    length(config.PV_Options));

coin_lights = zeros(length(config.Battery_Options),...
    length(config.PV_Options));

coin_cellphone = zeros(length(config.Battery_Options),...
    length(config.PV_Options));


pv_per_watt = zeros(length(config.PV_Options),1);
batt_per_kwh = zeros(length(config.Battery_Options),1);

pv_total_cost = zeros(length(config.PV_Options),1);
pv_number_final = zeros(length(config.PV_Options),1);
pv_module_final = zeros(length(config.PV_Options),1);


batt_total_cost = zeros(length(config.Battery_Options),1);
batt_number_final = zeros(length(config.Battery_Options),1);
batt_capacity_final = zeros(length(config.Battery_Options),1);


%% LEFT OFF HERE (leave as n<10 for now)
%should do same logic on number of panels / batteries, just answer will be
%1 for smaller n
%batt_number_final and pv_number_final, for num.consumers > 10


%% Develop Matrix Entries for Cost and Reliability Matrices
%Still have num.consumers qualifier for > or < n=10

for i=1:length(config.Battery_Options)
    
    batt = config.Battery_Options(i);
    
    if batt <=200
    nbatt = 1;
    batt_per_kwh(i) = cost_batt_per_kwh(batt,voltage);
    batt_size = batt*voltage;
    batt_cost = batt_per_kwh(i)*batt_size/1000;

    else
    nbatt = batt / increment_batt;
    batt_per_kwh(i) = cost_batt_per_kwh(increment_batt,voltage);
    batt_size = batt*voltage;
    batt_cost = batt_per_kwh(i)*increment_batt*nbatt*voltage/1000;
    end
    %batt = 180;
    %batt_per_kwh(i) = cost_batt_per_kwh(batt,voltage);
    
    %battery.F1 = config.Batt_sizes(i) / config.Batt_sizes_t1(i);  %1 hr cap / 10 hr cap
    %battery.F2 = config.Batt_sizes(i) / config.Batt_sizes_t2(i);  %1 hr cap / 20 hr cap
    
    %%determine c/k constants
    % qmax_battbank is calculated in kibam dispatch
    %fun = @kibam_c_k;
    %x0 = [0.5, 0.5];
    %x = fsolve(fun,x0);
    
    %battery.c_1 = x(1);
    %battery.k_1 = x(2);
    
    
    for j=1:length(config.PV_Options)
        
        pv = config.PV_Options(j);
        %pv = 250;
        batt_size = batt*voltage;
        
        if pv <= 80
            config.PVeff = .135;
        else
            config.PVeff = .15;
        end
        
        available_solar = singlehome_gen_sizing(config,solar,num,pv);
        %simple battery model
        %[reliability_critical_matrix(i,j),critical_amount_served(i,j), reliability_total_matrix(i,j), total_amount_served(i,j)] = connect_ulink_customer_gen_sizing_all_load(available_solar,battery,loadprof,num,batt_size);
        %kinetic battery model
        %[reliability_critical_matrix(i,j), critical_amount_served(i,j),reliability_total_matrix(i,j), total_amount_served(i,j)] = connect_ulink_customer_gen_sizing_kibam_pb_acid(available_solar,battery,loadprof,num,batt);
        %kibam with dispatch/curtailment
        %[reliability_critical_matrix(i,j), critical_amount_served(i,j),reliability_total_matrix(i,j), total_amount_served(i,j), cycles] = connect_ulink_customer_gen_sizing_kibam_load_dispatch(available_solar,battery,loadprof,num,batt,load);
        %kibam with dispatch/curtailment + SOC ranges (green / yellow /
        %red)
        %kibam with only vector (not doing per consumer.. only calculating
        %for nth consumer.. no matrix issues / floating point operations
        %minimized / memory required is less.
        %[reliability_critical_matrix(i,j), critical_amount_served(i,j),reliability_total_matrix(i,j), total_amount_served_mat(i,j), cycles] = connect_ulink_customer_gen_sizing_kibam_load_dispatch_8760x1(available_solar,battery,loadprof,num,batt,load);
        
        % new definition of reliability (hour count)
        loadprof1 = loadprof.final.critical;
        [reliability_critical_matrix(i,j), critical_amt_hours(i,j),cycles] = connect_ulink_customer_gen_sizing_kibam_load_dispatch_SAIDI(available_solar,battery,loadprof1,num,batt,load);
        loadprof2 = loadprof.final.total;
        [reliability_total_matrix(i,j), total_amt_hours(i,j), cycles] = connect_ulink_customer_gen_sizing_kibam_load_dispatch_SAIDI(available_solar,battery,loadprof2,num,batt,load);

        
        %cycles_critical(i,j) = cycles.final_critical; %count based on ~daily cycle (versus hourly)
        %cycles_total(i,j) = cycles.final_total;
        
        %solar / load coincidence factor
        %[coincidence] = solar_load_factor(available_solar,loadprof);
        %coin_fan(i,j) = coincidence.fan_avg;
        %coin_lights(i,j) = coincidence.lights_avg;
        %coin_cellphone(i,j) = coincidence.cell_avg;
        
        
        %pv_per_watt(j) = cost_pv_per_watt(pv);
        
        %total cost of pv/batt combination
        %cost_pv_batt_matrix(i,j) = batt_per_kwh(i)*batt_size/1000 + pv_per_watt(j)*pv ;
        %ratio_matrix_pv(i,j) = pv_per_watt(j)*pv / cost_pv_batt_matrix(i,j);
        
        if pv <=300
        pv_per_watt(j) = cost_pv_per_watt(pv);
        pv_cost = pv_per_watt(j)*pv;
        num_pv = 1;
        else
        num_pv = pv / increment_pv;
        pv_per_watt(j) = cost_pv_per_watt(increment_pv);
        pv_cost = pv_per_watt(j)*increment_pv*num_pv;
        end
        
        %total cost of pv/batt combination
        cost_pv_batt_matrix(i,j) =  batt_cost + pv_cost;
        
        
        
        
    end
    
    
end

% figure;
% pcolor(config.Battery_Options,config.PV_Options,ratio_matrix_pv');
% colorbar;
% xlabel('Battery Capacity (Ah)');
% ylabel('Solar Module (W)');
% legend('Proportion of Solar of Generation Asset Costs (%)');
% set(gca,'fontsize',15);



%% Check Reliability Thresholds

% % rel_crit = [.8, .9, .92, .93, .94, .95, .96, .97, .99];
% % rel_tot = [.8, .9, .92, .93, .94, .95, .96, .97];
% % 
% % dol_watt = zeros(length(rel_crit),length(rel_tot));
% % sol_batt = zeros(length(rel_crit),length(rel_tot));
% % lcoe = zeros(length(rel_crit),length(rel_tot));
% % sol = zeros(length(rel_crit),length(rel_tot));
% % bat = zeros(length(rel_crit),length(rel_tot));
% % 
% % for i=1:length(rel_crit)
% % 
% %     for j=1:length(rel_tot)

%%  RELIABILITY BASED ON NUMBER OF HOURS
reliability.critical = num.critical_hour_threshold; %night loads (lights)

reliability.total = num.total_hour_threshold; %incl. day loads (fan)

%% RELIABILITY BASED ON HOURS (CHANGE DECISION IN MATRIX INDEX SELECTION BELOW
reliability.critical_hours = 100;  % NUMBER OF HOURS
reliability.total_hours = 500;
%%


reliability.range_critical = 0; % - range of thresholds
reliability.range_total = 0;

%find indices of matrix that meet total reliability threshold
[ind.total_i,ind.total_j] = find(reliability_total_matrix > (reliability.total-reliability.range_total) & ...
    reliability_total_matrix < (reliability.total+reliability.range_total));

%turn into a matrix
%ind.total = [ind.total_i,ind.total_j];

%get subscripts of every matrix entry
ind.total_sub = sub2ind(size(reliability_total_matrix),ind.total_i,ind.total_j);

if isempty(ind.total_sub)
    [ind.total_i,ind.total_j] = find(reliability_total_matrix > (reliability.total-reliability.range_total));
    ind.total_sub = sub2ind(size(reliability_total_matrix),ind.total_i,ind.total_j);
    
    if isempty(ind.total_sub)
        error('No PV/Batt Combinations Meet Total Reliability Threshold. Increase Available State-Space and Review Total Reliability Values to Check for Memory/Stability Issues.');
    end
end



%find indices of matrix that meet critical reliability threshold
[ind.critical_i,ind.critical_j] = find(reliability_critical_matrix > (reliability.critical-reliability.range_critical) & ...
    reliability_critical_matrix < (reliability.critical+reliability.range_critical));

%ind.critical = [ind.critical_i,ind.critical_j];

%get subscripts of every matrix entry
ind.critical_sub = sub2ind(size(reliability_critical_matrix),ind.critical_i,ind.critical_j);

if isempty(ind.critical_sub)
    [ind.critical_i,ind.critical_j] = find(reliability_critical_matrix > (reliability.critical-reliability.range_critical));
    ind.critical_sub = sub2ind(size(reliability_critical_matrix),ind.critical_i,ind.critical_j);
    
    if isempty(ind.critical_sub)
        error('No PV/Batt Combinations Meet Critical Reliability Threshold. Increase Available State-Space and Review Critical Reliability Values to Check for Memory/Stability Issues.');
    end
end



% Overlap of Subscripts / Indices
%if there are more of one element versus the other
ind.cost_sub = [];

if length(ind.critical_sub) <= length(ind.total_sub)
    ind.cost_sub = find(ismember(ind.critical_sub,ind.total_sub));
    ind.cost_sub_use = ind.critical_sub(ind.cost_sub);
    
else
    ind.cost_sub = find(ismember(ind.total_sub,ind.critical_sub));
    ind.cost_sub_use = ind.total_sub(ind.cost_sub);
   
end
%get indices back

%find lowest cost combination
[ind.cost_value, ind.cost_index] = min(cost_pv_batt_matrix(ind.cost_sub_use));



ind.cost_index_use = ind.cost_sub_use(ind.cost_index);

%[ind.cost_i,ind.cost_j] = ind2sub(size(reliability_critical_matrix),ind.cost_sub);
%find combination of battery and pv
[ind.batt_final,ind.pv_final]= ind2sub(size(cost_pv_batt_matrix),ind.cost_index_use);

[tps.batt_final,tps.pv_final]= ind2sub(size(cost_pv_batt_matrix),ind.cost_sub_use);

a=[config.Battery_Options(tps.batt_final), config.PV_Options(tps.pv_final)];
b = [batt_per_kwh(tps.batt_final),pv_per_watt(tps.pv_final)];
c= a(:,1).*b(:,1)*12/1000 + a(:,2).*b(:,2);

a(:,3) = c;

%d=eye(length(a(:,1)),length(a(:,2))).*a(:,3);

%% ORIGINAL -- 
% ind.pv_final = 9;
% ind.batt_final = 3;

config.pv_final = config.PV_Options(ind.pv_final);
config.batt_final = config.Battery_Options(ind.batt_final);
critical_final = reliability_critical_matrix(ind.batt_final,ind.pv_final);
total_final = reliability_total_matrix(ind.batt_final,ind.pv_final);


%total_amount_served = critical_amount_served(ind.batt_final,ind.pv_final)+total_amount_served(ind.batt_final,ind.pv_final);
total_amount_served = total_amount_served_mat(ind.batt_final,ind.pv_final);


config.pv_per_watt = pv_per_watt(ind.pv_final);
config.batt_final_per_kwh = batt_per_kwh(ind.batt_final);

config.cycles_final_critical = ceil(cycles_critical(ind.batt_final,ind.pv_final)); %unused value right now. 
config.cycles_final_total = ceil(cycles_total(ind.batt_final,ind.pv_final)); %need integer value


% config.cycles_final_critical = 365;
% config.cycles_final_total = 365; 


% Cost, Revenue, and DCF Analysis
% note that cost is being updated to include replacement costs
%adding system lifetime

%careful with costs on battery/pv here
%need conditional in cost function

%% THESIS CHANGES -- SENSITIVITY ANALYSIS 

% conn = [0 5, 10, 15, 20];
% disc = [.08 .10 .12 .14 .16 .18 .20];
% payb = [2, 3, 5, 7, 10];
% ir = zeros(length(conn),length(disc));
% lco = zeros(length(conn),length(disc));
% mon = zeros(length(conn),length(disc));
% 
% ir = zeros(length(payb),length(disc));
% lco = zeros(length(payb),length(disc));
% mon = zeros(length(payb),length(disc));
% 
% for i=1:length(payb)
%     
%     for j=1:length(disc)
% config.discount_rate = disc(j); %world bank reference
% config.discount_rate_month = (1+config.discount_rate)^(1/12)-1;
% config.payback = payb(i); % years
% 
% config.connection = 0;
% 

% 
% ir(i,j) = revenue.irr;
% lco(i,j) = cost.lcoe;
% mon(i,j) = revenue.breakeven_monthly_connection;
% 
%     end
% end


% thesis figure
% cost.thesis_plot = [[cost.npv],[revenue.vector]];
% figure;
% plot(cost.thesis_plot);
% legend('Cost: Initial and Yearly','Revenue: Connection and Yearly','location','SouthEast');
% set(gca,'fontsize',15);
% xlabel('Years');
% ylabel('Dollars ($)');
% 
% dol_watt(i,j) = cost.initial_per_watt_ulink;
% sol_batt(i,j) = ind.cost_value;
% lcoe(i,j) = cost.lcoe;
% sol(i,j) = config.pv_final;
% bat(i,j) = config.batt_final;
% 
% 
% 
%     end
% end

% figure;
% subplot(3,1,1);
% a=contour(rel_crit,rel_tot,lcoe');
% clabel(a,'manual','FontSize',15);
% set(gca,'fontsize',15);
% xlabel('Critical Reliability (%)');
% ylabel('Total Reliability (%)');
% colorbar;
% legend('LCOE ($/kWh)');
% subplot(3,1,2);
% a=contour(rel_crit,rel_tot,dol_watt');
% clabel(a,'manual','FontSize',15);
% set(gca,'fontsize',15);
% xlabel('Critical Reliability (%)');
% ylabel('Total Reliability (%)');
% colorbar;
% legend('Capital Cost ($/W)');
% subplot(3,1,3);
% a=contour(rel_crit,rel_tot,sol_batt');
% clabel(a,'manual','FontSize',15);
% set(gca,'fontsize',15);
% xlabel('Critical Reliability (%)');
% ylabel('Total Reliability (%)');
% colorbar;
% legend('Solar and Battery Cost ($)');
% 
% figure;
% subplot(2,1,1);
% a=contourf(rel_crit,rel_tot,sol');
% clabel(a,'manual','FontSize',15);
% set(gca,'fontsize',15);
% xlabel('Critical Reliability (%)');
% ylabel('Total Reliability (%)');
% colorbar;
% legend('Solar Module (W)');
% 
% subplot(2,1,2);
% a=contourf(rel_crit,rel_tot,bat');
% clabel(a,'manual','FontSize',15);
% set(gca,'fontsize',15);
% xlabel('Critical Reliability (%)');
% ylabel('Total Reliability (%)');
% colorbar;
% legend('Battery Capacity (Ah)');

%% ORIGINAL 

[cost, revenue] = cost_uLink_csr_gen_reliability(num,demand,config, total_amount_served,voltage);
 network_cost = cost.wiring + cost.ethernet;
elapsed_time = toc;

%% OUTPUTS 
disp(strcat({'LCOE Over 5 Years is (in $/kWh) '},{' '},{num2str(cost.lcoe)}));
disp(strcat({'Breakeven Monthly (w/ No Fee)'},{' $'},{num2str(revenue.breakeven_monthly)}));
disp(strcat({'Breakeven Monthly (w/ Fee)'},{' $'},{num2str(revenue.breakeven_monthly_connection)}));

%  Display Results 
disp(strcat({'Final PV Size Watts '},{' '},{num2str(config.pv_final)}));
disp(strcat({'Final Batt Size Ah '},{' '},{num2str(config.batt_final)}));
%disp(strcat({'Final Cost $ of PV/Batt Combo '},{' $'},{num2str(ind.cost_value)}));


disp(strcat({'Final Critical Reliability % '},{' '},{num2str(critical_final)}));
disp(strcat({'Final Total Reliability % '},{' '},{num2str(total_final)}));

%disp(strcat({'Amount of Time to Run (Seconds)'},{num2str(elapsed_time)}));

% Display Values / Plot
%disp(strcat({'Average Yearly Reliability serving Total Load is '},{num2str(network.total.reliability(end))}));
%disp(strcat({'Wiring/Ethernet Cost '},{' $'},{num2str(network_cost)}));
%disp(strcat({'App Costs (Fan/LED), Not Incl. in LCOE'},{' $'},{num2str(cost.lighting+cost.fan)}));

disp(strcat({'Cap Cost $/W '},{' '},{num2str(cost.initial_per_watt_ulink)}));
disp(strcat({'NPV Based on Costs '},{' $'},{num2str(cost.npv_total)}));
disp(strcat({'Total Capital Cost '},{' $'},{num2str(cost.initial)}));
%disp(strcat({'PV Size is (Watts) '},{num2str(config.PV_Size)}));
%disp(strcat({'Battery Size is (Ah) '},{num2str(config.Batt_Size)}));

%config.cost_string = {strcat('Solar',{' $'},num2str(cost.solar_panel)),strcat('Battery',{' $'},num2str(cost.battery)),strcat('Wiring',{' $'},num2str(cost.wiring)),strcat('Ethernet',{' $'},num2str(cost.ethernet)),...
%    strcat('Poles',{' $'},num2str(cost.poles)),strcat('Installation',{' $'},num2str(cost.installation)),strcat('Appliances',{' $'},num2str(cost.fan+cost.lighting)),strcat('Network Devices',{' $'},num2str(cost.ABox+cost.BBox))};
config.cost_string = {strcat('Solar'),strcat('Battery'),strcat('Wiring'),strcat('Ethernet'),...
    strcat('Poles'),strcat('Installation'),strcat('Appliances'),strcat('Network Devices')};


config.explode = ones(length(config.cost_string),1);


%% Figures
figure(15);
pie([cost.solar_panel,cost.battery,cost.wiring,cost.ethernet,...
    cost.poles,cost.installation,cost.fan+cost.lighting,cost.ABox+cost.BBox],config.explode);
set(gca,'fontsize',15);
legend(config.cost_string,'Location','southoutside','Orientation','horizontal');

title(strcat('Distribution of Capital Cost of',{' $'},num2str(cost.initial_ulink)));


% figure(16);
% pie([cost.initial_ulink,(-cost.maint_replace_total/config.lifetime)],{strcat('Capital',{' $'},num2str(cost.initial_ulink)),strcat('Operation/Main',{' $'},num2str(-cost.maint_replace_total/config.lifetime))});
% title('Capital and Operating Cost (Yearly) Breakdown in %');
% set(gca,'fontsize',15);


 %% Matrix Outputs

% figure(17);
% %imagesc(fliplr(config.Battery_Options),config.PV_Options,reliability_critical_matrix);
% surf(reliability_critical_matrix);
% xlabel('Solar Module Options (Watts)');
% set(gca,'XTickLabel',{'50','100','125','150','200','300'});
% %set(gca,'XTick',{'50','75','80','100','120','125','150','200','250','300'})
% set(gca,'YTickLabel',{'60','100','120','165','200'});
% ylabel('Battery Capacity Options (Amp-Hours)');
% set(gca,'fontsize',15);
% zlabel('Reliability');
% 
% 
% figure(18);
% %imagesc(fliplr(config.Battery_Options),config.PV_Options,reliability_critical_matrix);
% surf(reliability_total_matrix);
% xlabel('Solar Module Options (Watts)');
% set(gca,'XTickLabel',{'50','100','125','150','200','300'});
% set(gca,'YTickLabel',{'60','100','120','165','200'});
% ylabel('Battery Capacity Options (Amp-Hours)');
% set(gca,'fontsize',15);
% zlabel('Reliability');
% 
% 
% figure(19);
% %imagesc(fliplr(config.Battery_Options),config.PV_Options,reliability_critical_matrix);
% surf(cost_pv_batt_matrix);
% xlabel('Solar Module Options (Watts)');
% set(gca,'XTickLabel',{'50','100','125','150','200','300'});
% set(gca,'YTickLabel',{'60','100','120','165','200'});
% ylabel('Battery Capacity Options (Amp-Hours)');
% set(gca,'fontsize',15);
% zlabel('Generation Asset Costs ($)');


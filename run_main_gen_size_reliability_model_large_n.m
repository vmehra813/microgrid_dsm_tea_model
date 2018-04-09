%% Main Script for Running Generation Asset / Reliability Model: N>10
%Larger Scale Version
% Updated 3/29/16

%% Initialize
clc
clearvars
fclose all;
close all;

%% Specify Number of Users

% uLink unit: 5 users; 1 generator (A+B box); 4 customers (B box)
%threshold for larger configurations is n=10 for now.
%note that load profiles are kept the same!
num.max_consumers = 100;

num.consumers = 20;
%should only enter n>10 here due to changes made in code!
num.A = 1;
num.B = num.consumers-num.A;
num.graph = num2str(num.consumers);
num.years = 5;

addpath('/Users/vmehra813/Dropbox (MIT)/[uLink]/Simulation/MatLab/uLink/Modeling Analysis/TPS');



%% PV and Battery Input Data: Costs

config.pv_batt = 'PVBatteryConfigurations.xlsx';
config.sheet_name = 'Data';
config.PanelSize = xlsread(config.pv_batt, config.sheet_name, 'PanelSize');
%config.PVArea = xlsread(config.pv_batt, config.sheet_name,'PVArea');
config.BattCap = xlsread(config.pv_batt, config.sheet_name, 'BattCap');
%config.PVEff = xlsread(config.pv_batt, config.sheet_name, 'PVEff');
%config.system = 8; % 7 is 75 watts // 10 is 250 watts, 180 ah 8 is 125
%Watts, 180 Ah battery [old version]

config.PVEff = 0.15; % 15 %

% Capital Cost Data
config.PVCosts = xlsread(config.pv_batt,config.sheet_name,'PVCosts'); %7 is 75W, 10 is 250
config.BattCosts = xlsread(config.pv_batt,config.sheet_name,'BattCosts'); % 7 is 75Ah, 10 is 180

% inputs taken from Jharkhand Trial Summer 2015
%config.WireCost = .28; % $ / meter
config.WireCost = .06; % $ / meter - Source: Wardah
%config.EthernetCost = 0.42; % $ / meter
config.EthernetCost = 0.03; % $ / meter - Source: Wardah
config.DistrVolt = 24; % volts
config.WireRes = 0.00658*3.28; %ohm / meter

config.MinDist = 10; % meters
config.Avg_Dist = 40; %40 meters distance per house 
config.MaxDist = (num.B*config.Avg_Dist-config.MinDist)+config.MinDist; % meters, between min and num.B*avg.dist


config.BattVolt = 12; % volts
config.LED_Cost = 1.67; % $ / unit - Source: General price for 5 W LED ~Rs 100.
%config.Fan_Cost = 24.2; % $/unit [original is $24.2 from June 2015 DC Fan]
config.Fan_Cost = 10; % $/unit - Source: Alibaba 500 piece order (12V, 15W)
config.LED_Num = 2; % number of units
config.Fan_Num = 1; % number of units
config.Install_Cost = 2.67; % $ / labor/household for installation (All the way to household wiring)
config.ABox_Cost = 50; % $ / unit
config.BBox_Cost = 20; % $ / unit
config.Pole_Cost = 3.33; % $ / pole - Source: Wardah
config.discount_rate = .18; %18 % discount rate in India
config.payback = 5; % years
%config.Charger = 3; % 200 Rs per charger
config.Charger = 1.5; % Source: Alibaba 100 Rs per charger



%% Timescale / Lifetime of Equipment
config.payback = 5; % years for revenue payback

config.lifetime = 20; %years, lifetime of system
%assume this is always greater than any of equipment lifetimes!

%config.maint_rate = .11; %maintenance cost in percent of initial capital cost (equipment replacment + labor for maintenance)
config.maint_rate = .025; %maintenance cost in percent of initial capital cost (labor for maintenance)

config.Battery_Lifetime = 3; %years (can do by cycle?)
%connect_ulink_customer_VM variables: SOC.total_cycles &
%SOC.critical.cycles
config.Battery_Cycles = 1500; %number of cycles

config.Poles_Lifetime = 5; %years
config.PV_Lifetime = 15; %years
config.Load_Lifetime = 2; %years
config.Box_Lifetime = 5; %years



%% Load Inputs
% NEED TO BE CONSISTENT W/ CREATE DEMAND PROFILES
config.hours_per_day = 4;
config.LED_Watt = 5; % Watts
config.Fan_Watt = 14; % Watts
config.Batt_Days = 2; %days w/ additional storage
config.MDOD = 0.5; %max discharge (Lead acid battery)
config.Charger_Watt = 2.5; % 2.5 watts

config.Batt_Eff = 0.85; %roundtrip efficiency


config.payment = 2; % $ 2 / month payment
config.connection = 30; % $60 / connection fee (JUSCO)


%% Generation Asset Sizing + Area of Panel
% size panel for charging + lighting + fans
%size battery for charging + lighting
% config.option = 1; % CHOOSE BETWEEN 1, 2 or 3
% % 1 is PV Size based on all loads
% % 2 is PV Size based on evening loads (should be daytime?)
% % 3 is PV Size as 250 watts
% [config.PV_Size, config.Batt_Size, config.PV_per_watt] = gen_asset_sizing(config,num);
% 
% %include logic for battery $/kWh
% 
% config.PV_area = pv_area(config);
% 
% if config.Batt_Size > 90
%     config.Batt_per_kwh = 110; %
% else
%     config.Batt_per_kwh = 150;
% end



%% Solar Panel Generation

solar.solaroutputs=xlsread('pvwatts_hourly_patamda_jharkhand.xlsx','DATA','A20:K8779');
%solar.irradiance_beam = solar.solaroutputs(:,4);    % plane of array irradiance (W/m^2) (beam = 4)
solar.irradiance_plane = solar.solaroutputs(:,8);    % plane of array irradiance (W/m^2) (beam = 4)
solar.tAmbient=solar.solaroutputs(:,6);             % ambient temperature in C
solar.hourofday=solar.solaroutputs(:,3);
solar.month=solar.solaroutputs(:,1);
%solar.temp_coeff = 0.06; % % degradation per degree C above STC // Tata Power Solar
solar.temp_coeff = 0.0041; % literature on mono-si
solar.temp_STC = 25; % 25 degrees C reference temperature

%visual plot of distances / nodes

% Call Function
%[a] = singlehome_generation_VM(config, solar,num);
%a.available_solar is output


%% Creating Demand Profiles

demand.inputs = 'demandInputs_uLink_gen_reliability.xlsx';
demand.type.residential = 'Residential';
demand.type.anchor = 'Anchor';
%two LED's, one fan, one cell phone charger

demand.distr_loss = 0.9;

%incorporating losses
demand.min_distance = config.MinDist; % meters
demand.max_distance = config.MaxDist; % meters
demand.wire_resistivity = 0.00659; %Ohm / foot
% NEED WIRE SIZE
%demand.distances_b = demand.wire_resistivity*((demand.max_distance - demand.min_distance)*rand(num.consumers,1)+demand.min_distance); %vector of resistances over length for b1-bn in a radial format from A box
demand.distances = randfixedsum(num.consumers,1,demand.max_distance,demand.min_distance,demand.max_distance);
demand.distances_b = demand.wire_resistivity*demand.distances;

demand.distr_voltage = config.DistrVolt;
demand.CNSE = 0.01;

%Call Function
[loadprof, load] = CreateDemandProfiles_VM(num,demand, solar)




%% Battery Balance

battery.charge_eff = 0.95;                                                      % Battery charging efficiency
battery.discharge_eff = 0.95;                                                   % Battery discharging efficieny
battery.voltage = 12;                                                         % Battery operating voltage [V]
%battery.size = config.Batt_Size*battery.voltage; % (7) is 75 Ah; (10) is 180 Ah
%battery.initial_capacity = battery.size / 1; %initial capacity
battery.MDOD = config.MDOD;

% converter efficiencies: not included yet
battery.a_dc_dc_converter_eff = 0.95; %stepping up voltage to 24 volts from ~12 at battery, A Box
battery.b_dc_dc_converter_eff_1 = 0.95; %can these be voltage ratios
battery.b_dc_dc_converter_eff_2 = 0.95;
battery.b_dc_dc_converter_eff_3 = 0.95;
battery.b_dc_dc_converter_eff_4 = 0.95;
battery.b_dc_dc_converter_eff_5 = 0.95;

voltage = battery.voltage;


%% Iterate on Generation Sizing
%read in spreadsheet to get vector of sizes / ranges

%  [50,60,80,100,130,150,160,210,240,250,260,280,290,300];
%max = 5kW
%combine with integer values 
%from tata power solar
config.PV_modules = [50; 60; 80; 100; 130; 150; 160; 210; 240; 250; 260; 280; 290; 300]; %atts
%back of the envelope calculation
config.PV_max = roundn(config.LED_Num*config.LED_Watt*num.max_consumers ...
     + config.Charger_Watt*num.max_consumers,1) ;
 
%from HBL batteries india
config.Batt_sizes = [80; 100; 135; 150; 180]; %amp-hours
%back of the envelope for max
config.Batt_max  = roundn((config.Batt_Days*config.hours_per_day*(config.LED_Num*config.LED_Watt*num.max_consumers ...
     + config.Charger_Watt*num.max_consumers) / ...
    (config.MDOD*config.Batt_Eff))/config.BattVolt,1);
 
increment.pv = 150; %pv increment 
increment.batt = 150; %ah increment



 
if num.consumers <= 10
    config.PV_Options = config.PV_modules;
    config.Battery_Options = config.Batt_sizes;
    
elseif num.consumers >=10 && num.consumers <= num.max_consumers
    

    config.PV_Options = max(config.PV_modules):increment.pv:config.PV_max;
    config.Battery_Options = max(config.Batt_sizes):increment.batt:config.Batt_max;
    
   
end


% create matrices / vectors to store calculated variables

reliability_critical_matrix = zeros(length(config.Battery_Options),...
    length(config.PV_Options));

reliability_noncritical_matrix = zeros(length(config.Battery_Options),...
    length(config.PV_Options));

cost_pv_batt_matrix = zeros(length(config.Battery_Options),...
    length(config.PV_Options));

critical_amount_served = zeros(length(config.Battery_Options),...
    length(config.PV_Options));

noncritical_amount_served = zeros(length(config.Battery_Options),...
    length(config.PV_Options));


pv_per_watt = zeros(length(config.PV_Options),1);
batt_per_kwh = zeros(length(config.Battery_Options),1);

pv_total_cost = zeros(length(config.PV_Options),1);
pv_number_final = zeros(length(config.PV_Options),1);
pv_module_final = zeros(length(config.PV_Options),1);


batt_total_cost = zeros(length(config.Battery_Options),1);
batt_number_final = zeros(length(config.Battery_Options),1);
batt_capacity_final = zeros(length(config.Battery_Options),1);


%% Develop Matrix Entries for Cost and Reliability Matrices
%Still have num.consumers qualifier for > or < n=10

if num.consumers <= 10

for i=1:length(config.Battery_Options)
    
    batt = config.Battery_Options(i);
    batt_per_kwh(i) = cost_batt_per_kwh(batt,voltage);
    
    
    for j=1:length(config.PV_Options)
        
        pv = config.PV_Options(j);
        batt_size = batt*voltage;
        
        available_solar = singlehome_gen_sizing(config,solar,num,pv);
        
        
        [reliability_critical_matrix(i,j),critical_amount_served(i,j), reliability_noncritical_matrix(i,j), noncritical_amount_served(i,j)] = connect_ulink_customer_gen_sizing_all_load(available_solar,battery,loadprof,num,batt_size);

        pv_per_watt(j) = cost_pv_per_watt(pv);
        
        %total cost of pv/batt combination
        cost_pv_batt_matrix(i,j) = batt_per_kwh(i)*batt_size/1000 + pv_per_watt(j)*pv ;       
        
    end
      
    
end


elseif num.consumers >=10 && num.consumers <= num.max_consumers
    
for i=1:length(config.Battery_Options)
    
    
    batt = config.Battery_Options(i);
    [batt_total_cost(i),batt_number_final(i),batt_capacity_final(i)] = cost_batt_total(batt,voltage,config);
    
    
    for j=1:length(config.PV_Options)
        
        pv = config.PV_Options(j);
        batt_size = batt*voltage;
        
        available_solar = singlehome_gen_sizing(config,solar,num,pv);
        
        [reliability_critical_matrix(i,j),critical_amount_served(i,j)] = connect_ulink_customer_gen_sizing_critical(available_solar,battery,loadprof,num,batt_size);
        [reliability_noncritical_matrix(i,j),noncritical_amount_served(i,j)] = connect_ulink_customer_gen_sizing_noncritical(available_solar,battery,loadprof,num,batt_size);
        [pv_total_cost(j),pv_number_final(j),pv_module_final(j)] = cost_pv_total(pv,config);
        
        %total cost of pv/batt combination
        cost_pv_batt_matrix(i,j) = batt_total_cost(i) + pv_total_cost(j) ;       
        
    end
      
    
end

   

end



%% Check Reliability Thresholds

reliability.critical = .99; %night loads (lights)
reliability.noncritical = .9; %incl. day loads (fan)

reliability.range_critical = 0.1; % - range of thresholds
reliability.range_noncritical = .075;

%find indices of matrix that meet noncritical reliability threshold
[ind.noncritical_i,ind.noncritical_j] = find(reliability_noncritical_matrix > (reliability.noncritical-reliability.range_noncritical)); %& ...
%reliability_noncritical_matrix < (reliability.noncritical+reliability.range_noncritical));

%turn into a matrix
%ind.noncritical = [ind.noncritical_i,ind.noncritical_j];

%get subscripts of every matrix entry
ind.noncritical_sub = sub2ind(size(reliability_noncritical_matrix),ind.noncritical_i,ind.noncritical_j);

if isempty(ind.noncritical_sub)
    error('No PV/Batt Combinations Meet NonCritical Reliability Threshold');
end

%find indices of matrix that meet critical reliability threshold
[ind.critical_i,ind.critical_j] = find(reliability_critical_matrix > (reliability.critical-reliability.range_critical)); %& ...
%reliability_critical_matrix < (reliability.critical+reliability.range_critical));

%ind.critical = [ind.critical_i,ind.critical_j];

%get subscripts of every matrix entry
ind.critical_sub = sub2ind(size(reliability_critical_matrix),ind.critical_i,ind.critical_j);

if isempty(ind.critical_sub)
    error('No PV/Batt Combinations Meet Critical Reliability Threshold');
end


%% Overlap of Subscripts / Indices
%if there are more of one element versus the other
ind.cost_sub = [];

if length(ind.critical_sub) <= length(ind.noncritical_sub)
    ind.cost_sub = find(ismember(ind.critical_sub,ind.noncritical_sub));
    ind.cost_sub_use = ind.critical_sub(ind.cost_sub);
    
else
    ind.cost_sub = find(ismember(ind.noncritical_sub,ind.critical_sub));
    ind.cost_sub_use = ind.noncritical_sub(ind.cost_sub);
    
    
end

%get indices back

%find lowest cost combination
[ind.cost_value, ind.cost_index] = min(cost_pv_batt_matrix(ind.cost_sub_use));

ind.cost_index_use = ind.cost_sub_use(ind.cost_index);

%[ind.cost_i,ind.cost_j] = ind2sub(size(reliability_critical_matrix),ind.cost_sub);
%find combination of battery and pv
[ind.batt_final,ind.pv_final]= ind2sub(size(cost_pv_batt_matrix),ind.cost_index_use);

%total amount of gen asset needed
config.pv_final = config.PV_Options(ind.pv_final);
config.batt_final = config.Battery_Options(ind.batt_final);

%costs
config.pv_final_cost = pv_total_cost(ind.pv_final);
config.batt_final_cost = batt_total_cost(ind.batt_final);

%size of each module/battery
config.pv_module_final = pv_module_final(ind.pv_final);
config.batt_capacity_final = batt_capacity_final(ind.batt_final);

%number of each module / battery
config.pv_number_final = pv_number_final(ind.pv_final);
config.batt_number_final = batt_number_final(ind.batt_final);


critical_final = reliability_critical_matrix(ind.batt_final,ind.pv_final);
noncritical_final = reliability_noncritical_matrix(ind.batt_final,ind.pv_final);

total_amount_served = critical_amount_served(ind.batt_final,ind.pv_final)+noncritical_amount_served(ind.batt_final,ind.pv_final);


%config.pv_per_watt = pv_per_watt(ind.pv_final);
%config.batt_final_per_kwh = batt_per_kwh(ind.batt_final);




%% Cost, Revenue, and DCF Analysis
% note that cost is being updated to include replacement costs
%adding system lifetime

%[cost, revenue] = cost_uLink_csr_gen_reliability(num,demand,config, total_amount_served,voltage)
[cost, revenue] = cost_uLink_csr_gen_reliability_large_n(num,demand,config, total_amount_served,voltage)


%% Display Values / Plot
%disp(strcat({'Average Yearly Reliability serving Total Load is '},{num2str(network.total.reliability(end))}));
disp(strcat({'Final PV Array Size, Watts '},{num2str(config.pv_final)}));
disp(strcat({'Final Batt Bank Size, Ah '},{num2str(config.batt_final)}));
disp(strcat({'Final PV Module W '},{num2str(config.pv_module_final)}));
disp(strcat({'Final PV Module No. '},{num2str(config.pv_number_final)}));
disp(strcat({'Final Batt Cap. Ah'},{num2str(config.batt_capacity_final)}));
disp(strcat({'Final Batt Cap No. '},{num2str(config.batt_number_final)}));


disp(strcat({'Final Cost $ '},{num2str(ind.cost_value)}));
disp(strcat({'Final Critical Reliability % '},{num2str(critical_final)}));
disp(strcat({'Final NonCritical Reliability % '},{num2str(noncritical_final)}));

disp(strcat({'LCOE Over 5 Years is (in $/kWh) '},{num2str(cost.lcoe)}));
disp(strcat({'Cap Cost $/W '},{num2str(cost.initial_per_watt)}));
disp(strcat({'NPV Based on Costs '},{num2str(cost.npv_total)}));
disp(strcat({'Capital Cost '},{num2str(cost.initial)}));
disp(strcat({'Breakeven Monthly (w/ No Fee) in $'},{num2str(revenue.breakeven_monthly)}));
disp(strcat({'Breakeven Monthly (w/ Fee) in $'},{num2str(revenue.breakeven_monthly_connection)}));



config.cost_string = {'Solar','Battery','Wiring','Ethernet',...
    'Lighting','Fans','Poles','Installation'};

config.explode = ones(length(config.cost_string),1);
figure(15);
pie([cost.solar_panel,cost.battery,cost.wiring,cost.ethernet,...
    cost.lighting,cost.fan,cost.poles,cost.installation],config.explode,config.cost_string);
%title('Distribution of Capital Cost w/o uLink Boxes');


figure(16);
pie([cost.initial,-cost.maintenance_total+cost.replacement_npv])
title('Capital and Operating Cost (Yearly) Breakdown in %');

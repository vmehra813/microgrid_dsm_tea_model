function [loadprof, load] = CreateDemandProfiles_1016_VM(num,demand, solar)

%% INPUT FILES


%% activity separation -- light load, cell phone load, fan load
%representing -- realistic operation // cant serve proportional load --
% change time span -- half hour


%% INITIALIZE LOAD PROFILE MATRICES for RESIDENTIAL


loadprof.critical=[];
loadprof.non_critical=[];
loadprof.total = [];

%NoTVNight = xlsread(config.demandInputs,'Residential','NoTVNight');

load.ProbOwn = xlsread(demand.inputs, demand.type.residential, 'ProbOwn');
load.NoTVNight = xlsread(demand.inputs, demand.type.residential, 'NoTVNight');
load.NoTVDay = xlsread(demand.inputs, demand.type.residential, 'NoTVDay');
load.NoLight3Hours = xlsread(demand.inputs, demand.type.residential, 'NoLight3Hours');
load.NoLight2Hours = xlsread(demand.inputs, demand.type.residential, 'NoLight2Hours');
load.NoLight1Hours = xlsread(demand.inputs, demand.type.residential, 'NoLight1Hours');
load.NoPhoneHours = xlsread(demand.inputs, demand.type.residential, 'NoPhoneHours');

%DemandGrowthRate = xlsread('demandInputs_uLink.xlsx', 'Residential', 'DemandGrowthRate');
%BuildingAmp = demandInputs.BuildingAmp;

load.DayDayVar = xlsread(demand.inputs, demand.type.residential, 'DayDayVar');
load.ApplianceVar = xlsread(demand.inputs, demand.type.residential, 'ApplianceVar');
load.Critical = xlsread(demand.inputs, demand.type.residential, 'Critical');

load.Power = 1000*(xlsread(demand.inputs, demand.type.residential, 'Power')); %Watts vs. kW
load.DailyDuration = xlsread(demand.inputs, demand.type.residential, 'DailyDuration');
load.EnabledWhen = xlsread(demand.inputs, demand.type.residential, 'EnabledWhen');
load.FanThresholdTemp = xlsread(demand.inputs, demand.type.residential, 'FanThresholdTemp');
load.AppNum = xlsread(demand.inputs, demand.type.residential, 'AppNum');

load.hour_year_vector = zeros(24,1);
load.hour_year_vector(1)=0;
load.hour_year_vector(2)=1;
load.hour_year_vector(3)=2;
load.hour_year_vector(4)=3;
load.hour_year_vector(5)=4;
load.hour_year_vector(6)=5;
load.hour_year_vector(7)=6;
load.hour_year_vector(8)=7;
load.hour_year_vector(9)=8;
load.hour_year_vector(10)=9;
load.hour_year_vector(11)=10;
load.hour_year_vector(12)=11;
load.hour_year_vector(13)=12;
load.hour_year_vector(14)=13;
load.hour_year_vector(15)=14;
load.hour_year_vector(16)=15;
load.hour_year_vector(17)=16;
load.hour_year_vector(18)=17;
load.hour_year_vector(19)=18;
load.hour_year_vector(20)=19;
load.hour_year_vector(21)=20;
load.hour_year_vector(22)=21;
load.hour_year_vector(23)=22;
load.hour_year_vector(24)=23;

load.hour_year_vector_final = repmat(load.hour_year_vector,[365,1]);



%load.HourofDay = xlsread(demand.inputs, demand.type.residential, 'HourofDay');

%DurationStdDev = xlsread(config.demandInputs,'Residential_Base','DurationStdDev');
%Appliances = xlsread(config.demandInputs,'Residential_Base','Appliances');

%%% assign variable names to each value (example Light1_Critical,
%%% Light1_Power, Light1DailyDuration

loadprof.light1_criticality = load.Critical(1);
loadprof.light2_criticality = load.Critical(2);
loadprof.light3_criticality = load.Critical(3);
loadprof.fan_criticality = load.Critical(4);
loadprof.tv_night_criticality = load.Critical(5);
loadprof.tv_day_criticality = load.Critical(6);
loadprof.tv_standby_criticality = load.Critical(7);
loadprof.phone_criticality = load.Critical(10);


%average daily duration use of each appliance
light1_dailydur = load.DailyDuration(1);
light2_dailydur = load.DailyDuration(2);
light3_dailydur = load.DailyDuration(3);
%fan_dailydur = DailyDuration(4);
tv_night_dailydur = load.DailyDuration(5);
tv_day_dailydur = load.DailyDuration(6);
phone_dailydur = load.DailyDuration(10);
%tv_standby_dailydur = DailyDuration(7);

%power draw per appliance
loadprof.light1_power = load.Power(1);
loadprof.light2_power = load.Power(2);
loadprof.light3_power = load.Power(3);
loadprof.fan_power = load.Power(4);
loadprof.tv_night_power = load.Power(5);
loadprof.tv_day_power = load.Power(6);
loadprof.tv_standby_power = load.Power(7);
loadprof.phone_power = load.Power(10);

light1_probown = load.ProbOwn(1);
light2_probown = load.ProbOwn(2);
light3_probown = load.ProbOwn(3);
fan_probown = load.ProbOwn(4);
tv_night_probown = load.ProbOwn(5);
tv_day_probown = load.ProbOwn(6);
tv_standby_probown = load.ProbOwn(7);
phone_probown = load.ProbOwn(10);


light1_appnum = load.AppNum(1);
light2_appnum  = load.AppNum(2);
light3_appnum  = load.AppNum(3);
fan_appnum  = load.AppNum(4);
tv_night_appnum  = load.AppNum(5);
tv_day_appnum  = load.AppNum(6);
tv_standby_appnum  = load.AppNum(7);
phone_appnum = load.AppNum(10);


light1_enabledwhen = load.EnabledWhen(1);
light2_enabledwhen = load.EnabledWhen(2);
light3_enabledwhen = load.EnabledWhen(3);
%fan_enabledwhen = load.EnabledWhen(4);
fan_availhours_threshold_temp = load.EnabledWhen(4);
fan_nhours_threshold_temp = load.FanThresholdTemp;
%hourofday = load.HourofDay;



%% CREATE DEMAND PROFILES FOR RESIDENTIAL CUSTOMERS

%-1 to 1 is uniformly randomly distributed number from -1 to 1
%adding variability to load profile


for iconsumer = 1:num.consumers
    daily_variations = load.DayDayVar(1)*(2*rand(365,1)-1);
    light1_appvar = load.ApplianceVar(1)*(2*rand(365,1)-1);
    light2_appvar = load.ApplianceVar(2)*(2*rand(365,1)-1);
    light3_appvar = load.ApplianceVar(3)*(2*rand(365,1)-1);
    fan_appvar = load.ApplianceVar(4)*(2*rand(365,1)-1);
    tv_night_appvar = load.ApplianceVar(5)*(2*rand(365,1)-1);
    tv_day_appvar = load.ApplianceVar(6)*(2*rand(365,1)-1);
    tv_standby_appvar = load.ApplianceVar(7)*(2*rand(365,1)-1);
    phone_appvar = load.ApplianceVar(10)*(2*rand(365,1)-1);
    
    light1_daydayvar = daily_variations;
    light2_daydayvar = daily_variations;
    light3_daydayvar = daily_variations;
    fan_daydayvar = daily_variations;
    tv_night_daydayvar = daily_variations;
    tv_day_daydayvar = daily_variations;
    tv_standby_daydayvar = daily_variations;
    phone_daydayvar = daily_variations;
    
    loadprof.fandemand = [];
    loadprof.lightingdemand1 = [];
    loadprof.lightingdemand2 = [];
    loadprof.lightingdemand3 = [];
    loadprof.tv_night_demand = [];
    loadprof.tv_day_demand = [];
    loadprof.tv_standby_demand = [];
    loadprof.phone_demand = [];
    
    
    %each appliance power draw
    
    
    
    for day = 1:365
        hours = (((day-1)*24) + 1):(((day-1)*24) + 24);
        
        %lighting 1
        lighting_hours_avail_logical = solar.irradiance_plane(hours) < light1_enabledwhen;
        
        %where is no light hour defined initially?
        %what is for loop doing
        for nolighthour = load.NoLight1Hours'
            lighting_hours_avail_logical(load.hour_year_vector_final(hours)==nolighthour) = false;
        end
        lighting_hours_avail = find(lighting_hours_avail_logical);
        lighting_nhours1 = light1_dailydur.*(1+light1_daydayvar(day)+light1_appvar(day));
        lighting_chanceon_perhour1 = lighting_nhours1./length(lighting_hours_avail);
        lighting_hours1 = lighting_hours_avail(rand(size(lighting_hours_avail))<lighting_chanceon_perhour1);
        lighting_demand_today1 = zeros(24,1);
        lighting_demand_today1(lighting_hours1) = loadprof.light1_power*light1_appnum*light1_probown;
        loadprof.lightingdemand1 = [loadprof.lightingdemand1; lighting_demand_today1];
        
        %lighting 2
        lighting_hours_avail_logical = solar.irradiance_plane(hours) < light2_enabledwhen;
        for nolighthour = load.NoLight2Hours'
            lighting_hours_avail_logical(load.hour_year_vector_final(hours)==nolighthour) = false;
        end
        lighting_hours_avail = find(lighting_hours_avail_logical);
        lighting_nhours2 = light2_dailydur.*(1+light2_daydayvar(day)+light2_appvar(day));
        lighting_chanceon_perhour2 = lighting_nhours2./length(lighting_hours_avail);
        lighting_hours2 = lighting_hours_avail(rand(size(lighting_hours_avail))<lighting_chanceon_perhour2);
        lighting_demand_today2 = zeros(24,1);
        lighting_demand_today2(lighting_hours2) = loadprof.light2_power*light2_appnum*light2_probown;
        loadprof.lightingdemand2 = [loadprof.lightingdemand2; lighting_demand_today2];
        
        %lighting 3
        lighting_hours_avail_logical = solar.irradiance_plane(hours) < light3_enabledwhen;
        for nolighthour = load.NoLight3Hours'
            lighting_hours_avail_logical(load.hour_year_vector_final(hours)==nolighthour) = false;
        end
        lighting_hours_avail = find(lighting_hours_avail_logical);
        lighting_nhours3 = light3_dailydur.*(1+light3_daydayvar(day)+light3_appvar(day));
        lighting_chanceon_perhour3 = lighting_nhours3./length(lighting_hours_avail);
        lighting_hours3 = lighting_hours_avail(rand(size(lighting_hours_avail))<lighting_chanceon_perhour3);
        lighting_demand_today3 = zeros(24,1);
        lighting_demand_today3(lighting_hours3) = loadprof.light3_power*light3_appnum*light3_probown;
        loadprof.lightingdemand3 = [loadprof.lightingdemand3; lighting_demand_today3];
        
        %fan -- updated
        fan_hours_avail = find(solar.tAmbient(hours) > fan_availhours_threshold_temp);
        fan_dailydur = sum(solar.tAmbient(hours) > fan_nhours_threshold_temp);
        fan_nhours = fan_dailydur.*(1+fan_daydayvar(day)+fan_appvar(day));
        fan_chanceon_perhour = fan_nhours/length(fan_hours_avail);
        fan_hours = fan_hours_avail(rand(size(fan_hours_avail))<fan_chanceon_perhour);
        fan_demand_today = zeros(24,1);
        fan_demand_today(fan_hours) = loadprof.fan_power*fan_appnum*fan_probown;
        loadprof.fandemand = [loadprof.fandemand; fan_demand_today];
        
        
        %         %tv night - updated
        %         tv_hours_avail_logical = true(24,1);
        %         for notvhour = NoTVNight' %hours when TVnight is off
        %             %NOTE: geodata.hourofday starts at 0 for hour 1
        %             tv_hours_avail_logical(geodata.hourofday(hours)==notvhour) = false;
        %         end
        %         tv_hours_avail = find(tv_hours_avail_logical);
        %         tv_nhours = tv_night_dailydur.*(1+tv_night_appvar(day)+tv_night_daydayvar(day));
        %         tv_chanceon_perhour = tv_nhours./length(tv_hours_avail);
        %         tv_hours = tv_hours_avail(rand(size(tv_hours_avail))<tv_chanceon_perhour);
        %         tv_night_demand_today = zeros(24,1); %*tv_standbypower*tv_number*tv_probown;
        %         tv_night_demand_today(tv_hours) = loadprof.tv_night_power*tv_night_appnum*tv_night_probown;
        %         tv_night_demand = [tv_night_demand; tv_night_demand_today];
        %
        %         %tv day - updated
        %         tv_hours_avail_logical = true(24,1);
        %
        %         for notvhour = NoTVDay' %hours when TVDay is off
        %             tv_hours_avail_logical(geodata.hourofday(hours)==notvhour) = false;
        %         end
        %         tv_hours_avail = find(tv_hours_avail_logical);
        %         tv_nhours = tv_day_dailydur.*(1+tv_day_appvar(day)+tv_day_daydayvar(day));
        %         tv_chanceon_perhour = tv_nhours./length(tv_hours_avail);
        %         tv_hours = tv_hours_avail(rand(size(tv_hours_avail))<tv_chanceon_perhour);
        %         tv_day_demand_today = zeros(24,1);%tv_demand_today = ones(24,1)*tv_standbypower*tv_number*tv_probown;
        %         tv_day_demand_today(tv_hours) = loadprof.tv_day_power*tv_day_appnum*tv_day_probown;
        %         tv_day_demand = [tv_day_demand; tv_day_demand_today];
        %
        %         %tv standby - updated
        %         tv_day_and_night_today = tv_night_demand_today + tv_day_demand_today;
        %         tv_standby_demand_today = zeros(24,1);
        %         tv_standby_demand_today(tv_day_and_night_today==0)=loadprof.tv_standby_power*tv_standby_appnum*tv_standby_probown;
        %         tv_standby_demand = [tv_standby_demand; tv_standby_demand_today];
        
        
        %%
        %cell phone charger
        phone_hours_avail_logical = true(24,1);
        
        for nophonehours = load.NoPhoneHours'
            phone_hours_avail_logical(hours==nophonehours) = false;
        end
        
        
        phone_hours_avail = find(phone_hours_avail_logical);
        phone_nhours = phone_dailydur.*(1+phone_appvar(day)+phone_daydayvar(day));
        phone_chanceon_perhour = phone_nhours./length(phone_hours_avail);
        phone_hours = phone_hours_avail(rand(size(phone_hours_avail))<phone_chanceon_perhour);
        phone_demand_today = zeros(24,1);
        phone_demand_today(phone_hours) = loadprof.phone_power*phone_appnum*phone_probown;
        loadprof.phone_demand = [loadprof.phone_demand; phone_demand_today];
        
    end
    
    %critical and non critical load separation
    
    %per appliance
    
    %multiply by 1000 to get in W -- DONE BEFORE NOW.
    loadprof.non_critical(:, iconsumer)=(loadprof.fandemand*~loadprof.fan_criticality+loadprof.lightingdemand1*~loadprof.light1_criticality+loadprof.lightingdemand2*~loadprof.light2_criticality+loadprof.lightingdemand3*~loadprof.light3_criticality+loadprof.phone_demand*~loadprof.phone_criticality);
    
    loadprof.non_critical_fan(:,iconsumer) = (loadprof.fandemand*~loadprof.fan_criticality);
    loadprof.non_critical_lighting(:,iconsumer) = (loadprof.lightingdemand1*~loadprof.light1_criticality+loadprof.lightingdemand2*~loadprof.light2_criticality+loadprof.lightingdemand3*~loadprof.light3_criticality);
    loadprof.non_critical_phone(:,iconsumer) = (loadprof.phone_demand*~loadprof.phone_criticality);
    
    loadprof.critical(:, iconsumer)=(loadprof.fandemand*loadprof.fan_criticality+loadprof.lightingdemand1*loadprof.light1_criticality+loadprof.lightingdemand2*loadprof.light2_criticality+loadprof.lightingdemand3*loadprof.light3_criticality+loadprof.phone_demand*loadprof.phone_criticality);
    
    
    %per appliance
    loadprof.critical_fan(:,iconsumer) = (loadprof.fandemand*loadprof.fan_criticality);
    loadprof.critical_lighting(:,iconsumer) = (loadprof.lightingdemand1*loadprof.light1_criticality+loadprof.lightingdemand2*loadprof.light2_criticality+loadprof.lightingdemand3*loadprof.light3_criticality);
    loadprof.critical_phone(:,iconsumer) = (loadprof.phone_demand*loadprof.phone_criticality);
    
    
    
end


loadprof.total = loadprof.critical + loadprof.non_critical;

% first a columns are a profiles

loadprof.a.total = loadprof.total(:,1:num.A);
loadprof.a.critical = loadprof.critical(:,1:num.A);

% thereafter, are b profiles

loadprof.b.total = loadprof.total(:,1:num.consumers);
loadprof.b.critical = loadprof.critical(:,1:num.consumers);


%% Aggregate B Loads // Add Up Recursive Individual Appliance Demands
% separate all B loads from matrix


loadprof.cumulative.total = zeros(length(loadprof.total),num.consumers);
loadprof.cumulative.critical = zeros(length(loadprof.critical),num.consumers);


loadprof.cumulative.non_critical_fan = zeros(length(loadprof.total),num.consumers);
loadprof.cumulative.critical_fan = zeros(length(loadprof.critical),num.consumers);

loadprof.cumulative.non_critical_light = zeros(length(loadprof.total),num.consumers);
loadprof.cumulative.critical_light = zeros(length(loadprof.critical),num.consumers);

loadprof.cumulative.non_critical_phone = zeros(length(loadprof.total),num.consumers);
loadprof.cumulative.critical_phone = zeros(length(loadprof.critical),num.consumers);


loadprof.cumulative.total(:,1) = loadprof.b.total(:,1);
loadprof.cumulative.critical(:,1) = loadprof.b.critical(:,1);

loadprof.cumulative.non_critical_fan = loadprof.non_critical_fan(:,1);
loadprof.cumulative.critical_fan = loadprof.critical_fan(:,1);

loadprof.cumulative.non_critical_light = loadprof.non_critical_lighting(:,1);
loadprof.cumulative.critical_light = loadprof.critical_lighting(:,1);

loadprof.cumulative.non_critical_phone = loadprof.non_critical_phone(:,1);
loadprof.cumulative.critical_phone = loadprof.critical_phone(:,1);



% create cumulative vectors recursively, including at appliance level. 
for k=2:num.consumers
    
    loadprof.cumulative.total(:,k) = loadprof.total(:,k) + loadprof.cumulative.total(:,k-1);
    loadprof.cumulative.critical(:,k) = loadprof.critical(:,k) + loadprof.cumulative.critical(:,k-1);
    
    loadprof.cumulative.non_critical_fan(:,k) = loadprof.non_critical_fan(:,k) + loadprof.cumulative.non_critical_fan(:,k-1);
    loadprof.cumulative.critical_fan(:,k) = loadprof.critical_fan(:,k) + loadprof.cumulative.critical_fan(:,k-1);
    
    loadprof.cumulative.non_critical_light(:,k) = loadprof.non_critical_lighting(:,k) + loadprof.cumulative.non_critical_light(:,k-1);
    loadprof.cumulative.critical_light(:,k) = loadprof.critical_lighting(:,k) + loadprof.cumulative.critical_light(:,k-1);
    
    loadprof.cumulative.non_critical_phone(:,k) = loadprof.non_critical_phone(:,k) + loadprof.cumulative.non_critical_phone(:,k-1);
    loadprof.cumulative.critical_phone(:,k) = loadprof.critical_phone(:,k) + loadprof.cumulative.critical_phone(:,k-1);
    
     
    
end




%% Including Distrubtion Losses


loadprof.b_critical_current = zeros(length(loadprof.total),num.consumers);
loadprof.b_total_current = zeros(length(loadprof.critical),num.consumers);



loadprof.b_total_current = loadprof.b.total / demand.distr_voltage;
loadprof.b_critical_current = loadprof.b.critical / demand.distr_voltage;

%i^2*r
loadprof.b_total_current_losses = bsxfun(@times,loadprof.b_total_current.^2,demand.distances_b.');
loadprof.b_critical_current_losses = bsxfun(@times,loadprof.b_critical_current.^2,demand.distances_b.'); % I^2 R

loadprof.b_cumulative_total_current_losses = zeros(length(loadprof.b_total_current_losses),num.consumers);
loadprof.b_cumulative_critical_current_losses = zeros(length(loadprof.b_critical_current_losses),num.consumers);


loadprof.b_cumulative_total_current_losses(:,1) = loadprof.b_total_current_losses(:,1);
loadprof.b_cumulative_critical_current_losses(:,1)= loadprof.b_critical_current_losses(:,1);

for l=2:num.consumers
    loadprof.b_cumulative_total_current_losses(:,l) = loadprof.b_total_current_losses(:,l)+loadprof.b_cumulative_total_current_losses(:,l-1);
    loadprof.b_cumulative_critical_current_losses(:,l)= loadprof.b_critical_current_losses(:,l) + loadprof.b_cumulative_critical_current_losses(:,l-1);
    
end



%% Include DC / DC Converter Losses

%Load converter: each household

demand.load_fcn_load_critical = polyval(demand.load_fcn,loadprof.critical);
demand.load_fcn_load_total = polyval(demand.load_fcn,loadprof.total);

for i=1:length(demand.load_fcn_load_critical)
    
    for j=1:num.consumers
        
        if demand.load_fcn_load_critical(i,j) <.67
            demand.load_fcn_load_critical(i,j) = .67;
        elseif demand.load_fcn_load_critical(i,j) > 0.95
            demand.load_fcn_load_critical(i,j) = 0.95;
            
        end
        
        if demand.load_fcn_load_total(i,j) < .67
            demand.load_fcn_load_total(i,j) = .67;
        elseif demand.load_fcn_load_total(i,j) > .95
            demand.load_fcn_load_total(i,j) = 0.95;
        end
        
    end
end

demand.load_conv_loss_critical = (loadprof.critical).*(1-demand.load_fcn_load_critical);
demand.load_conv_loss_total = (loadprof.total).*(1-demand.load_fcn_load_total);



loadprof.b_cumulative_total_load_losses = zeros(length(loadprof.b_total_current_losses),num.consumers);
loadprof.b_cumulative_critical_load_losses = zeros(length(loadprof.b_critical_current_losses),num.consumers);

loadprof.b_cumulative_total_load_losses(:,1) = demand.load_conv_loss_total(:,1);
loadprof.b_cumulative_critical_load_losses(:,1) = demand.load_conv_loss_critical(:,1);

for l=2:num.consumers
    loadprof.b_cumulative_total_load_losses(:,l) = demand.load_conv_loss_total(:,l)+loadprof.b_cumulative_total_load_losses(:,l-1);
    loadprof.b_cumulative_critical_load_losses(:,l)= demand.load_conv_loss_critical(:,l) + loadprof.b_cumulative_critical_load_losses(:,l-1);
    
end


%Link converter: entire load profile
if num.consumers <=10
    
    demand.link_fcn_load_critical = polyval(demand.link_fcn,loadprof.cumulative.critical);
    demand.link_fcn_load_total = polyval(demand.link_fcn,loadprof.cumulative.total);
    
    for i=1:length(demand.load_fcn_load_critical)
        
        for j=1:num.consumers
            
            if demand.link_fcn_load_critical(i,j) <.58
                demand.link_fcn_load_critical(i,j) = .58;
            elseif demand.link_fcn_load_critical(i,j) > 0.95
                demand.link_fcn_load_critical(i,j) = 0.95;
                
            end
            
            if demand.link_fcn_load_total(i,j) < .58
                demand.link_fcn_load_total(i,j) = .58;
            elseif demand.link_fcn_load_total(i,j) > 0.95
                demand.link_fcn_load_total(i,j) = 0.95;
            end
            
        end
    end
    
    
    demand.link_conv_loss_critical =(loadprof.cumulative.critical).*(1-demand.link_fcn_load_critical);
    demand.link_conv_loss_total =(loadprof.cumulative.total).*(1-demand.link_fcn_load_total);
    
    
    %% Including 3 Losses: Wiring, Load, Link Converters
    
    loadprof.cumulative.total = loadprof.cumulative.total + loadprof.b_cumulative_total_current_losses + loadprof.b_cumulative_total_load_losses  + demand.link_conv_loss_total; %included losses in demand as energy that needs to be served
    loadprof.cumulative.critical = loadprof.cumulative.critical + loadprof.b_cumulative_critical_current_losses + loadprof.b_cumulative_critical_load_losses  + demand.link_conv_loss_critical;
    
else
    
    loadprof.cumulative.total = 1.05*(loadprof.cumulative.total + loadprof.b_cumulative_total_current_losses + loadprof.b_cumulative_total_load_losses); %included losses in demand as energy that needs to be served
    loadprof.cumulative.critical = 1.05*(loadprof.cumulative.critical + loadprof.b_cumulative_critical_current_losses + loadprof.b_cumulative_critical_load_losses);
    
end




%% Final Vectors to Use
%concatenate two matrices
loadprof.final.total = loadprof.cumulative.total;

loadprof.final.critical = loadprof.cumulative.critical;

loadprof.final.non_critical = loadprof.final.total - loadprof.final.critical;


%% Plot Load Profiles
% figure;
% subplot(2,1,1), plot(loadprof.final.total(:,end),'c');
% xlabel('Hour of Year');
% ylabel('Watts');
% set(gca,'fontsize',15);
% legend('Total Load (5 Houses)');
% subplot(2,1,2), plot(loadprof.final.critical(:,end),'r');
% xlabel('Hour of Year');
% ylabel('Watts');
% set(gca,'fontsize',15);
% legend('Critical Load (5 Houses)');
% 
% fan_tot = loadprof.cumulative.non_critical_fan(:,end)+loadprof.cumulative.critical_fan(:,end);
% light_tot = loadprof.cumulative.non_critical_light(:,end)+loadprof.cumulative.critical_light(:,end);
% 
% figure;
% subplot(2,1,1), plot(1000:4500,fan_tot(1000:4500));
% hold on;
% plot(1000:4500,solar.tAmbient(1000:4500));
% xlabel('Hour of Year');
% %ylabel('Temperature (Celsius) / Watts');
% ylim([0,75]);
% set(gca,'fontsize',15);
% legend('Fan Load (Watts)','Temperature (C)');
% subplot(2,1,2), plotyy(2000:2180,light_tot(2000:2180),2000:2180,solar.irradiance_plane(2000:2180));
% xlabel('Hour of Year');
% %ylabel('Temperature (Celsius) / Watts');
% set(gca,'fontsize',15);
% legend('Light Load (Watts)','Solar Irradiance (Watts/m^2)');
% ylim([0,45]);

%% Plotting Converter Losses 

% figure;
% subplot(2,1,1), plot(demand.link_power,polyval(demand.link_fcn,demand.link_power));
% %title('Fitted Efficiency of Link Converter');
% xlabel('Power Level, W');
% ylabel('Network Boost Converter Efficiency, %');
% set(gca,'fontsize',16)
% 
% subplot(2,1,2), plot(demand.link_conv_loss_total(:,end));
% hold on;
% plot(demand.link_conv_loss_critical(:,end));
% hold off;
% %title('Link Converter Losses, All Houses');
% xlabel('Hour of Year');
% ylabel('Losses, W');
% legend('Network Total Losses', 'Network Critical Losses');
% set(gca,'fontsize',16)
% 
% 
% figure;
% subplot(2,1,1), plot(demand.load_power, polyval(demand.load_fcn,demand.load_power));
% %title('Fitted Efficiency of Load Converter');
% xlabel('Power Level, W');
% ylabel('Household Buck Converter Efficiency, %');
% set(gca,'fontsize',16)
% % 
% subplot(2,1,2), plot(loadprof.b.cumulative.total.load_losses(:,1));
% hold on;
% plot(loadprof.b.cumulative.critical.load_losses(:,1));
% hold off;
% % title('Load Converter Losses, One House');
% xlabel('Hour of Year');
% ylabel('Losses, W');
% legend('Household Total Losses', 'Household Critical Losses');
% set(gca,'fontsize',16)

%% Plotting Wiring Losses
% 
% figure;
% plot(loadprof.b.cumulative.total.current_losses(:,end));
% hold on;
% plot(loadprof.b.cumulative.critical.current_losses(:,end));
% hold off;
% xlabel('Hour of Year');
% ylabel('Losses (Watts)');
% legend('Wiring Total Losses', 'Wiring Critical Losses');
% set(gca,'fontsize',15);





end


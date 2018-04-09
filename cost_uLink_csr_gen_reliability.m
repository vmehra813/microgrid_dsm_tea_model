function [cost, revenue] = cost_uLink_csr_gen_reliability(num,demand,config, total_amount_served,voltage)

%% Capital costs // based on equipment size and number of users

% cost.solar_panel = num.A*config.PV_Size*config.PV_per_watt;
% cost.battery = num.A*config.BattVolt*(1/1000)*config.Batt_per_kwh*config.Batt_Size;
cost.solar_panel = config.pv_final*config.pv_per_watt;
cost.battery = config.batt_final*voltage*config.batt_final_per_kwh/1000;


cost.wiring = sum(config.WireCost*demand.max_distance);
cost.ethernet = sum(config.EthernetCost*demand.max_distance);
cost.lighting = num.consumers*config.LED_Cost*config.LED_Num;  %2 lights per house
cost.fan = num.consumers*config.Fan_Num*config.Fan_Cost; % 1 fan per house
cost.poles = num.consumers*config.Pole_Cost;
cost.installation = config.Install_Cost*num.consumers; %installation
cost.ABox = config.ABox_Cost*num.A;
cost.BBox = config.BBox_Cost*(num.consumers-num.A);


cost.initial = cost.solar_panel + cost.battery + cost.wiring + cost.ethernet ...
    + cost.lighting + cost.fan + cost.poles + cost.installation;
%not including loads

%cost.initial = cost.solar_panel + cost.battery + cost.wiring + cost.ethernet ...
%+ cost.poles + cost.installation;

% including loads
cost.initial_ulink = cost.solar_panel + cost.battery + cost.wiring + cost.ethernet ...
    + cost.lighting + cost.fan + cost.poles + cost.installation + ...
    cost.ABox + cost.BBox;

% cost.initial_ulink = cost.solar_panel + cost.battery + cost.wiring + cost.ethernet ...
%      + cost.poles + cost.installation + ...
%     cost.ABox + cost.BBox;

cost.initial_per_watt = cost.initial / (config.pv_final);

cost.initial_per_watt_ulink = cost.initial_ulink / (config.pv_final);

%% Equipment Replacement Cost: system lifetime is used here
% batt replacement now based on rainflow cycle counting method (daily)
cost.poles_replacement = zeros(config.lifetime,1);
cost.pv_replacement = zeros(config.lifetime,1);
cost.load_replacement = zeros(config.lifetime,1);
cost.box_replacement = zeros(config.lifetime,1);
cost.batt_replacement = zeros(config.lifetime,1);

cost.batt_replace_num = round(config.Battery_Cycles / config.cycles_final_total);
% battery replacement depends on cycles 

j=1;

while j<= config.lifetime
    
    %     if mod(j,config.Battery_Lifetime)==0
    %         cost.batt_replacement(j)= cost.battery;
    %     end
    
    if mod(j,cost.batt_replace_num) == 0
        cost.batt_replacement(j) = cost.battery;
    end
    
    if mod(j,config.Poles_Lifetime)==0
        cost.poles_replacement(j)=cost.poles;
    end
    
    if mod(j,config.PV_Lifetime) ==0
        cost.pv_replacement(j) = cost.solar_panel;
    end
    
    if mod(j,config.Load_Lifetime) == 0
        cost.load_replacement(j) = cost.fan+cost.lighting;
    end
    
    if mod(j,config.Box_Lifetime)==0
        cost.box_replacement(j) = cost.ABox+cost.BBox;
    end
    
    j=j+1;
end


%w/ loads
cost.total_replacement = cost.batt_replacement+cost.poles_replacement + ...
    cost.pv_replacement + cost.load_replacement + cost.box_replacement;

%w/o loads
%cost.total_replacement = cost.batt_replacement+cost.poles_replacement + ...
%    cost.pv_replacement + cost.box_replacement;

cost.replacement_npv = pvvar(cost.total_replacement,config.discount_rate);



%% Maintenance cost, monthly

% THIS IS USED IN GRAPH
cost.maintenance = config.maint_rate * cost.initial_ulink;
% discount over lifetime of system

cost.maintenance_npv = repmat(-cost.maintenance,config.lifetime,1);

cost.npv = vertcat(-cost.initial_ulink,cost.maintenance_npv-cost.total_replacement);

cost.maintenance_total = pvvar(cost.maintenance_npv,config.discount_rate);

cost.maint_replace_total = pvvar(cost.npv,config.discount_rate);
% over 5 years

% LCOE -- cost above  / lifetime of years of what you generated



%% old code

%cost.replacement = num.A*(cost.solar_panel/cost.solar_lifetime + ...
%     cost.battery/cost.battery_lifetime + (cost.lighting+cost.fan)/cost.load_lifetime ...
%     + cost.ABox / cost.ABox_lifetime + cost.BBox / cost.BBox_lifetime);


%cost.maintenance = num.A*200/60; % 200 Rs / month per installation

%LCOE = cost per month ($) /  consumption (kWh)

% pvvar(cashflow matrix, daily rate, cfdates)
%irr(cash stream)
%india discount rate and interest rate

%distinction between amount served and energy produced

%if using more electricity from same source, LCOE goes down for N
%increasing

% cost.LCOE_fiveyear = [];
% cost.LCOE_monthly = [];
%
% %discount rate
% cost.discount_rate = 0.07;
% cost.inflation = 0.015;


%
% for i = 1:num.consumers
% cost.LCOE_fiveyear(i) = (cost.initial + cost.replacement*num.years*12 + cost.maintenance*12*num.years) ...
%     / (sum(network.total.amountserved.cumulative(:,i))*num.years/1000);
%
% cost.LCOE_monthly(i) = (cost.replacement + cost.maintenance) ...
%     / (sum(network.total.amountserved.cumulative(:,i))/(12*1000));
%
% end



%% Basic revenue model: payback period is used here

revenue.yearly = config.payment*12*num.consumers; % 2$ month payment per b customer (5) if entrepreneur owns asset

revenue.connection = config.connection*num.consumers; % 3500 Rs / connection charge by JUSCO
%revenue.connection = config.connection*cost.initial_ulink; % 20% capital cost 

%revenue.vector_npv = repmat(revenue.yearly,config.lifetime,1);
%revenue.npv = vertcat(revenue.connection,revenue.vector_npv);
%revenue.npv_cost = pvvar(revenue.vector_npv+cost.maintenance_npv,config.discount_rate);


revenue.npv_connection = vertcat(-cost.initial_ulink+revenue.connection,cost.maintenance_npv-cost.total_replacement);

cost.npv_total_connection = pvvar(revenue.npv_connection,config.discount_rate);
cost.npv_total = pvvar(cost.npv,config.discount_rate);


revenue.breakeven_monthly = -payper(config.discount_rate_month,config.payback*12,cost.npv_total)/(num.consumers);

revenue.breakeven_monthly_connection = -payper(config.discount_rate_month,config.payback*12,cost.npv_total_connection)/(num.consumers);

%revenue.vector = vertcat(revenue.connection,repmat(revenue.breakeven_monthly*12*num.consumers,config.lifetime,1));
revenue.vector = vertcat(repmat(revenue.breakeven_monthly*12*num.consumers,config.lifetime,1));
revenue.vector_connection = vertcat(repmat(revenue.breakeven_monthly_connection*12*num.consumers,config.lifetime,1));

cost.lcoe = 1000*(-cost.npv_total / (config.lifetime*total_amount_served));

revenue.irr = irr([cost.npv_total; revenue.vector]);

revenue.irr_connection = irr([cost.npv_total_connection; revenue.vector_connection]);


end



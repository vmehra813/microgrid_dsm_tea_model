function [ batt_per_kwh ] = cost_batt_per_kwh( batt,voltage )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


% batt_per_wh = 1.3073*(batt*voltage)^(-.307);
% batt_per_kwh = batt_per_wh *1000;

%from jamshedpur distributor data

%batt_per_kwh = 288.8*exp(-7.235*batt*voltage/1000) ...
%                + 125.3*exp(-0.09339*batt*voltage/1000);

            
batt_per_kwh = 288.8*exp(-.08682*batt) ...
                + 125.3*exp(-0.001121*batt);


end


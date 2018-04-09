function [ pv_per_watt ] = cost_pv_per_watt( pv )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%pv_per_watt = 34.882*pv^(-.592);
%divide by 2? tarahuja is 0.67 / watt for 250 watt panel

%pv_per_watt = 34.882*pv^(-.592) / 2;

%based on jharkhand distributor prices 
pv_per_watt = 0.5689*exp(-0.04871*pv) + 0.8098*exp(-0.0001196*pv);




end


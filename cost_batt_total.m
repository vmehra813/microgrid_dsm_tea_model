function [ batt_total_cost,batt_number_final,batt_capacity_final ] = cost_batt_total( batt,voltage,config )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

batt_vector = config.Batt_sizes;
batt_number = zeros(length(batt_vector),1);
batt_costs = zeros(length(batt_vector),1);

for i=1:length(batt_vector)
    batt_number(i) = ceil(batt/batt_vector(i));
end
    
for j=1:length(batt_vector)
    batt_costs(j) = batt_number(j)*voltage*batt_vector(j)*1.3073*(batt_vector(j)*voltage)^(-.307);
    
end

[batt_total_cost, ind_batt] = min(batt_costs);
batt_number_final = batt_number(ind_batt);
batt_capacity_final = batt_vector(ind_batt);



end
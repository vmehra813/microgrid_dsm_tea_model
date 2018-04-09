function [ pv_total_cost, pv_number_final,pv_module_final ] = cost_pv_total( pv,config )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

pv_vector = config.PV_modules;
pv_number = zeros(length(pv_vector),1);
pv_costs = zeros(length(pv_vector),1);


for i=1:length(pv_vector)
    pv_number(i) = ceil(pv/pv_vector(i));
end
    
for j=1:length(pv_vector)
    pv_costs(j) = pv_number(j)*pv_vector(j)*(34.882*pv_vector(j)^(-.592) /2);
    
end

[pv_total_cost, ind_pv]= min(pv_costs);


pv_number_final = pv_number(ind_pv);
pv_module_final = pv_vector(ind_pv);



end
function [reliability_noncritical,noncritical_amount_served] = connect_ulink_customer_gen_sizing_noncritical(available_solar,battery,loadprof,num,batt)

%% Given remaining SOC from A load, ability to serve another customer.

%right now code can take in multiple a generators but is only using 1 a box
%remaining SOC
%incorporate reliability preferences
%% Total Served for A&B, Hour by Hour
%CHANGED THIS TO BE FOR NON-CRITICAL LOADS
%initialize vectors
network.total.reliability = zeros(num.consumers,1);
network.total.amountserved.cumulative = zeros(length(loadprof.final.non_critical),num.consumers);
network.total.nonserved = zeros(length(loadprof.final.non_critical),num.consumers);
%network.total.amountserved.cumulative(1) = 20;
network.total.solar_load_coincidence = zeros(length(loadprof.final.non_critical),num.consumers);


SOC.total = zeros(length(loadprof.final.non_critical), num.consumers);
SOC.total_spillage = zeros(length(loadprof.final.non_critical),num.consumers);
SOC.total_cycles = zeros(num.consumers,1);

battery.initial_capacity = batt*0.5;


for j=1:num.consumers %column is B customer
    
   
    
    for i=1  %special case for i=1 to introduce initial state of charge of battery
        if available_solar(i,1) >= loadprof.final.non_critical(i,j)
            SOC.total(i,j) = battery.initial_capacity + (available_solar(i,1)-loadprof.final.non_critical(i,j))*battery.charge_eff;
            
        elseif available_solar(i,1) < loadprof.final.non_critical(i,j)
            SOC.total(i,j) = battery.initial_capacity - (loadprof.final.non_critical(i,j)-available_solar(i,1))*battery.discharge_eff;
        end
    end
    
    for i=2:length(loadprof.final.non_critical)  %remaining interations look at previous hour
        if available_solar(i,1) >= loadprof.final.non_critical(i,j)
            
            SOC.total(i,j) = SOC.total(i-1,j)+(available_solar(i,1)-loadprof.final.non_critical(i,j))*battery.charge_eff;
            network.total.solar_load_coincidence(i,j) = loadprof.final.non_critical(i,j)/available_solar(i,1); % solar / load coincidence factor
            
            if SOC.total(i,j) < batt 
                SOC.total(i,j) = SOC.total(i,j);
            else
                SOC.total(i,j) = SOC.total(i-1,j);
                SOC.total_spillage(i,j) = available_solar(i,1) - loadprof.final.non_critical(i,j); %spillage
                
            end
            
        elseif available_solar(i,1) < loadprof.final.non_critical(i,j)
            
            SOC.total(i,j) = SOC.total(i-1,j)-(loadprof.final.non_critical(i,j)-available_solar(i,1))*battery.discharge_eff;
            network.total.solar_load_coincidence(i,j) = available_solar(i,1) / loadprof.final.non_critical(i,j); % solar load coincidence
            
            if SOC.total(i,j) > battery.MDOD*batt %MDOD constraint
                SOC.total(i,j)= SOC.total(i,j);
            else SOC.total(i,j) = SOC.total(i-1,j); %nonserved
                network.total.nonserved(i,j) = loadprof.final.non_critical(i,j) - available_solar(i,1);
                
            end
        end
        
            network.total.amountserved.cumulative(i,j) = available_solar(i,1) + SOC.total(i-1,j) - SOC.total(i,j);
        
            if network.total.amountserved.cumulative(i,j) > loadprof.final.non_critical(i,j)
                network.total.amountserved.cumulative(i,j) = loadprof.final.non_critical(i,j);
            end
            
    
    end
    
    network.total.reliability(j) = (sum(network.total.amountserved.cumulative(:,j)))/ (sum(loadprof.final.non_critical(:,j)));
  
    SOC.total_cycles(j) = length(findpeaks(SOC.total(:,j)/batt,'MinPeakHeight',0.7,'MinPeakDistance',1));
    
    
end

reliability_noncritical = network.total.reliability(end);
noncritical_amount_served = sum(network.total.amountserved.cumulative(:,end));

figure(14);
plot(loadprof.final.non_critical(:,end)-network.total.amountserved.cumulative(:,end));
xlabel('Hour of Year');
ylabel('Watts');
title('Amount of Non Served Energy Across Network, NonCritical Load')





function [reliability_critical, critical_amount_served] = connect_ulink_customer_gen_sizing_critical(available_solar,battery,loadprof,num,batt)

%% Given remaining SOC from A load, ability to serve another customer.

%right now code can take in multiple a generators but is only using 1 a box
%remaining SOC
%incorporate reliability preferences


%% Critical Served for B
%
network.critical.reliability = zeros(num.consumers,1);
network.critical.amountserved.cumulative = zeros(length(loadprof.final.critical),num.consumers);
%network.critical.amountserved.cumulative(1) = 20;

network.critical.percentserved.cumulative = zeros(length(loadprof.final.critical),num.consumers);
network.critical.percentserved.cumulative(1) = 1;
network.critical.nonserved = zeros(length(loadprof.final.critical),num.consumers);
network.critical.solar_load_coincidence = zeros(length(loadprof.final.critical),num.consumers);

SOC.critical = zeros(length(loadprof.final.total), num.consumers);
SOC.critical_spillage = zeros(length(loadprof.final.critical),num.consumers);
SOC.critical_cycles = zeros(num.consumers,1);

battery.initial_capacity = batt*0.5;

for j=1:num.consumers
for i=1  %special case for i=1 to introduce initial state of charge of battery
if available_solar(i,1) >= loadprof.final.critical(i,j)
                SOC.critical(i,j) = battery.initial_capacity + (available_solar(i,1)-loadprof.final.critical(i,j))*battery.charge_eff;
                
elseif available_solar(i,1) < loadprof.final.critical(i,j)
            SOC.critical(i,j) = battery.initial_capacity - (loadprof.final.critical(i,j)-available_solar(i,1))*battery.discharge_eff;
end
end

for i=2:length(loadprof.final.critical)  %remaining interations look at previous hour
    if available_solar(i,1) >= loadprof.final.critical(i,j)
        
        SOC.critical(i,j) = SOC.critical(i-1,j)+(available_solar(i,1)-loadprof.final.critical(i,j))*battery.charge_eff;
        network.critical.solar_load_coincidence(i,j) = loadprof.final.critical(i,j)/available_solar(i,1);
        if SOC.critical(i,j) < batt
            SOC.critical(i,j) = SOC.critical(i,j);
        else
            SOC.critical(i,j) = SOC.critical(i-1,j);
            SOC.critical_spillage(i,j) = available_solar(i,1) - loadprof.final.critical(i,j); %spillage

        end
    
    elseif available_solar(i,1) < loadprof.final.critical(i,j)
         
         SOC.critical(i,j) = SOC.critical(i-1,j)-(loadprof.final.critical(i,j)-available_solar(i,1))*battery.discharge_eff;
         network.critical.solar_load_coincidence(i,j) = available_solar(i,1) / loadprof.final.critical(i,j);
         if SOC.critical(i,j) > battery.MDOD*batt
            SOC.critical(i,j)= SOC.critical(i,j);
         else
             SOC.critical(i,j) = SOC.critical(i-1,j);
             network.critical.nonserved(i,j) = loadprof.final.critical(i,j) - available_solar(i,1); %non served
            
         end
         
    end
    
    network.critical.amountserved.cumulative(i,j) = (available_solar(i,1) + SOC.critical(i-1,j) - SOC.critical(i,j));
    
     if network.critical.amountserved.cumulative(i,j) > loadprof.final.critical(i,j)
                network.critical.amountserved.cumulative(i,j) = loadprof.final.critical(i,j);
     end


    
end

     network.critical.reliability(j) = (sum(network.critical.amountserved.cumulative(:,j)))/ (sum(loadprof.final.critical(:,j)));
     SOC.critical_cycles(j) = length(findpeaks(SOC.critical(:,j)/batt,'MinPeakHeight',0.7,'MinPeakDistance',1));

end

reliability_critical = network.critical.reliability(end);
critical_amount_served = sum(network.critical.amountserved.cumulative(:,end));

%% PLOTTING

figure(13);
plot(loadprof.final.critical(:,end)-network.critical.amountserved.cumulative(:,end));
xlabel('Hour of Year');
ylabel('Watts');
title('Amount of Non Served Energy Across Network, Critical Load')


% figure(9);
% subplot(2,2,3,'align'), plotyy(1:num.consumers,network.total.reliability,1:num.consumers,network.critical.reliability);
% xlabel('Number of B Customers');
% ylabel('Avg Yearly Reliability, Percent');
% title(strcat('Load Served with Decreasing Reliability, A+B Customers, for N=',num.graph))
% legend('Total Reliability', 'Critical Reliability');
% 
% subplot(2,2,1,'align'),plot(1:length(SOC.total),SOC.total(:,num.consumers)/battery.size);
% xlabel('Hour');
% ylabel('Remaining SOC');
% title(strcat(' Total SOC Remaining for N=', num.graph))
% 
% subplot(2,2,2,'align'),plot(1:length(SOC.critical),SOC.critical(:,num.consumers)/battery.size);
% xlabel('Hour');
% ylabel('Remaining SOC');
% title(strcat('Critical SOC Remaining for N=',num.graph))
% 
% 
% 
% 
% figure(10);
% subplot(1,2,1,'align')
% plot(1:length(SOC.critical_spillage),SOC.critical_spillage(:,num.consumers), 'm');
% xlabel('Time');
% ylabel('Amt of Critical Spillage in Wh');
% title(strcat('Spilliage for Serving Critical Load for N=',num.graph));
% 
% subplot(1,2,2,'align')
% plot(1:length(SOC.total_spillage),SOC.total_spillage(:,num.consumers), 'c');
% xlabel('Time');
% ylabel('Amt of Total Spillage in Wh');
% title(strcat('Spillage for Serving Total Load for N=',num.graph));
% 
% 
% figure(11);
% subplot(1,2,1,'align')
% 
% plot(1:length(network.total.nonserved),network.total.nonserved(:,num.consumers), 'c');
% xlabel('Time');
% ylabel('Amt of Total NonServed in Wh');
% title(strcat('Total Non Served Energy, Cumulative for N=',num.graph));
% 
% subplot(1,2,2,'align')
% 
% plot(1:length(network.critical.nonserved),network.critical.nonserved(:,num.consumers), 'm');
% xlabel('Time');
% ylabel('Amt of Critical NonServed in Wh');
% title(strcat('Critical Non Served Energy, Cumulative for N=',num.graph));
% 
% figure(12);
% plotyy(1:num.consumers,SOC.total_cycles,1:num.consumers,SOC.critical_cycles);
% legend('Total Load', 'Critical Load');
% title('Number of Cycles when SOC was above 70% serving Critical and Total Load for Increasing N');





end




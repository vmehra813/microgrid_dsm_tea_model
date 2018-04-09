%% PV / Battery Scale / Curve Fitting Script

%exide, sukam, hbl costs
%from jharkhand distributor in jamshedpur

rs_dollars = 66; %conversion rate rs to dollars
batt_volt = 12; %volts

%% Exide 
%C/20 discharge rate, ah

exide.capacity = [40;
                  40;
                  20;
                  20;
                  40;
                  60;
                  75;
                  75;
                  100;
                  105;
                  120;
                  150;
                  150;
                  200];
                  %40;
                  %75;
                  %100;
                  %150;
                  %200;];
              
              %last 5 are gel VRLA batteries
           
% cost is in rupees 
exide.rs_cost = [ 4183;
                  3676;
                  2662;
                  3163;
                  4703;
                  5862;
                  6478;
                  7281;
                  9448;
                  8663;
                  10680;
                  11821;
                  13662;
                  18414];
                  %6557;
                  %11278;
                  %13560;
                  %19918;
                  %22810;];

              %last 5 are gel VRLA batteries
       
%warranty in years 

exide.warranty = [3;
                  3;
                  3;
                  5;
                  5;
                  3;
                  3;
                  5;
                  5;
                  5;
                  5;
                  3;
                  5;
                  5];
                  %5;
                  %5;
                  %5;
                  %5;
                  %5;];
             

exide.data = [exide.capacity, exide.rs_cost/rs_dollars, exide.warranty];


exide.batt_per_ah = (exide.rs_cost/rs_dollars)./exide.capacity; %$/Ah

%% Amara Raja Batteries
%C/10 discharge rate
amaraja.capacity = [100; %solar battery 
                    120;
                    150;
                    100; %tall tubular 
                    150;
                    165;
                    180;];
                
amaraja.rs_cost = [10750; %rupees 
                  12455;
                  13765;
                   9495;
                   12295;
                   13055;
                   14150;];
              

amaraja.batt_data = [amaraja.capacity, amaraja.rs_cost/rs_dollars];


                   
amaraja.batt_per_ah = (amaraja.rs_cost/rs_dollars)./amaraja.capacity; %$/ah

%% Su Kam Batteries

sukam.batt_capacity = [20;
                       40;
                       75;
                       100;
                       120;
                       150;
                       200;];
                   
sukam.batt_cost = [2415;
                  3743;
                  6158;
                  7245;
                  8453;
                  10505;
                  14490;];
                  
sukam.batt_per_ah = (sukam.batt_cost/rs_dollars)./sukam.batt_capacity;



%% Plotting Battery Costs versus Scale ($/kWh vs. Ah)  

% figure(1);
% scatter(exide.capacity,exide.data(:,2));
% hold on;
% scatter(amaraja.capacity, amaraja.batt_data(:,2));
% scatter(sukam.batt_capacity,sukam.batt_cost/rs_dollars);
% title('Exide/Amaraja/SuKam Battery Costs vs. Capacity');
% legend('Exide','Amaraja','SuKam');
% xlabel('Capacity in Ah');
% ylabel('Costs in $');
% set(gca,'fontsize',14);


all_batt = [exide.capacity; amaraja.capacity; sukam.batt_capacity;]; %ah
all_batt_cost = [1000*exide.batt_per_ah./batt_volt; 1000*amaraja.batt_per_ah./batt_volt; 1000*sukam.batt_per_ah./batt_volt;]; % $/kWh

%conduct exponential fit // OLS 
exp_fit = fit(all_batt,all_batt_cost,'exp2');  % ah capacity vs. $/kWh

figure(2);
scatter(exide.capacity,1000*exide.batt_per_ah./batt_volt);
hold on;
scatter(amaraja.capacity,1000*amaraja.batt_per_ah./batt_volt);
scatter(sukam.batt_capacity,1000*sukam.batt_per_ah./batt_volt);
legend('Exide','Amaraja','SuKam','ExpFit');
title('Exide/Amaraja/SuKam $/kWh Costs vs. Capacity');
xlabel('Capacity in kWh');
ylabel('$/KWh');
set(gca,'fontsize',15);

figure(3);
plot(exp_fit,all_batt,all_batt_cost);
xlabel('Capacity (Ah)');
ylabel('Unit Cost ($/kWh)');
%legend('Data from Exide + Amara Raja + SuKam', 'ExpFit');
%title('Exponential Fit for Battery Costs + Size');
set(gca,'fontsize',15);


%% PV Data

amaraja.pv_capacity = [100;
                       120;
                       150;
                       250;];         

amaraja.pv_cost = [5132;
                   6550;
                   7700;
                   12725;];
               
amaraja.pv_per_watt = (amaraja.pv_cost/rs_dollars)./amaraja.pv_capacity;


sukam.pv_capacity = [10;
                     20;
                     40;
                     50;
                     75;
                     80;
                     100;
                     120;
                     125;
                     150;
                     250;
                     250;
                     250;
                     10;
                     20;
                     40;
                     50;
                     75;
                     80;
                     100;
                     125;
                     150;
                     250;];
                 
sukam.pv_cost = [726;
                 1386;
                 2112;
                 2640;
                 3877;
                 4136;
                 5170;
                 6204;
                 6462;
                 7755;
                 12650;
                 12650;
                 12650;
                 770;
                 1474;
                 2332;
                 2915;
                 4208;
                 4488;
                 5610;
                 7012;
                 8415;
                 13750;];
 %sukam warranty = 25 years           
 %polycrystalline and then monocrystalline
 sukam.pv_per_watt = (sukam.pv_cost/rs_dollars)./sukam.pv_capacity;

 
 %% Plotting PV Scale 
 
figure(4);
scatter(amaraja.pv_capacity,amaraja.pv_cost/rs_dollars);
hold on;
scatter(sukam.pv_capacity,sukam.pv_cost/rs_dollars);
%title('Amaraja/SuKam PV Costs vs. Capacity');
legend('Exide','SuKam');
xlabel('Module Capacity (Watts)');
ylabel('Costs in $');
set(gca,'fontsize',15);

all_pv = [amaraja.pv_capacity; sukam.pv_capacity;];
all_pv_cost = [amaraja.pv_per_watt; sukam.pv_per_watt;];

%conduct exponential fit // OLS 
exp_fit_pv = fit(all_pv,all_pv_cost,'exp2');

figure(5);
scatter(amaraja.pv_capacity,amaraja.pv_per_watt);
hold on;
scatter(sukam.pv_capacity,sukam.pv_per_watt);
legend('Amaraja','SuKam','ExpFit');
%title('Amaraja/SuKam $/W Costs vs. Capacity');
xlabel('Module Capacity (Watts)');
ylabel('Unit Cost ($/Watt)');
set(gca,'fontsize',15);

figure(6);
plot(exp_fit_pv,all_pv,all_pv_cost);
xlabel('Module Capacity (Watts)');
ylabel('Unit Cost ($/Watt)');
%legend('Data from Amara Raja + SuKam', 'ExpFit');
%title('Exponential Fit for PV Costs + Size');
set(gca,'fontsize',15);


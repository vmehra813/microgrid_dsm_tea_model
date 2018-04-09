%Main Script for Cost Analysis of a 5 node basic uLink system
clc
clear all 

%Inputs
num.hours=8760; %24 for 1 day; 8760 for 1 year
num.Households=5;
Battery.Capacity=35*12; %Wh 
PVpanel.Rating=200;
Battery.MinDOD=0.05*Battery.Capacity; %Minimum depth of discharge
Battery.Initial=0.5*Battery.Capacity;
Battery.ChargeEff=0.95;
Battery.DischargeEff=0.95;
PVpanel.RatingAtIrradiance=1000;
Wiring.Length=40; 
Wiring.Rho=4e-8;
Wiring.Size=2.5e-6;%m
Network.Voltage=24;
Network.Inefficiency=0.7;


%% Input data 
Generation.data=xlsread('pvwatts_hourly_patamda_jharkhand','A20:K8779'); %Averged 10 years of data from NREL
Generation.Irradiance=Generation.data(:,8);
Generation.Temperature=Generation.data(:,6);
Generation.Hour=Generation.data(:,3);

load.NoLight1Hours = xlsread('demandInputs_uLink.xlsx','Residential', 'NoLight1Hours');
load.NoLight2Hours = xlsread('demandInputs_uLink.xlsx','Residential', 'NoLight2Hours');
load.NoPhoneHours = xlsread('demandInputs_uLink.xlsx','Residential', 'I27:I50');
load.Power = xlsread('demandInputs_uLink.xlsx','Residential', 'Power');

%% Generation data
Generation.Output= PVpanel.Rating*Generation.Irradiance/PVpanel.RatingAtIrradiance; %Output of solar panel with 1000 W/m2 rating

%% Create Load profile
for j=1:num.Households
[loadprof] = CreateDemand_withRandomness(Generation,load);
Load.data(1:num.hours,j)=loadprof(1:num.hours,1);
end
Load.demand_all=sum(Load.data,2);
Load.Current=Load.data/Network.Voltage;

%% Load Converter Loss
Converter.data_eff=xlsread('PowerElectronics_Eff','LoadConverter','A2:B11');
Converter.load_power=Converter.data_eff(:,1);
Converter.load_Eff=Converter.data_eff(:,2);
Function_LoadConverter_Eff=polyfit(Converter.load_power,Converter.load_Eff,4);
Function_LoadConverter_Eff_withLoad= polyval(Function_LoadConverter_Eff,Load.data);
%%To check 
%  figure
%  plot(Converter.load_power,Converter.load_Eff,'o');
%  hold on
%  plot(Converter.load_power,Converter.load_Eff);
%  f1=polyval(Function_LoadConverter_Eff,Converter.load_power);
%  plot(Converter.data_power,f1,'r--')
Loss.Converter_load=(Load.data).*(100-Function_LoadConverter_Eff_withLoad)/100;
Load.data_andConverter=Load.data+Loss.Converter_load;

%% Distribution Loss
Current_data_andCovnerter=Load.data_andConverter/Network.Voltage;
rho_al= 4*10^-8;%Resititive of alum
Wiring.Res=(Wiring.Rho*Wiring.Length)/Wiring.Size;
Loss.Wiring=(Current_data_andCovnerter.^2).*Wiring.Res;
Load.data_andConverter_andWiring=Load.data_andConverter+Loss.Wiring;

%% Link Converter Loss 
Load.data_all=sum(Load.data_andConverter_andWiring,2); %<<<<<<<<<This has to be changed
Converter.data_eff=xlsread('PowerElectronics_Eff','LinkConverter','A2:B15');
Converter.link_power=Converter.data_eff(:,1);
Converter.link_Eff=Converter.data_eff(:,2);
Function_LinkConverter_Eff=polyfit(Converter.link_power,Converter.link_Eff,10);
Function_LinkConverter_Eff_withLoad= polyval(Function_LinkConverter_Eff,Load.data_all);
% %To check 
%  figure
%  plot(Converter.link_power,Converter.link_Eff,'o');
%  hold on
%  plot(Converter.link_power,Converter.link_Eff);
%  f1=polyval(Function_LinkConverter_Eff,Converter.link_power);
%  plot(Converter.link_power,f1,'r--')
%figure
Loss.Converter_link=(Load.data_all).*(100-Function_LinkConverter_Eff_withLoad)/100;
Load.data_all_andLinkConverter=Load.data_all+Loss.Converter_link;

%% Battery loss
Load.data_all_Battery=Load.data_all_andLinkConverter/0.95;

%% Charge Controller Loss 
Converter.data_eff=xlsread('PowerElectronics_Eff','ChargeController','A2:B17');
Converter.charge_power=Converter.data_eff(:,1);
Converter.charge_Eff=Converter.data_eff(:,2);
Function_chargeConverter_Eff=polyfit(Converter.charge_power,Converter.charge_Eff,4);

% Function_chargeConverter_Eff_withLoad= polyval(Function_chargeConverter_Eff,Load.data_all_Battery);
% Loss.Converter_charge=(Load.data_all_Battery).*(100-Function_chargeConverter_Eff_withLoad)/100;
% Load.data_all_charge=Load.data_all_Battery+Loss.Converter_charge;

Function_chargeConverter_Eff_withGeneration= polyval(Function_chargeConverter_Eff,Generation.Output);
Loss.Generation_Converter_charge=(Generation.Output).*(100-Function_chargeConverter_Eff_withGeneration)/100;
Effective.Powertobattery=Generation.Output+Loss.Generation_Converter_charge;
%To check 
%  figure
%  plot(Converter.charge_power,Converter.charge_Eff,'o');
%  hold on
%  plot(Converter.charge_power,Converter.charge_Eff);
%  f1=polyval(Function_chargeConverter_Eff,Generation.Output);
%  plot(Generation.Output,f1,'b--')
%% Battery Energy 
Battery.Energy(1)=Battery.Initial;
Battery.Spillage=zeros(num.hours,1);
Load.notserved=zeros(num.hours,1);
Battery.state_charging(1)=0;
Battery.state_discharging(1)=0;
for i=1:num.hours-1
    if Effective.Powertobattery(i)>=Load.data_all_andLinkConverter(i) %Charging 
        Battery.state_charging(i+1)=1;%charging 
        Battery.Energy(i+1)= (Effective.Powertobattery(i)-Load.data_all_andLinkConverter(i))*Battery.ChargeEff+Battery.Energy(i);
        if Battery.Energy(i+1)>Battery.Capacity
            Battery.Energy(i+1)=Battery.Capacity;
            Battery.Spillage(i+1)=(Effective.Powertobattery(i)-Load.data_all_andLinkConverter(i));
        end
    else %Discharging
        Battery.state_discharging(i+1)=-1;%charging 
        Battery.Energy(i+1)= (Effective.Powertobattery(i)-Load.data_all_andLinkConverter(i))*Battery.DischargeEff+Battery.Energy(i);
        if Battery.Energy(i+1)<Battery.MinDOD
            Battery.Energy(i+1)=Battery.MinDOD;
            Load.notserved(i+1)= (Load.data_all_andLinkConverter(i)-Effective.Powertobattery(i));
          
        end
    end
end

Reliability=(1-sum(Load.notserved)/sum(Load.data_all_andLinkConverter))*100;

%% Cost Analysis 
num.Years=15;
Battery.CostperWh=0.5;% $/Wh
Battery.Cost=Battery.CostperWh*Battery.Capacity;
Battery.Life=5; %years
Battery.Amount=num.Years/Battery.Life; 

PVpanel.CostperW=0.7;
PVpanel.Cost=PVpanel.CostperW*PVpanel.Rating;

Wiring.Costperm=0.06+0.03; %Cost of Aluminum and Cost of Cat3 cable
Wiring.Cost=Wiring.Costperm*Wiring.Length*num.Households;
Wiring.PolesCostperpole=3.33; 
Wiring.PolesCost=Wiring.PolesCostperpole*num.Households;%1 pole per household
Wiring.PolesLife=5; %years
Wiring.PolesAmount=num.Years/Wiring.PolesLife;

Electronics.CostA=40;
Electronics.CostB=20*num.Households;

Cost.Initial=Battery.Cost+ PVpanel.Cost+ Wiring.Cost+ Wiring.PolesCost+ Electronics.CostA+ Electronics.CostB;
Cost.Total= Battery.Cost*Battery.Amount+PVpanel.Cost+ Wiring.Cost+ Wiring.PolesCost*Wiring.PolesAmount+ Electronics.CostA+ Electronics.CostB; 

Load.notservedIndex=find(Load.notserved);
for j=1:num.hours
    Load.demand_all(j==Load.notservedIndex)= false;
end
Energy_Served=sum(Load.demand_all); %This should be the load not served without losses
Total_Energy=num.Years*Energy_Served/1000; %kWh
LCOE=Cost.Total/Total_Energy

LCOE2=Cost.Total/(num.Years*sum(Load.data_all_andLinkConverter)/1000)

%% Plots 

% %%Generation plot 
x_axis=1:1:num.hours;
% p1=plot(x_axis,Generation.Output(x_axis,1));
% xlabel('Hours'); ylabel('PV panel output power (W)'); 
% set(p1,'color', [0.6 0.6 0.6],'LineWidth',1);
% set(gca,'linewidth',1,'fontsize', 20)
% axis([0 num.hours 0 PVpanel.Rating+50]);
% save2pdf('GenerationData.pdf')

%%Load plot
% figure
% x_axis=1:1:num.hours;
% p1=plot(x_axis,Load.data(:,1));
% xlabel('Hours'); ylabel('Load Profile of a single unit (W)'); 
% set(p1,'color', [0.6 0.6 0.6],'LineWidth',1);
% set(gca,'linewidth',1,'fontsize', 15)
% axis([0 num.hours 0 max(Load.data(:,1))+2]);
% save2pdf('LoadProfile.pdf')


plot(Battery.state_charging,'k')
hold on
plot(Battery.state_discharging,'r')
hold on 
plot(Load.data(:,1))
%plot(Load.demand_all)




%hold on 
%plot(x_axis,Load.demand_all)
 
% plot(x_axis,(Load.data_all_andLinkConverter-Load.demand_all))
% %legend('Load including Losses','Load Demand');
% xlabel('Hour (h)');
% ylabel('W');
% figure
%  plot(x_axis,Generation.Temperature(1:num.hours,1))
%  xlabel('Hour (h)');
%  ylabel('Temperature (C)');
 
%   plot(x_axis,Generation.Irradiance(1:num.hours,1))
%  xlabel('Hour (h)');
%  ylabel('Irradiance (W/m2)');
%figure
% plot(x_axis,Load.notserved)
% hold on 
% plot(x_axis,Battery.Spillage)
% legend('Load not served','Battery spillage');
%xlabel('Hour (h)');
% ylabel('W');
% figure
% 
% plot(x_axis,Battery.Energy/12)
% xlabel('Hour (h)');
% ylabel('Ah')


%This script tries to find the lowest cost system design for the required
%reliability
clc 
clear all

%Inputs
Accepted_ReliabilityPer=0.97; %Percentage availability 
CosttoMinimize=0; % 1=Lifetime cost; 0=InitialCost;

num.hours=8760; %24 for 1 day; 8760 for 1 year
num.Households=5;
Battery.MinDODPer=0.5;
Battery.InitialPer=0.5;
Battery.ChargeEff=0.85;
Battery.DischargeEff=0.85;
PVpanel.RatingAtIrradiance=1000;
Wiring.Length=40; 
Wiring.Rho=4e-8;
Wiring.Size=2.5e-6;%m
Network.Voltage=24;
Network.Inefficiency=0.7;


PVpanelRating_Available=150:10:300;
BatteryCapacity_Available=[10:5:100]*12;


%% Generation data

Generation.data=xlsread('pvwatts_hourly_patamda_jharkhand','A20:K8779'); %Averged 10 years of data from NREL
Generation.Irradiance=Generation.data(:,8);
Generation.Temperature=Generation.data(:,6);
Generation.Hour=Generation.data(:,3);

%%Load Data
load.NoLight1Hours = xlsread('demandInputs_uLink_VM.xlsx','Residential', 'NoLight1Hours');
load.NoLight2Hours = xlsread('demandInputs_uLink_VM.xlsx','Residential', 'NoLight2Hours');
load.NoPhoneHours = xlsread('demandInputs_uLink_VM.xlsx','Residential', 'I27:I50');
load.Power = xlsread('demandInputs_uLink_VM.xlsx','Residential', 'Power');


%% Load Profile
[loadprof] = CreateDemand2(Generation,load);
Load.data=loadprof(1:num.hours,1);
Load.demand_all=Load.data*num.Households;
Load.Current=Load.data/Network.Voltage;

%% Load Converter Loss
Converter.data_eff=xlsread('PowerElectronics_Eff','LoadConverter','A2:B11');
Converter.load_power=Converter.data_eff(:,1);
Converter.load_Eff=Converter.data_eff(:,2);
Function_LoadConverter_Eff=polyfit(Converter.load_power,Converter.load_Eff,4);
Function_LoadConverter_Eff_withLoad= polyval(Function_LoadConverter_Eff,Load.data);
Loss.Converter_load=(Load.data).*(100-Function_LoadConverter_Eff_withLoad)/100;
Load.data_andConverter=Load.data+Loss.Converter_load;
% 
% Load.data_andConverter=Load.data;
%% Distribution Loss
Current_data_andCovnerter=Load.data_andConverter/Network.Voltage;
rho_al= 4*10^-8;%Resititive of alum
Wiring.Res=(Wiring.Rho*Wiring.Length)/Wiring.Size;
Loss.Wiring=(Current_data_andCovnerter.^2).*Wiring.Res;
Load.data_andConverter_andWiring=Load.data_andConverter+Loss.Wiring;

%% Link Converter Loss 
Load.data_all=Load.data_andConverter_andWiring*num.Households; %<<<<<<<<<This has to be changed
Converter.data_eff=xlsread('PowerElectronics_Eff','LinkConverter','A2:B15');
Converter.link_power=Converter.data_eff(:,1);
Converter.link_Eff=Converter.data_eff(:,2);
Function_LinkConverter_Eff=polyfit(Converter.link_power,Converter.link_Eff,10);
Function_LinkConverter_Eff_withLoad= polyval(Function_LinkConverter_Eff,Load.data_all);
Loss.Converter_link=(Load.data_all).*(100-Function_LinkConverter_Eff_withLoad)/100;
Load.data_all_andLinkConverter=Load.data_all+Loss.Converter_link;

%% Battery loss
Load.data_all_Battery=Load.data_all_andLinkConverter/0.85;

%% Charge Controller Loss 
Converter.data_eff=xlsread('PowerElectronics_Eff','ChargeController','A2:B17');
Converter.charge_power=Converter.data_eff(:,1);
Converter.charge_Eff=Converter.data_eff(:,2);
Function_chargeConverter_Eff=polyfit(Converter.charge_power,Converter.charge_Eff,4);

%%
Flag=0;
for numBatteries=1:1:length(BatteryCapacity_Available)
    Battery.Capacity=BatteryCapacity_Available(numBatteries);
    Battery.MinDOD=Battery.MinDODPer*Battery.Capacity; %Minimum depth of discharge
    Battery.Initial=Battery.InitialPer*Battery.Capacity;
    
   
    
    for numPV=1:1:length(PVpanelRating_Available)
        PVpanel.Rating=PVpanelRating_Available(numPV);
        Generation.Output= PVpanel.Rating*Generation.Irradiance/PVpanel.RatingAtIrradiance;%Output of solar panel with 1000 W/m2 rating
    
        %Charge Controller Loss 
        Function_chargeConverter_Eff_withGeneration= polyval(Function_chargeConverter_Eff,Generation.Output);
        Loss.Generation_Converter_charge=(Generation.Output).*(100-Function_chargeConverter_Eff_withGeneration)/100;
        Effective.Powertobattery=Generation.Output+Loss.Generation_Converter_charge;
        
        
          %Battery Energy 
          Battery.Energy(1)=Battery.Initial;
          Battery.Spillage=zeros(num.hours,1);
          Load.notserved=zeros(num.hours,1);
          for i=1:num.hours-1
              if Effective.Powertobattery(i)>=Load.data_all_andLinkConverter(i) %Charging
                  Battery.Energy(i+1)= (Effective.Powertobattery(i)-Load.data_all_andLinkConverter(i))*Battery.ChargeEff+Battery.Energy(i);
                  if Battery.Energy(i+1)>Battery.Capacity
                      Battery.Energy(i+1)=Battery.Capacity;
                      Battery.Spillage(i+1)=(Effective.Powertobattery(i)-Load.data_all_andLinkConverter(i));
                  end
              else %Discharging
                  Battery.Energy(i+1)= (Effective.Powertobattery(i)-Load.data_all_andLinkConverter(i))*Battery.DischargeEff+Battery.Energy(i);
                  if Battery.Energy(i+1)<Battery.MinDOD
                      Battery.Energy(i+1)=Battery.MinDOD;
                      Load.notserved(i+1)= (Load.data_all_andLinkConverter(i)-Effective.Powertobattery(i));
                      
                  end
              end
              
          end
          
         Reliability(numBatteries,numPV)=(1-sum(Load.notserved)/sum(Load.data_all_andLinkConverter))*100;
         [TotalCost(numBatteries,numPV), InitialCost(numBatteries,numPV)]=FindCost(BatteryCapacity_Available(numBatteries), PVpanelRating_Available(numPV),Wiring,num,Load);
         if Reliability(numBatteries,numPV)<Accepted_ReliabilityPer*100;
         TotalCost(numBatteries,numPV)=inf;
         InitialCost(numBatteries,numPV)=inf;
         end
         
    end  
end

if CosttoMinimize==1
[MinCost,I]=min(TotalCost(:));
[rowIndex, colIndex]=ind2sub(size(TotalCost),I);
else
 [MinCost,I]=min(InitialCost(:));
[rowIndex, colIndex]=ind2sub(size(InitialCost),I);
end

% 
 fprintf('Selected Battery=%d Ah and PV panel=%d W.\n',BatteryCapacity_Available(rowIndex)/12, PVpanelRating_Available(colIndex))
 fprintf('Initial cost of the system= $%g and lifetime cost=$%g\n',InitialCost(rowIndex,colIndex),TotalCost(rowIndex,colIndex))
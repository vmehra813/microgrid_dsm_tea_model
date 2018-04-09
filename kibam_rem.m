%function [dispatch,savedBalance, battOps] = singlenodesim(hindex,demand,PV,solarCap,BESS,battECpctyNom,genset,gensetIndex,distLossFctr,dieselFCost,plotYesNo)

%%KIBAM Script / Function
%Adapted from REM 
%what does calcEnergyBalance do in REM?



nbatt = battECpctyNom / BESS.energy;%# of batteries in battery bank 

%battECpctyNom = 50; % kWh
% [OLD CODE] battSoCMax = BESS.SOCmax;
% [OLD CODE] battSoCMin = BESS.SOCmin;
% deltaT = 1; %BESS.deltaT; % set equal to 1 (assumes delta T is 1 hour everywhere right now).
% eff_batt_d_c = .89; %BESS.eff_batt_d_c; %Homer assumes same effiency for charging and discharging - sqrt of roundtrip eff
% alpha_c = 1; %BESS.alpha_c; %Max charge rate in A/Ah
% 
% Imax = 16.5; % BESS.Imax; %Max charPVge current in A
% Vnom = 12; %BESS.Vnom;%Nominal voltage in V
% LifeTimeThroughPut = 256; %BESS.LifeTimeThroughPut; %Lifetime Throughput of battery in KWH
% c = .305; %BESS.c; %capacity ratio of battery
% k = 2.12; %BESS.k; %rate constant of battery 
% Qmax_batt = 56.1*Vnom/1000; %BESS.Qmax_batt; % originally given in Amp hours from Homer from battery window - Qmax_batt has been converted to KWH

deltaT = 1; % set equal to 1 (assumes delta T is 1 hour everywhere right now).
eff_batt_d_c = BESS.eff_batt_d_c; %Homer assumes same effiency for charging and discharging - sqrt of roundtrip eff
alpha_c = BESS.alpha_c; %Max charge rate in A/Ah

Imax = BESS.Imax; %Max charge current in A
Vnom = BESS.Vnom;%Nominal voltage in V
LifeTimeThroughPut = BESS.LifeTimeThroughPut; %Lifetime Throughput of battery in KWH
c = BESS.c; %capacity ratio of battery
k = BESS.k; %rate constant of battery 
Qmax_batt = BESS.Qmax_batt; % originally given in Amp hours from Homer from battery window - Qmax_batt has been converted to KWH
Qmax_battbank = Qmax_batt*nbatt; % Qmax_battbank is max capacity of the battery bank in KWH
battECumLim = LifeTimeThroughPut*nbatt; 
e_raised = exp(-k*deltaT);
battECpctyNom = Qmax_battbank;


battPrice = BESS.cost; %100; % $/kWh
%battCapacity = 5;   % kWh
%battCAPEX = battPrice*nbatt;
battCAPEX = battPrice*nbatt;%battCapPrice * battCapacity;  % $
% [OLD CODE] battECumLim = BESS.lifetimeThroughput*nbatt; %kWh --- ***** --- update to read this from catalog...
%Create matrices to store the values for the state variables (these will be
%updated in each interation

Q_1 = zeros(size([demand.HiPty;1])); %some have an added a row bc of t+1 iteration later on
Q_2 = zeros(size([demand.HiPty;1]));
battEcum = zeros(size(demand.HiPty));
battEoutMax = zeros(size([demand.HiPty;1]));
battEinMax = zeros(size([demand.HiPty;1]));
Ebatt_max_in_lim1 = zeros(size(demand.HiPty));
Ebatt_max_in_lim2 = zeros(size(demand.HiPty));
Ebatt_max_in_lim3 = zeros(size(demand.HiPty));
battSOC = zeros(size(demand.HiPty));
Q_total = zeros(size([demand.HiPty;1]));
Power_in_out_batt = zeros(size([demand.HiPty;1]));
%these are the starting values of Q_2 for the first iteration - this is a property of the battery bank 
%define the variable values in the first iteration
Initial_battSOC = BESS.SOCinit;
Q_1(1) = c*Initial_battSOC*Qmax_battbank; %these are the starting values of Q_1 and Q_2 for the first iteration - this is a property of the battery bank 
Q_2(1) = (1-c)*Initial_battSOC*Qmax_battbank; 

battEcum(1) = 0;
Q_total(1) = Q_1(1) + Q_2(1);
battSOC(1) = Q_total(1)/(Qmax_battbank);
%max_cap*Vnom/1000 converts max_cap, which is in Amp*hours to KWH

%Max discharge Energy per hour - how do you tell that this is for the first
%iteration?
% Ebatt_dmax_kbm and Ebatt_cmax_kbm is the energy that can be discharged
% and charged at any deltaT into and out of the battery BANK
%Ebatt_dmax_kbm(t) 
%Both have already taken into account efficiency
%battEoutMax(1)= -1*eff_batt_d_c*deltaT*((-k*c*Qmax_battbank + k*Q_1(1)*e_raised + Q_total(1)*k*c*(1-e_raised)/...
    %(1-e_raised+c*(k*deltaT-1+e_raised))));
EoutMax_kineticmodellimit = -1*eff_batt_d_c*deltaT*((-k*c*Qmax_battbank + k*Q_1(1)*e_raised + Q_total(1)*k*c*(1-e_raised))/...
    (1-e_raised+c*(k*deltaT-1+e_raised)));
EoutMax_SOCminlimit = (battSOC(1)-BESS.SOCmin)*Qmax_battbank;
battEoutMax(1)= min(EoutMax_kineticmodellimit,EoutMax_SOCminlimit);



%Max charge Energy per hour - the minimum of three limitations for the
%battery bank 
%Ebatt_cmax_kbm(t) is battEinMax - this is the maximum amount of energy
%that the battery can charge in the hour (t)
%Ebatt_max_in_lim1(1) = deltaT*(k*Q_1(1)*e_raised+Q_total(1)*k*c*(1-e_raised))/(1-e_raised+c*(k*deltaT-1+e_raised));
%Ebatt_max_in_lim2(1) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-Q_total(1))/deltaT;
%Ebatt_max_in_lim3(1) = deltaT*nbatt*Imax*Vnom/1000;

%battEinMax(1) = min([Ebatt_max_in_lim1(1), Ebatt_max_in_lim2(1), Ebatt_max_in_lim3(1)])/eff_batt_d_c;

Ebatt_max_in_lim1(1) = deltaT*k*Q_1(1)*e_raised+Q_total(1)*k*c*(1-e_raised)/(1-e_raised+c*(k*deltaT-1+e_raised));
Ebatt_max_in_lim2(1) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-Q_total(1))/deltaT;
Ebatt_max_in_lim3(1) = deltaT*nbatt*Imax*Vnom/1000;
Ebatt_max_in_limSOCmax = (BESS.SOCmax-battSOC(1))*Qmax_battbank;

battEinMax(1) = min([Ebatt_max_in_lim1(1), Ebatt_max_in_lim2(1), Ebatt_max_in_lim3(1), Ebatt_max_in_limSOCmax])/eff_batt_d_c;

% battLossIn %= @(battEIn,battSoC) battEIn*battEffInNom;
% battLossOut %= @(battEOut,battSoC) battEOut*battEffOutNom;
battEInMaxNom = 0.198*nbatt;%BESS.maxchargeratenom * nbatt; % adjust later for variable timestep
battEOutMaxNom = 6.6*nbatt;%BESS.maxdischargeratenom * nbatt; % adjust later for variable timestep
% [OLD CODE] f_battEInMax = @(battSoC) min(battECpctyNom*(battSoCMax-battSoC),...
    % [OLD CODE]battEInMaxNom*(1-BESS.IO*(exp((1-battSoC)/BESS.NKT)-1)));  % = @(battSoC)
% [OLD CODE] f_battEOutMax = @(battSoC) min(battECpctyNom*(battSoC-battSoCMin),...
    % [OLD CODE]battEOutMaxNom*(1-BESS.IO*(exp((battSoC)/BESS.NKT)-1))); % = @(battSoC)

    %Updating Q1 and Q2 equations at the end of of the time step
%NEED move this further down
% P in the equation = ((battEinMax(t)-battEoutMax(t))/deltaT) and refers to
% the power (net)[KW) that is charged or discharged into battery
% Q_1 is the energy in tank 1 (kwH) available energy at the end of the
% iteration (t), so it will be the starting Q_1 for the next iteration
%Q_1(t+1) = Q_1(t)*E^(-k*deltaT)+(Q_total(t)*k*c-((battEinMax(t)-battEoutMax(t))/deltaT)*(1-e^-k*deltaT)/k+P*c*(k*deltaT-1+e^(-k*deltaT)/k;
%Q_2(t+1) = Q_2(t)*E^(-k*deltaT)+(Q_total(t)*(1-c))*(1-e^-k*deltaT)+((battEinMax(t)-battEoutMax(t))/deltaT)*(1-c)*(k*deltaT-1+e^(-k*deltaT)/k;




if battECpctyNom > 0
         %Q_total(t) = Q_1(t) + Q_2(t);
         
         %11.20 - add this back in if it doesn't work!! battSOC(t+1) = battSOC(t) + (dispatch.battin(t)/Qmax_battbank) - (dispatch.battout(t)/Qmax_battbank);
         
%        [OLD CODE] battEOutMax(t+1) = f_battEOutMax(battSoC(t+1));
%        [OLD CODE]  battEInMax(t+1) = f_battEInMax(battSoC(t+1));
        battEcum(t+1) = battEcum(t) + dispatch.battin(t) + dispatch.battout(t); %is this right, or should we should battery losses count too (or just losses into battery)?
   

        
    %Updating Q1 and Q2 equations at the end of of the time step
    %NEED move this further down
    % E_in_out_batt = ((battEinMax(t)-battEoutMax(t))/deltaT) and refers to
    % the power (net)[KW) that is charged or discharged into batt bank(?)
    % Q_1(t+1) is the energy in tank 1 (kwH) available energy at the end of the
    % iteration (t), so it will be the starting Q_1 for the next iteration
    % (Q_1(t+1))
    %E_in_out_batt = ((battEinMax(t)-battEoutMax(t))/deltaT);
    %e_raised = exp(-k*deltaT);
    %Q_1(t+1) = Q_1(t)*e_raised+Q_total(t)*k*c-E_in_out_batt*(1-e_raised)/k+E_in_out_batt*c*(k*deltaT-1+e_raised)/k;
    %Q_2(t+1) = Q_2(t)*e_raised+(Q_total(t)*(1-c))*(1-e_raised)+E_in_out_batt*(1-c)*(k*deltaT-1+e_raised)/k;
    Power_in_out_batt(t) = ((dispatch.battin(t)- dispatch.battout(t))/deltaT);     
    e_raised = exp(-k*deltaT);
    
    %error in Q_1 - fixed
    %should add in charge discharge efficiencies -- if charge, then *
    %power_in_outbatt by charge eff etc. -- check if it's been accounted
    %for elsewhere
    Q_1(t+1) = Q_1(t)*e_raised+(Q_total(t)*k*c+Power_in_out_batt(t))*(1-e_raised)/k+Power_in_out_batt(t)*c*(k*deltaT-1+e_raised)/k;
    Q_2(t+1) = Q_2(t)*e_raised+(Q_total(t)*(1-c))*(1-e_raised)+Power_in_out_batt(t)*(1-c)*(k*deltaT-1+e_raised)/k;
    Q_total(t+1) = Q_1(t+1) + Q_2(t+1);
    battSOC(t+1) = Q_total(t+1)/Qmax_battbank;
   %this is battSOC at the beginning of battSOC(t+1) timestep
   % battEoutMax and battEinMax represent the max energy that can be discharged
    % and charged at any deltaT into and out of the battery BANK
    %Both have already taken into account efficiency
    %battEoutMax(t+1)= -1*eff_batt_d_c*deltaT*((-k*c*Qmax_battbank + k*Q_1(t+1)*e_raised + Q_total(t+1)*k*c*(1-e_raised)/...
    %(1-e_raised+c*(k*deltaT-1+e_raised))));
    
%     battEoutMax(t+1)= -1*eff_batt_d_c*deltaT*((-k*c*Qmax_battbank + k*Q_1(t+1)*e_raised + Q_total(t+1)*k*c*(1-e_raised))/...
%         (1-e_raised+c*(k*deltaT-1+e_raised)));
    
    EoutMax_kineticmodellimit = -1*eff_batt_d_c*deltaT*((-k*c*Qmax_battbank + k*Q_1(t+1)*e_raised + Q_total(t+1)*k*c*(1-e_raised))/...
        (1-e_raised+c*(k*deltaT-1+e_raised)));
    EoutMax_SOCminlimit = (battSOC(t+1)-BESS.SOCmin)*Qmax_battbank;
    battEoutMax(t+1)= min(EoutMax_kineticmodellimit,EoutMax_SOCminlimit);
    
    %Max charge Energy per hour - the minimum of three limitations for the
    %battery bank 
    %battEinMax - this is the maximum amount of energy
    %that the battery bank can charge in the hour (t)
    %Ebatt_max_in_lim1(t+1) = deltaT*(k*Q_1(t+1)*e_raised+Q_total(t+1)*k*c*(1-e_raised))/(1-e_raised+c*(k*deltaT-1+e_raised));
    %Ebatt_max_in_lim2(t+1) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-Q_total(t+1))/deltaT;
    %Ebatt_max_in_lim3(t+1) = deltaT*nbatt*Imax*Vnom/1000;
    %battEinMax(t+1) = min([Ebatt_max_in_lim1(t+1), Ebatt_max_in_lim2(t+1), Ebatt_max_in_lim3(t+1)])/eff_batt_d_c;
    Ebatt_max_in_lim1(t+1) = deltaT*(k*Q_1(t+1)*e_raised+Q_total(t+1)*k*c*(1-e_raised))/(1-e_raised+c*(k*deltaT-1+e_raised));
    Ebatt_max_in_lim2(t+1) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-Q_total(t+1))/deltaT;
    Ebatt_max_in_lim3(t+1) = deltaT*nbatt*Imax*Vnom/1000;
    Ebatt_max_in_limSOCmax = (BESS.SOCmax-battSOC(t+1))*Qmax_battbank;
    battEinMax(t+1) = min([Ebatt_max_in_lim1(t+1), Ebatt_max_in_lim2(t+1), Ebatt_max_in_lim3(t+1), Ebatt_max_in_limSOCmax])/eff_batt_d_c;
  
   
    %what does this mean?
    % dynamic adjustment to batt value
        imax = size(srtdDispStackNoBatt,1);
        %battSOCFrac = (battSoC(t)-BESS.SOCmin)/(BESS.SOCmax-BESS.SOCmin);
        battSOCFrac = (BESS.SOCmax-battSOC(t+1))/(BESS.SOCmax-BESS.SOCmin);
        i = ceil(imax*battSOCFrac);
        if i == 0
            i = 1;
        end
        if i > imax
           i = imax; 
        end
        vBattOutAdj(t+1) = srtdDispStackNoBatt{i,1}*(battOutEffNom*battPEOutEffNom*invrtrEffNom*distEffNom)...
            +battCAPEX/battECumLim;
    end








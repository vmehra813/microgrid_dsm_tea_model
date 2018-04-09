function [dispatch,savedBalance, battOps] = singlenodesim(hindex,demand,PV,solarCap,BESS,battECpctyNom,genset,gensetIndex,distLossFctr,dieselFCost,plotYesNo)
%Current version 11/15/2014_11:15AM 
%%%TESTING - VIVIAN 10_05_2015
%%%Vivian - Added pmin limit for the generator
dbstop if error

% load(fullfile('C:\REM\runs\Kirambi_20150220\input\region1','genCatalog.mat'));
% load(fullfile('C:\REM\runs\Kirambi_20150220\input\region1', 'settings.mat'));
% load(fullfile('C:\REM\runs\Kirambi_20150220\input\region1','demandLibrary.mat'));
% load(fullfile('C:\REM\runs\Kirambi_20150220\input\region1', 'geodata.mat'));
global genCatalog
global settings
global geodata


plotYesNo = 0;

%Battery parameters are saved under genCatalog - solvegenCatalog writes raw
%info into this. Under genCatalog, there are n batteries, saved
%under BESS(1)... BESS(n)
%Below -- simply assigning variables values which are saved under
%genCatalog
%% Parameters and functions not changing over time

%% Select subset hours of year
demand.HiPty = demand.HiPty(hindex);
demand.LoPty = demand.LoPty(hindex);

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

%This is defined later in the script
% battEcum(t+1) = battEcum(t) + dispatch.battin(t) + dispatch.battout(t);
%dieselFCost = settings.dieselFCost; %1; % $/l  % -- pull this from elsewhere?
%dieselED = 10; % kWh/l  % -- pull this from elsewhere?

%is this (below) where the variables, values, structures, etc. are updated?
genEMax = genset.pMax(gensetIndex); % adjust later for variable timestep
genEMin = genset.pMin(gensetIndex);
%genEff = @(genE) .3; %= @(genE) genEffNom;
genFuelCostPerKWHFun = @(genE) polyval(genset.LPerKWHCurve(gensetIndex,:),genE)*dieselFCost;

%genEff = @(genE) .3; %= @(genE) genEffNom; % adjust to pull variable genset eff...
%genSUF = genset.startupFuel(gensetIndex);



%distLossFctr

% %% Select subset hours of year
% demand.HiPty = demand.HiPty(hindex);
% demand.LoPty = demand.LoPty(hindex);


%% State Tracking:
solarE = geodata.resources.solar.sFactorDC(hindex) * solarCap;
%solarE = zeros(24,1);

% [OLD CODE] battSoC = zeros(size(demand.HiPty)); %add 1 to some of these??
% [OLD CODE] battEOutMax = zeros(size(demand.HiPty));% zeros(24,1);
% [OLD CODE] battEInMax = zeros(size(demand.HiPty)); %zeros(24,1);
% [OLD CODE] battECum = zeros(size(demand.HiPty)); %zeros(24,1);

dispatch.solar = zeros(size(demand.HiPty)); %zeros(24,1);
dispatch.gen = zeros(size(demand.HiPty)); %zeros(24,1);
dispatch.crtllo = zeros(size(demand.HiPty)); %zeros(24,1);
dispatch.crtlhi = zeros(size(demand.HiPty)); %zeros(24,1);
dispatch.battout = zeros(size(demand.HiPty)); %zeros(24,1);
dispatch.battin = zeros(size(demand.HiPty)); %zeros(24,1);
dispatch.genFlag = zeros(size(demand.HiPty)); %zeros(24,1);

%demand.HiPty = zeros(24,1);
%demand.LoPty = zeros(24,1);

savedBalance = zeros(size(demand.HiPty)); %zeros(24,1);
savedLosses = zeros(size(demand.HiPty)); %zeros(24,1);
savedNonGridLosses = zeros(size(demand.HiPty));

% dispatchinit.solar = zeros(24,1);
% dispatchinit.gen = zeros(24,1);
% dispatchinit.crtllo = zeros(24,1);
% dispatchinit.crtlhi = zeros(24,1);
% dispatchinit.battout = zeros(24,1);
% dispatchinit.battin = zeros(24,1);
% 
% savedBalanceInit= zeros(24,1);
% savedLossesInit = zeros(24,1);

%% Data
%solarE(6:16) = 10;

%demand.HiPty = 3*rand(size(demand.HiPty));
%demand.LoPty = 5*rand(size(demand.LoPty));

% [OLD CODE] battSoC(1) = BESS.SOCinit;
% OLD CODE] battEOutMax(1) =  f_battEOutMax(battSoC(1));
% [OLD CODE] battEInMax(1) =  f_battEInMax(battSoC(1));
% [OLD CODE] battECum(1) = 0;

dispatch.genFlag(1) = 1;

%% State Definition:
% solarE(t)
%
% battSoC(t)
% battEOutMax(t)
% battEInMax(t)
% battECum(t)
%
% genFlag(t)
%
% dmndEHiPty(t)
% dmndELoPty(t)

% Dispatch/UC/BatterySchedule Decisions

% DbattEOut(t)
% DbattEIn(t)
% DgenEOut(t)
% DdmndEHiPtyCrtl(t)
% DdmndELoPtyCrtl(t)
% DsolarESpill(t)
availGen = ones(size(demand.HiPty)); %ones(24,1); %@(SOC) or received from "UC"
reqGen = zeros(size(demand.HiPty));%ones(size(demand.HiPty)); %ones(24,1);   %= 0 or 1  %go back and figure out gen on/off availability vs requirement

% ---- adjust for losses
vBattOutAdj = zeros(size(demand.HiPty)); %zeros(24,1);

% ---- adjust for losses to grid
cSolar = 0; %from solar VC
%cGenHiLoad = .2; %from gen eff curve
cGenMidLoad = genFuelCostPerKWHFun((genEMax+genEMin)/2);
%dieselFCost/(dieselED*genEff((genEMax+genEMin)/2)); %from gen eff curve
%cGenLoLoad = .3; %from gen eff curve
cCrtlLoPty = settings.penalty(1) + settings.price(1); %from CNSE
cCrtlHiPty = settings.penalty(2) + settings.price(2); %from CNSE

solarPEEffNom = genCatalog.PV_PE.effNom; %0.95;
battPEOutEffNom = genCatalog.BESS_PE.outEffNom; %0.95;
battPEInEffNom = genCatalog.BESS_PE.inEffNom; %0.95;
% battInEffNom = BESS.chargeeff; %0.9; CHANGE THIS BACK IF CREATEGENCATALOG
% FOR NEW BATTERY DOESNT WORK
% battOutEffNom = BESS.dischargeeff; %0.9;
battInEffNom = BESS.eff_batt_d_c; %0.9;
battOutEffNom = BESS.eff_batt_d_c; %0.9;
invrtrEffNom = genCatalog.ACDCConverter.invrtrEffNom; %0.9;
rctfrEffNom = genCatalog.ACDCConverter.rctfrEffNom; %0.9;
%distLossFctr = settings.distLossFctr; %0.05;

solarPELossF = 1 - solarPEEffNom;
battPEOutLossF = 1 - battPEOutEffNom;
battPEInLossF = 1 - battPEInEffNom;
battInLossF = 1 - battInEffNom;
battOutLossF = 1 - battOutEffNom;
invrtrLossF = 1 - invrtrEffNom;
rctfrLossF = 1 - rctfrEffNom;
distEffNom = 1 - distLossFctr;

% cost of sources to grid
cSolarAdj = cSolar/(solarPEEffNom*invrtrEffNom*distEffNom);
cGenMidLoadAdj = cGenMidLoad/(distEffNom);
cCrtlLoPtyAdj = cCrtlLoPty; 
cCrtlHiPtyAdj = cCrtlHiPty;
%vBattOutAdj = vBattOut(:).*(battOutEffNom*battPEOutEffNom*invrtrEffNom*distEffNom);

% cost of sources to charge battery
cSolarAdj2 = cSolar/(solarPEEffNom*battPEInEffNom*battInEffNom) + battCAPEX/battECumLim;
cGenMidLoadAdj2 = cGenMidLoad/(rctfrEffNom*battPEInEffNom*battInEffNom) + battCAPEX/battECumLim;
cCrtlLoPtyAdj2 = cCrtlLoPty/(distEffNom*rctfrEffNom*battPEInEffNom*battInEffNom) + battCAPEX/battECumLim;
cCrtlHiPtyAdj2 = cCrtlHiPty/(distEffNom*rctfrEffNom*battPEInEffNom*battInEffNom) + battCAPEX/battECumLim;

    if availGen(1) && genEMax > 0 %genEMax referring to generator (genEMax > 0 means if there is a generator?) - this is for t=1 
        dispStackNoBatt = {cSolarAdj,'a_solar';
            cGenMidLoadAdj, 'b_gen';
            cCrtlLoPtyAdj, 'c_dlo'; %cost of non served lo demand?
            cCrtlHiPtyAdj, 'd_dhi'}; %cost of non served high demand?
    else
        dispStackNoBatt = {cSolarAdj,'a_solar';
            cCrtlLoPtyAdj, 'c_dlo';
            cCrtlHiPtyAdj, 'd_dhi'};
    end
    srtdDispStackNoBatt = sortrows(dispStackNoBatt); %starting srtDispStack?
    
    % dynamic adjustment to batt value - first alue is assigned here
    if battECpctyNom ~= 0
        imax = size(srtdDispStackNoBatt,1);
        battSOCFrac = (battSOC(1)-BESS.SOCmin)/(BESS.SOCmax-BESS.SOCmin);
        i = ceil(imax*battSOCFrac);
        vBattOutAdj(1) = srtdDispStackNoBatt{i,1}*(battOutEffNom*battPEOutEffNom*invrtrEffNom*distEffNom)...
            +battCAPEX/battECumLim;
    end

for t = 1:size(solarE,1)
    %% Establish dispatch priority
    
    maxed.solar = 0;
    maxed.gen = 0;
    maxed.crtllo = 0;
    maxed.crtlhi = 0;
    maxed.battout = 0;
    
    if availGen(t) && genEMax > 0 %dispStack = costs to charge grid (aka meet demand)??
        dispStack = {cSolarAdj,'a_solar';
            cGenMidLoadAdj, 'b_gen';
            cCrtlLoPtyAdj, 'c_dlo';
            cCrtlHiPtyAdj, 'd_dhi';
            vBattOutAdj(t), 'e_battout'};
    else %without generator?
        dispStack = {cSolarAdj,'a_solar'; 
            cCrtlLoPtyAdj, 'c_dlo';
            cCrtlHiPtyAdj, 'd_dhi';
            vBattOutAdj(t), 'e_battout'};
    end
    srtdDispStack = sortrows(dispStack); %sorted by increasing costs to charge grid??
    
    if availGen(t) && genEMax > 0 % dispStack2 = costs to charge battery??
        dispStack2 = {cSolarAdj2,'a_solar';
            cGenMidLoadAdj2, 'b_gen';
            cCrtlLoPtyAdj2, 'c_dlo';
            cCrtlHiPtyAdj2, 'd_dhi';
            vBattOutAdj(t), 'e_battout'};
    else
        dispStack2 = {cSolarAdj2,'a_solar';
            cCrtlLoPtyAdj2, 'c_dlo';
            cCrtlHiPtyAdj2, 'd_dhi';
            vBattOutAdj(t), 'e_battout'};
    end % note generator cannot be used to charge battery???
    srtdDispStack2 = sortrows(dispStack2);
    
    %first run to meet grid demand (srtdDispStack), next will run with
    %remaining energy to charge battery (srtdDispStack2)
    for j = 1:size(srtdDispStack,1) %going sequentially through each avail resource from cheapest to most expensive
        if (j == 1) && reqGen(t) %if only generator use ?
            dispatch.gen(t) = genEMin;
            [balance,~] = calcEnergyBalance(dispatch,demand,t,BESS,distLossFctr);
            if balance >=0
                lastin = 'b_gen';
                break
            end
        end
        
        % if the cheapest cost is a, then dispatch will be solar avial at
        % time t (later if there is excess demand not met, will figure out
        % which resource next to use (remember dispatch order is decided
        % above)
        if srtdDispStack{j,2}(1) == 'a' %refers to first letter 
            dispatch.solar(t) = solarE(t);
            maxed.solar = 1;
        elseif srtdDispStack{j,2}(1) == 'b'
            dispatch.gen(t) = genEMax;
            maxed.gen = 1;
        elseif srtdDispStack{j,2}(1) == 'c'
            dispatch.crtllo(t) = demand.LoPty(t);
            maxed.crtllo = 1;
        elseif srtdDispStack{j,2}(1) == 'd'
            dispatch.crtlhi(t) = demand.HiPty(t);
            maxed.crtlhi = 1;
        elseif srtdDispStack{j,2}(1) == 'e'
            dispatch.battout(t) = battEoutMax(t); 
            maxed.battout = 1;
        end
        
        %check if dispatch can cover demand. If yes, break.  If no, try
        %next dispatch option in the stacks.
        [balance,~] = calcEnergyBalance(dispatch,demand,t,BESS,distLossFctr);
        if balance >=-1e-10 %used this instead of 0 because sometimes math is slightly off
            lastin = srtdDispStack{j,2};
            break
        else
            continue
        end
    end
    
    % balance dispatch - this adjusts the dipatched energy to match the
    % remainder of energy that is needed
    switch lastin %figuring out which resource was the last one used to meet demand in time t (for grid)
        case 'a_solar'
            adj = balance/(invrtrEffNom*solarPEEffNom);
            dispatch.solar(t) = dispatch.solar(t) - adj;
            maxed.solar = 0;
        case 'b_gen'
            adj = balance;
            dispatch.gen(t) = max((dispatch.gen(t) - adj), genEMin); %%%%% ADDED THIS
            maxed.gen = 0;
        case 'c_dlo'
            adj = balance/(1+distLossFctr);
            dispatch.crtllo(t) = dispatch.crtllo(t) - adj;
            maxed.crtllo = 0;
        case 'd_dhi'
            adj = balance/(1+distLossFctr);
            dispatch.crtlhi(t) = dispatch.crtlhi(t) - adj;
            maxed.crtlhi = 0;
        case 'e_battout'
            adj = balance/(battOutEffNom*battPEOutEffNom*invrtrEffNom);
            dispatch.battout(t) = dispatch.battout(t) - adj; %taking into consideration losses
            maxed.battout = 0;
    end
    
    balanced = 1;
    
%     dispatchinit.solar(t) = dispatch.solar(t);
%     dispatchinit.gen(t) = dispatch.gen(t);
%     dispatchinit.crtllo(t) = dispatch.crtllo(t);
%     dispatchinit.crtlhi(t) = dispatch.crtlhi(t);
%     dispatchinit.battout(t) = dispatch.battout(t);
%     dispatchinit.battin(t) = dispatch.battin(t);
%     
%     [savedBalanceInit(t),savedLossesInit(t)] = calcEnergyBalance(dispatch,demand,t,BESS);
    
%now deciding whether or not to charge battery????
    for j = 1:size(srtdDispStack2,1)
        % check if next option is cheaper than battery value
        if srtdDispStack2{j,1} < vBattOutAdj(t)
            %if it is, check if it is maxed. If maxed, continue
            switch srtdDispStack2{j,2}
                case 'a_solar'
                    if maxed.solar
                        continue
                    else
                        dispatch.solar(t) = solarE(t);
                        maxed.solar = 1;
                    end
                case 'b_gen'
                    if maxed.gen
                        continue
                    else
                        dispatch.gen(t) = genEMax;
                        maxed.gen = 1;
                    end
                case 'c_dlo'
                    if maxed.crtllo
                        continue
                    else
                        dispatch.crtllo(t) = demand.LoPty(t);
                        maxed.crtllo = 1;
                    end
                case 'd_dhi'
                    if maxed.crtlhi
                        continue
                    else
                        dispatch.crtlhi(t) = demand.HiPty(t);
                        maxed.crtlhi = 1;
                    end
                case 'e_battout'
                    if maxed.battout
                        continue
                    else
                        dispatch.battout(t) = battEOutMax(t);
                        maxed.battout = 1;
                    end
            end
            
            % If not maxed, max it and max batt in.  Check balance.
            dispatch.battin(t) = battEinMax(t);
            [balance,~,totalDCOut] = calcEnergyBalance(dispatch,demand,t,BESS,distLossFctr);
            balanced = 0;
            % If balance neg, set last in to batt in; continue
            % If balance pos, set last in to jth item; break
            if balance <=0
                lastin = 'battin';
                continue
            else
                lastin = srtdDispStack2{j,2};
                break
            end
        else
            % If next dispatch item is more expensive than batt value,
            % break
            break
        end
    end
    
    % balance dispatch
    if ~balanced
        switch lastin
            case 'a_solar'
%                 adj = balance/(invrtrEffNom*solarPEEffNom);
%                 dispatch.solar(t) = dispatch.solar(t) - adj;
                
                oldtotalDCOut = totalDCOut;
                count = 0;
                while 1
                    if oldtotalDCOut >= 0
                        adj = balance/(invrtrEffNom*solarPEEffNom);
                        dispatch.solar(t) = dispatch.solar(t) - adj;
                        [balance, ~, newtotalDCOut] = calcEnergyBalance(dispatch,demand,t,BESS,distLossFctr);
                    else
                        adj = balance*(rctfrEffNom/solarPEEffNom);
                        dispatch.solar(t) = dispatch.solar(t) - adj;
                        [balance, ~, newtotalDCOut] = calcEnergyBalance(dispatch,demand,t,BESS,distLossFctr);
                    end
                    if abs(balance) < 1e-10
                        break
                    end
                    oldtotalDCOut = newtotalDCOut;
                    count = count+1;
                    if count > 10
                        error('solar in not converging');
                    end
                end
                
            case 'b_gen'
                adj = balance;
                dispatch.gen(t) = max((dispatch.gen(t) - adj), genEMin); %%% ADDED THIS
            case 'c_dlo'
                adj = balance/(1+distLossFctr);
                dispatch.crtllo(t) = dispatch.crtllo(t) - adj;
            case 'd_dhi'
                adj = balance/(1+distLossFctr);
                dispatch.crtlhi(t) = dispatch.crtlhi(t) - adj;
            case 'e_battout'
                adj = balance/(battOutEffNom*battPEOutEffNom*invrtrEffNom);
                dispatch.battout(t) = dispatch.battout(t) - adj;
            case 'battin'
                oldtotalDCOut = totalDCOut;
                count = 0;
                while 1
                    if oldtotalDCOut >= 0
                        adj = -balance*((battInEffNom*battPEInEffNom)/invrtrEffNom);
                        dispatch.battin(t) = dispatch.battin(t) - adj;
                        [balance, ~, newtotalDCOut] = calcEnergyBalance(dispatch,demand,t,BESS,distLossFctr);
                        %                         if newtotalDCOut >= 0
                        %
                        %                             break
                        %                         end
                    else
                        adj = -balance*(battInEffNom*battPEInEffNom*rctfrEffNom);
                        dispatch.battin(t) = dispatch.battin(t) - adj;
                        [balance, ~, newtotalDCOut] = calcEnergyBalance(dispatch,demand,t,BESS,distLossFctr);
                        %                         if newtotalDCOut < 0
                        %                             break
                        %                         end
                    end
                    if abs(balance) < 1e-10
                        break
                    end
                    %                     if oldtotalDCOut - newtotalDCOut < 1e-10
                    %                         break
                    %                     end
                    oldtotalDCOut = newtotalDCOut;
                    count = count+1;
                    if count > 10
                        error('batt in not converging');
                    end
                end
        end
    end
    
     %% NEW PART
    if strcmp(lastin,'b_gen') && (dispatch.gen(t)<genEMin && dispatch.gen(t)>0)
        dispatch_old=dispatch;
        %% OPTION 1 - P_gen=0
            dispatch.gen(t)=0;
        
            
            [~, index] = ismember('b_gen',srtdDispStack(:,2))      ;
             cost_option_1=0;

                for j = index+1:size(srtdDispStack,1)
                    lastin = srtdDispStack(j,2);
                               if strcmp(lastin,'c_dlo')
                                       [balance, ~, ~] = calcEnergyBalance(dispatch,demand,t,BESS,distLossFctr);
                                       adj = balance/(1+distLossFctr);
                                       dispatch.crtllo(t) = min(dispatch.crtllo(t) - adj,demand.LoPty(t));
                                       
                                       [~, index_temp] = ismember('c_dlo',srtdDispStack(:,2));  
                                       cost_option_1=cost_option_1+dispatch.crtllo(t)*cell2mat(srtdDispStack(index_temp,1));
                               end
                               if strcmp(lastin,'d_dhi')
                                       [balance, ~, ~] = calcEnergyBalance(dispatch,demand,t,BESS,distLossFctr);
                                       adj = balance/(1+distLossFctr);
                                       dispatch.crtlhi(t) = min(dispatch.crtlhi(t) - adj,demand.LoPty(t));
                                      
                                        [~, index_temp] = ismember('d_dhi',srtdDispStack(:,2));
                                        cost_option_1=cost_option_1+dispatch.crtlhi(t)*cell2mat(srtdDispStack(index_temp,1));
                               end
                                if strcmp(lastin,'e_battout')
                                      [balance, ~, ~] = calcEnergyBalance(dispatch,demand,t,BESS,distLossFctr);
                                       adj = balance/(battOutEffNom*battPEOutEffNom*invrtrEffNom);
                                       dispatch.battout(t) = min(dispatch.battout(t) - adj, battEoutMax(t));
                                       
                                       [~, index_temp] = ismember('e_battout',srtdDispStack(:,2));  
                                       cost_option_1=cost_option_1+dispatch.battout(t)*cell2mat(srtdDispStack(index_temp,1));
                               end
                end
                
                dispatch_option1=dispatch;
        
        %% OPTION 2 - P_gen=genEMIN
                dispatch=dispatch_old;
                dispatch.gen(t)=genEMin;
            
                cost_option_2=dispatch.gen(t)*cell2mat(srtdDispStack(index,1));
                
                lastin='e_battout';
                [balance, ~, ~] = calcEnergyBalance(dispatch,demand,t,BESS,distLossFctr);

                adj = balance/(battOutEffNom*battPEOutEffNom*invrtrEffNom);
                
                
                if dispatch.battout(t) - adj<0;
                lastin='battin';
                dispatch.battout(t)=0;
                oldtotalDCOut = totalDCOut;
                count = 0;
                while 1
                    if oldtotalDCOut >= 0
                        adj = -balance*((battInEffNom*battPEInEffNom)/invrtrEffNom);
                        dispatch.battin(t) = dispatch.battin(t) - adj;
                        [balance, ~, newtotalDCOut] = calcEnergyBalance(dispatch,demand,t,BESS,distLossFctr);
                        %                         if newtotalDCOut >= 0
                        %
                        %                             break
                        %                         end
                    else
                        adj = -balance*(battInEffNom*battPEInEffNom*rctfrEffNom);
                        dispatch.battin(t) = dispatch.battin(t) - adj;
                        [balance, ~, newtotalDCOut] = calcEnergyBalance(dispatch,demand,t,BESS,distLossFctr);
                        %                         if newtotalDCOut < 0
                        %                             break
                        %                         end
                    end
                    if abs(balance) < 1e-10
                        break
                    end
                    %                     if oldtotalDCOut - newtotalDCOut < 1e-10
                    %                         break
                    %                     end
                    oldtotalDCOut = newtotalDCOut;
                    count = count+1;
                    if count > 10
                        error('batt in not converging');
                    end
                end
                
                if dispatch.battin(t)>battEinMax(t)
                    spilled(t)=dispatch.battin(t)-battEinMax(t);
                    dispatch.battin(t)=battEinMax(t);
                    
                   [~, index_temp] = ismember('b_gen',srtdDispStack2(:,2));  
                    cost_option_2=dispatch.battin(t)*cell2mat(srtdDispStack2(index_temp,1));

                end
                
                else
                    dispatch.battout(t)=dispatch.battout(t) - adj;
                    
                    [~, index_temp] = ismember('e_battout',srtdDispStack(:,2));  
                    cost_option_2=cost_option_2-adj*cell2mat(srtdDispStack(index_temp,1));
                end
                
                dispatch_option2=dispatch;
                
                if cost_option_1<cost_option_2
                    dispatch=dispatch_option1;
                end
                
                 
    end

    
    %% Save Energy Balance
    [savedBalance(t),savedLosses(t),~,savedNonGridLosses(t)] = calcEnergyBalance(dispatch,demand,t,BESS,distLossFctr);
    
    %% Update State Variables
    %solarE(t+1)
    
    dispatch.genFlag(t+1) = (dispatch.gen(t)>0);
    
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
%     valbatt = vBattOut(t+1)
%     soc = battSoC(t+1)
%     srtdDispStack2
    
    %dmndEHiPty(t+1)
    %dmndELoPty(t+1)
    
    % Dispatch Decisions
    %DbattEOut(t+1)
    %DbattEIn(t+1)
    %DgenEOut(t+1)
    %DdmndEHiPtyCrtl(t+1)
    %DdmndELoPtyCrtl(t+1)
    
    
    %%%INSERT SIMULATION HERE
end
% % area([-savedLosses,-dispatch.battin,dispatch.battout,dispatch.solar,dispatch.gen,dispatch.crtllo,dispatch.crtlhi]);
% % M = {'losses';'battin';'battout';'solar';'gen';'crtllo';'crtlhi'};

dispatch.losses = savedLosses;
dispatch.nonGridLosses = savedNonGridLosses;
dispatch.solarAvail = solarE;

if plotYesNo

    figure(2);
    %area([zeros(size(dispatch.battout)),dispatch.battout,dispatch.solar,dispatch.gen,dispatch.crtllo,dispatch.crtlhi,solarE-dispatch.solar]);
    
    areaPlot = area([zeros(size(dispatch.battout)),dispatch.battout,dispatch.solar,dispatch.gen,dispatch.crtllo+dispatch.crtlhi,solarE-dispatch.solar]);
    
    %set(areaPlot(1),'FaceColor',);
    %set(areaPlot(2),'FaceColor',);
    set(areaPlot(3),'FaceColor','y');
    set(areaPlot(4),'FaceColor','g');
    set(areaPlot(5),'FaceColor','r');
    set(areaPlot(6),'FaceColor',[.5 .5 .5]);
    %bar([zeros(size(dispatch.battout)),dispatch.battout,dispatch.solar,dispatch.gen,dispatch.crtllo,dispatch.crtlhi,solarE-dispatch.solar],...
    %    1,'stack');
    M = {'Battery Charging';'Battery Discharging';'Solar';'Diesel';'Nonserved Demand';'Solar Spill'};
    
    legend(M,'Location','SouthEastOutside');
    hold on
    area(-dispatch.battin);
    
    %bar(-dispatch.battin,1,'stack');
    plot(1:length(hindex),demand.HiPty+demand.LoPty,'k','LineWidth',5);
    hold off
%     plot(1:size(demand.HiPty,1),demand.HiPty+demand.LoPty+savedLosses+dispatch.battin);
%     plot(1:size(demand.HiPty,1),demand.HiPty+demand.LoPty+savedLosses);
%     plot(1:size(demand.HiPty,1),demand.HiPty+demand.LoPty);
%     plot(1:size(demand.HiPty,1),demand.HiPty);
    %plot(1:size(battSoC,1),((battSoC-BESS.SOCmin)./(BESS.SOCmax-BESS.SOCmin) - 1).*5,'k--o');
    figure(3);
    plot(1:size(battSOC,1),battSOC,'k');
    %pause(.2);
    
    
end
% 
% 
% 
% 
% 
% figure(2);
% area([-dispatchinit.battin,dispatchinit.battout,dispatchinit.solar,dispatchinit.gen,dispatchinit.crtllo,dispatchinit.crtlhi]);
% M = {'battin';'battout';'solar';'gen';'crtllo';'crtlhi'};
% 
% legend(M,'Location','SouthEastOutside');
% hold on
% %plot(1:24,demand.HiPty+demand.LoPty);
% plot(1:size(demand.HiPty,1),demand.HiPty+demand.LoPty+savedLossesInit);
% plot(1:size(demand.HiPty,1),demand.HiPty);
% hold off
belowPmin = dispatch.gen(find((dispatch.gen < genEMin))); 
belowPmin = belowPmin(find(belowPmin > 0));

battOps.battSOC = battSOC;
battOps.battEInMax = battEinMax;
battOps.battEOutMax = battEoutMax;
battOps.Q_total = Q_total;

end

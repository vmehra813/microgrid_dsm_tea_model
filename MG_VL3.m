function [dispatch ] = CC_VL( hindex, genCatalog, REMsettings, demand, Ppv, ES_size, gensetIndex_1, time_step)
%% Notes:
%1. Prioritize battery usage over diesel usage
%2. Simone's battery losses convetion: Pbatt out refers to the usable batt
%out (after losses removed), while Pbatt in refers to the total Power in
%(before losses)

%Logic: prioritize  PV and battery usage first

%If battery SOC is below the threshold and the generator is already on in
%the previous time step:
%Then the generator is turned on max (essentially, for ops reasons).

%If the generator is on at max for ops reasons: First, PV is subtracted from the load to determine the remaining power
%needed (balance)the battery then follows the load (tries to
%meet any load that the generator can't meet)
%If the generator is not turned on for ops reasons, then again, PV is subtracted from the total load. The battery
%tries to follow the remaining load. If it can't, then the generator will
%try to meet it.

%Next, the power into the battery is calculated: 
%If the battery is already discharging in the hour,then it cannot accept any power input. 
%Otherwise, power into the battery from DC sources is first calculated.
%Then, any remaining capacity is filled by excess AC power. 

%% Initializing Variables
BESS=genCatalog.BESS(1);
% Inverter=genCatalog.ACDCConverter;
% genset=genCatalog.genset;



%%Define Inputs
BESS=genCatalog.BESS(1); %the first battery
Inverter=genCatalog.ACDCConverter;
genset=genCatalog.genset;
SOC_set = 0.6; %[genCatalog.BESS.SOCmax];

%code not written for ES_size = 0
if ES_size == 0
    ES_size = 0.0001;
    SOC_set = 0;
end 

dt=time_step/60; %hourly?

% gen = genCatalog.genset; %%NEED TO CHANGE TO INCLUDE THE COST FUNCTION, LINEARIZED
originalDemand = demand; %save the original demand as a variable
expectedLosses = REMsettings.distLossFctr;
% expLoss(1) = expectedLosses;
solarPEEffNom = 1;%genCatalog.PV_PE.effNom; %0.95;


Ppv=[0; Ppv];

demand=[0; originalDemand(hindex)]*(1+expectedLosses);
nt= length(demand);%length(demand);
% for iDistrLoss = 1:2 %iterate to calculate the distribution loss factor
    
demand=[0; originalDemand(hindex)]*(1+expectedLosses);   
% expLoss(iDistrLoss) = expectedLosses;
l=demand;

Q=zeros(1,nt);
bil_ac=zeros(1,nt);
Pexcess_ac=zeros(1,nt);
Pexcess_dc=zeros(1,nt);
Pexcess=zeros(1,nt);
% P_lf=zeros(1,nt);
bilancioI=zeros(1,nt); %balance?

P_gen1=zeros(1,nt); %power of gen1?

Pout_bat=zeros(1,nt);
Pin_bat_max=zeros(1,nt);
Pout_bat_max=zeros(1,nt);
Pout_max_kbm=zeros(1,nt);
Pout_max_mdr=zeros(1,nt);
Pin_bat=zeros(1,nt);
Pin_bat_ac=zeros(1,nt);
Pin_bat_dc=zeros(1,nt);
Pin_max_kbm=zeros(1,nt);
Pin_max_mcr=zeros(1,nt);
P_bat=zeros(1,nt);
E1=zeros(1,nt);
E2=zeros(1,nt);

C_cc=zeros(1,nt);
E_cc=zeros(1,nt);
LL=zeros(1,nt);          % perdita percentuale di vita utile (per ciascun istante temporale)
Ah=zeros(1,nt);          % pacchetto di energia ceduto/ricevuto dalle batterie
Unmet_Load=zeros(1,nt);
Served_L=zeros(1,nt);
Pleft_ac=zeros(1,nt);
SOC=zeros(1,nt);
Ediss=zeros(1,nt);

Pdiss_AC=zeros(1,nt);
Pdiss_DC=zeros(1,nt);

%Vivian: Added these Below
PgenI_ops = zeros(1,nt);
% bilancio_temp = zeros(1,nt);
Pbatt = zeros(1,nt);
%% BATTERY

Cap_max=BESS.energy;                 % capacit? nominale della singola batteria [kWh]
Nbat=ES_size/Cap_max;

C_inv_bat=BESS.cost;
C_inv=C_inv_bat*Nbat;          % ?

Emax_bat=Cap_max*Nbat;         % capacit? massima dell'accumulo

SOC_min=BESS.SOCmin;                   % SOC stato di carica della batteria: SOC=Eacc(i)/Enom_bat
SOC_max=BESS.SOCmax;
c=BESS.c;                       % capacity ratio
k=BESS.k;                       % rate constant 1/h
alpha_c=BESS.alpha_c;                     % maximum charge rate A/Ah
Imax=BESS.Imax;                             % maximun charge current A
Imax_d= BESS.I_discharge_max;%BESS.Imax_d; ????????
Vnom=BESS.Vnom;                        % nominal voltage V

Pmax_mcc=Nbat*Imax*Vnom/1000;
Pmax_mdc=Nbat*Imax_d*Vnom/1000;
eta_rt=BESS.eff_batt_d_c^2;
eta_bat=BESS.eff_batt_d_c;

Qlife=BESS.LifeTimeThroughPut;                    % kWh
A_tot=Qlife*Nbat; 

% INVERTER
eta_inv=Inverter.invrtrEffNom;
eta_rect=Inverter.rctfrEffNom;


% Initial conditions
SOC(1)=BESS.SOCinit;
E1(1)=BESS.SOCinit*Emax_bat*c;
E2(1)=BESS.SOCinit*Emax_bat*(1-c);
E_bat(1)=E1(1)+E2(1);

C_cc(1)=0; %somma dei C_cc dalla prima ora dell'anno fino all'ultima ora del giorno precedente
E_cc(1)=0; %somma delle E_cc dalla prima ora dell'anno fino all'ultima ora del giorno precedente

%% GEN1

gen = genCatalog.genset; 
CF=REMsettings.dieselFCost;


Pnom1=gen.pMax(gensetIndex_1); %max of gen1

Pmin1= gen.pMin(gensetIndex_1); %original LF

%% GenCosts

pMax = Pnom1; % from above
pMin = Pmin1; % from above
genFuelCostPerKWHFun = @(genE) polyval(genset.LPerKWHCurve(gensetIndex_1,:),genE)*CF*genE;

%input genE into function genFuelCostPerKWHFun as pMin and pMax
genFuelCost_pMin = genFuelCostPerKWHFun(pMin);
genFuelCost_pMax = genFuelCostPerKWHFun(pMax);

m = (genFuelCost_pMax-genFuelCost_pMin)/(pMax - pMin);
y_int = genFuelCost_pMin - m*pMin;

if isnan(m) 
    m = 1;
    y_int = 1;
end 

%y = mx(-m*x1 +y1)
%%%%%%% y_int should be something like y_int=genFuelCost_pMin-m*pMin;
%m = operational costs (marginal costs)
%y-intercept = fixed costs


% Cfix_gen1= y_int;%gen.fix(gensetIndex_1);
% Cmar_fuel1=m; %gen.mar(gensetIndex_1);

%%

for i=2:nt

%% Define Battery Limits
    % limite sulla massima potenza di carica della batteria -- limit on the maximum power of the battery charge

    Pin_max_kbm(i)=abs((-k*c*Emax_bat+k*E1(i-1)*exp(-k*dt)+Q(i-1)*k*c*(1-exp(-k*dt)))/(1-exp(-k*dt)+c*(k*dt-1+exp(-k*dt))));
    Pin_max_mcr(i)=(1-exp(-alpha_c*dt))*(Emax_bat-Q(i-1))/dt;
    Pin_bat_max(i)=min([Pin_max_kbm(i)/eta_bat,Pin_max_mcr(i)/eta_bat,Pmax_mcc]);

    % limite sulla massima potenza di scarica della batteria -- limit on the maximum discharge power of the battery

    Pout_max_kbm(i)=abs((k*E1(i-1)*exp(-k*dt)+Q(i-1)*k*c*(1-exp(-k*dt)))/(1-exp(-k*dt)+c*(k*dt-1+exp(-k*dt))));
    Pout_max_mdr(i)=max((SOC(i-1)-SOC_min)*Emax_bat/dt,0);
    Pout_bat_max(i)=eta_bat*min([Pout_max_kbm(i),Pout_max_mdr(i),Pmax_mdc]);
    %Refers to the USABLE battout of the battery after losses
    % NOTE: Pout_bat(i) (vs. Pout_bat_max) refers to the USABLE battout (after battery losses) that needs to be
    % outputted

    %% If the gen is NOT on for operational purposes
%     if PgenI_ops(i) <= 0
        P_gen1(i) = 0; %Pgen-Pload <0 so assume generator is off first
        %calculate whether power is needed from the battery
        %assume load is DC
        Pbatt_bal(i) = Ppv(i)-l(i);
        
        if Pbatt_bal(i) > 0 %if there is excess energy
            Pbatt(i) = Pbatt_bal(i);
            Unmet_Load(i) = 0;
        elseif Pbatt_bal(i) < 0 %means that there is't excess energy from PV
            %check to see if generator is needed
            if -Pbatt_bal(i) >= Pout_bat_max(i)%*eta_inv
                Pbatt(i) = -Pout_bat_max(i); %Pout_bat(i) as the max possible (keep sign convention)
                Pgen_bal(i) = (Pbatt(i)-Ppv(i))+l(i); %should always be positive
                Unmet_Load(i) = Pgen_bal(i);
                Pdiss_AC(i) = 0;

            elseif -Pbatt_bal(i) < Pout_bat_max(i)%*eta_inv
                Pbatt(i) = Pbatt_bal(i);
                Unmet_Load(i) = 0;
                Pdiss_DC(i) = 0;
       
                
            end
        end 
%     end              
             
%% Calculate energy input into battery 
%Pout_bat and Pin_bat convention is both (+)

if Pbatt(i) >= 0 %energy going into battery
    Pin_bat(i) = min(Pbatt(i),Pin_bat_max(i)); %Pin_bat is always +
    Pin_bat_test = Pin_bat(i);
    Pbatt_test = Pbatt(i);
    Pin_bat_max_test = Pin_bat_max(i);
    Pdiss_DC(i) = Pbatt(i)-Pin_bat(i);
else %energy going out of battery 
    Pout_bat(i) = -Pbatt(i);
end 


%% Update battery model

        if Pin_bat(i)>0
            P_bat(i)=-Pin_bat(i)*eta_bat; %the actual amount going into the battery
        elseif Pout_bat(i)>0
            P_bat(i)=Pout_bat(i)/eta_bat; %the actual amount leaving the battery
        end         

         E1(i)=E1(i-1)*exp(-k*dt)+(Q(i-1)*k*c-P_bat(i))*(1-exp(-k*dt))/k-P_bat(i)*c*(k*dt-1+exp(-k*dt))/k;
         E2(i)=E2(i-1)*exp(-k*dt)+Q(i-1)*(1-c)*(1-exp(-k*dt))-P_bat(i)*(1-c)*(k*dt-1+exp(-k*dt))/k; 
         Q(i)=E1(i)+E2(i);
         SOC(i)=Q(i)/Emax_bat;
         

%% Error testing
        if Pin_bat(i) >0 && Pout_bat(i) > 0
             error('Battery charging and discharging in i');
        end
        
        if Pin_bat(i) < 0 || Pout_bat(i) < 0
            error('Power in or out in i is negative');
        end
        
%% testing
T = i;
if Pbatt(T) > 0 && P_gen1(i) > 0 %Pbatt is + (charging)
    eBal(T) = Pin_bat(T)-Pout_bat(T) - Ppv(T)+Pdiss_DC(T) + l(T); %(P_gen1(T)+Unmet_Load(T)-l(T)'-Pdiss_AC(T))*eta_rect;
end 

end %end of iterations through hours

 %iterations through grid losses
%% Save dispatch
dispatch.gen1=P_gen1(2:nt); %ok
dispatch.solar=Ppv(2:nt)'; %ok

dispatch.battout=Pout_bat(2:nt); %ok
dispatch.battin=Pin_bat(2:nt); %ok

dispatch.battE=Q(2:nt); %ok

dispatch.unmet=Unmet_Load(2:nt)'; %ok
dispatch.demand=originalDemand;%l(2:nt); %ok
dispatch.adjDemand = l(2:nt);
%%%INSERT HERE
dispatch.Pdiss_AC = Pdiss_AC(2:nt)'; %ok

dispatch.Pdiss_DC = Pdiss_DC(2:nt)'; %ok


index = find(dispatch.gen1 > 0); %ok
dispatch.genFlag(index) = 1; %ok
dispatch.genFlag(1) = 1; %ok
dispatch.genFlag = dispatch.genFlag'; %ok
dispatch.gen = dispatch.gen1'; %ok
dispatch.crtlhi = min(dispatch.unmet,originalDemand); %dispatch.unmet; %ok
dispatch.crtllo = zeros(length(dispatch.crtlhi),1); %ok

%%%edited

dispatch.battout = dispatch.battout'/eta_bat; %ok
dispatch.battin = dispatch.battin'*eta_bat; %turn to convention of Doug %ok

dispatch.solarAvail = Ppv; %ok
dispatch.solar = dispatch.solar'; %ok
dispatch.SOC = SOC'; %ok
dispatch.spill = (Pdiss_AC+Pdiss_DC)'; %ok

%now we are in doug's convention for battout and battin


% % % DON'T NEED THESE
dispatch.gen_batt = zeros(1,nt); %Pin_bat_ac(2:nt)';
dispatch.gen_spill = zeros(1,nt); %dispatch.Pdiss_AC;
dispatch.gen_direct = zeros(1,nt); %dispatch.gen - dispatch.gen_batt - dispatch.gen_spill;

dispatch.solar_batt = zeros(1,nt); %Pin_bat_dc(2:nt)'; %Pin_bat_dc still in Simone's convention so do not have losses removed
dispatch.solar_spill = zeros(1,nt); %dispatch.Pdiss_DC;
dispatch.solar_direct = zeros(1,nt); %dispatch.solar - dispatch.solar_batt-dispatch.solar_spill; %should show up as a little more than demand bc need to input a little more to make up for inverter losses

%
solarPEEffNom = genCatalog.PV_PE.effNom; %0.95;
% battPEOutEffNom = genCatalog.BESS_PE.outEffNom; %0.95;
% battPEInEffNom = genCatalog.BESS_PE.inEffNom; %0.95;
% solarPELoss = dispatch.solar.*(1-solarPEEffNom); %losses from solar - ok
% battout = dispatch.battout - dispatch.battin; %ok -- note: it is only one or the other at every hour
% battLossOut = dispatch.battout.*(1-eta_bat); %ok
% battLossIn = dispatch.battin.*((1-eta_bat)/eta_bat);%ok
% battLoss = battLossOut + battLossIn; %ok
% battPELoss = (dispatch.battout-battLossOut).*(1-battPEOutEffNom) +...
%     (dispatch.battin+battLossIn).*((1-battPEInEffNom)/battPEInEffNom); %ok
% 
% % 
% % netDCOut = (dispatch.solar-dispatch.Pdiss_DC)-solarPELoss+battout-battLoss-battPELoss; %battout includes inv losses
% netDCOut = (dispatch.solar-dispatch.Pdiss_DC)-solarPELoss+battout-battLoss-battPELoss; %battout includes inv losses
% 
% 
% dispatch.invLoss = zeros(nt-1,1);%ok
% dispatch.rectLoss = zeros(nt-1,1); %ok
% 
% %if totalDCOut > 0 %E going out of batt (must pass thro inverter)
% dispatch.invLoss(netDCOut > 0) = netDCOut(netDCOut > 0)/eta_inv*(1-eta_inv); 
% 
% 
% %if totalDCOut < 0 %E going into batt (must pass thro rectifier)
% dispatch.rectLoss(netDCOut < 0) = -netDCOut(netDCOut < 0)/eta_rect*(1-eta_rect); %ok

% 
% dispatch.gen_batt = zeros(nt-1,1); %ok
% % netAC= dispatch.gen-dispatch.demand-dispatch.Pdiss_AC-dispatch.unmet;
% netAC= dispatch.gen-dispatch.demand-dispatch.Pdiss_AC+dispatch.unmet; %Unmet should be subratracted from demand
% 
% 
% dispatch.gen_batt(netAC>0) = netAC(netAC>0); %ok
% dispatch.gen_direct = dispatch.gen-dispatch.gen_batt-dispatch.Pdiss_AC; %ok
% dispatch.gen_spill = dispatch.Pdiss_AC; %ok
% 
% 
% dispatch.solar_batt = zeros(nt-1,1);
% 
% 
% for i=1:nt-1
%     if battout(i)<0 %where battin occurs
%         dispatch.solar_batt(i)= Pin_bat_dc(i);%-battout(i)/eta_bat-dispatch.gen_batt(i);
% 
%     end
% end
% %dispatch.solar(netAC<=0)    - dispatch.Pdiss_DC(netAC<=0)
% dispatch.solar_direct = dispatch.solar - dispatch.solar_batt-dispatch.Pdiss_DC;
% dispatch.solar_spill = dispatch.Pdiss_DC;

%losses - ok
% dispatch.nongridLosses = battLoss;%-dispatch.battout/eta_bat*(1-eta_bat)+dispatch.battin*(1-eta_bat);%losses from running th battery

% demandMet = sum(demand) - sum(dispatch.crtllo) - sum(dispatch.crtlhi); %how much power was actually generated and sent through wires
% distLosses = demandMet*REMsettings.distLossFctr; %the losses from the power generated using expected dist loss factor
% glosses = demandMet - demandMet/(1+REMsettings.distLossFctr);
% dispatch.expectedLosses = glosses/sum(originalDemand) ; %update the losses for the next iteration (add expected losses onto demand)
% 
% 
% dispatch.distLosses = glosses; %re-calculate losses from sending power through wires
% 

% end  %end of iterations through grid losses

end
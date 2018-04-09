function [reliability_critical, critical_amount_served,reliability_total, total_amount_served, cycles] = connect_ulink_customer_gen_sizing_kibam_load_dispatch_8760x1(available_solar,battery,loadprof,num,batt,load)

%% Given remaining SOC from A load, ability to serve another customer.

%right now code can take in multiple a generators but is only using 1 a box
%remaining SOC
%incorporate reliability preference
addpath('/Users/vmehra813/Dropbox (MIT)/[uLink]/Simulation/MatLab/uLink/Modeling Analysis/Generation + Reliability Model/rainflow');
addpath('/Users/ramatya/Dropbox (MIT)/[uLink]/Simulation/MatLab/uLink/Modeling Analysis/Generation + Reliability Model/rainflow');


%% Battery Inputs

if num.consumers <= 10
    num_batt = 1;
else
    num_batt = batt_number_final;
    
end

% deltaT = 1;   % Time step? 1 hour?
% alpha_c = 1;  % max charge rate in A/Ah i.e same as saying 1 Ah capacity battery charging at 1 hour rate
% Imax = 16.5;
deltaT = battery.deltaT;
alpha_c = battery.alpha_c;
Imax = battery.Imax;

%battery.voltage = 12;

%define these for HBL battery! have inputs in main script under battery
%data structure
c = battery.c;
k = battery.k;

% c= 0.305; % capacity ratio
% k = 2.12; % rate constant

%roundtrip efficiency sqrt is what HOMER uses for charge/discharge eff
%Qmax_batt = batt/1;                 % currently have this in Ah (maximum battery capacity)
t=battery.t2;
Qmax_batt = (batt*(1-exp(-t*k)*(1-c)+t*c*k))/(t*c*k);
Qmax_battbank = Qmax_batt*num_batt; % Qmax_battbank is max capacity of the battery bank in Ah


e_raised = exp(-k*deltaT);


%Initialize Matrices
total.Q_1 = zeros(length(loadprof.final.total),1); %some have an added a row bc of t+1 iteration later on
total.Q_2 = zeros(length(loadprof.final.total),1);
%battEcum = zeros(length(loadprof.final.total),1);
total.battEoutMax = zeros(length(loadprof.final.total),1);
total.battEinMax = zeros(length(loadprof.final.total),1);
total.Ebatt_max_in_lim1 = zeros(length(loadprof.final.total),1);
total.Ebatt_max_in_lim2 = zeros(length(loadprof.final.total),1);
total.Ebatt_max_in_lim3 = zeros(length(loadprof.final.total),1);
total.Ebatt_max_in_limSOCmax = zeros(length(loadprof.final.total),1);
total.EoutMax_kineticmodellimit = zeros(length(loadprof.final.total),1);
total.battSOC = zeros(length(loadprof.final.total),1);
total.Q_total = zeros(length(loadprof.final.total),1);


critical.Q_1 = zeros(length(loadprof.final.total),1); %some have an added a row bc of t+1 iteration later on
critical.Q_2 = zeros(length(loadprof.final.total),1);
%battEcum = zeros(length(loadprof.final.total),1);
critical.battEoutMax = zeros(length(loadprof.final.total),1);
critical.battEinMax = zeros(length(loadprof.final.total),1);
critical.Ebatt_max_in_lim1 = zeros(length(loadprof.final.total),1);
critical.Ebatt_max_in_lim2 = zeros(length(loadprof.final.total),1);
critical.Ebatt_max_in_lim3 = zeros(length(loadprof.final.total),1);
critical.Ebatt_max_in_limSOCmax = zeros(length(loadprof.final.total),1);
critical.EoutMax_kineticmodellimit = zeros(length(loadprof.final.total),1);
critical.battSOC = zeros(length(loadprof.final.total),1);
critical.Q_total = zeros(length(loadprof.final.total),1);


battery.initial_capacity = .6;

%% Defining Battery Charge / Discharge [TOTAL vs NON CRITICAL]
%these are matrices
%avail_solar = repmat(available_solar,1,num.consumers);
avail_solar = available_solar;



loadprof.final.non_critical = loadprof.final.total;

%power.in_out_batt_total = avail_solar - loadprof.final.total;
power.in_out_batt_total =  loadprof.final.total(:,end) - avail_solar; %negative is charging
power.in_out_batt_critical =  loadprof.final.critical(:,end) - avail_solar;
current.in_out_batt_total = power.in_out_batt_total / battery.voltage;
current.in_out_batt_critical = power.in_out_batt_critical / battery.voltage;


%% Initialize Cycle Counter 

cycles.num = 365; %assume max 2 charge/discharge cycles per day 356*2
cycles.total_mean = [];
cycles.critical_mean = [];
cycles.total_use = [];


%% Initial Condition
%first hour is the same
%total
total.Q_1(1) = c*battery.initial_capacity*Qmax_battbank; %these are the starting values of Q_1 and Q_2 for the first iteration - this is a property of the battery bank
total.Q_2(1) = (1-c)*battery.initial_capacity*Qmax_battbank;

%battEcum(1) = 0;
total.Q_total(1) = total.Q_1(1) + total.Q_2(1);
total.battSOC(1) = total.Q_total(1)/(Qmax_battbank);

%removed -1* 
% total.EoutMax_kineticmodellimit(1) = -1*battery.discharge_eff*deltaT*((-k*c*Qmax_battbank + k*total.Q_1(1)*e_raised + total.Q_total(1)*k*c*(1-e_raised))/...
%     (1-e_raised+c*(k*deltaT-1+e_raised)));

%EoutMax_SOCminlimit = (total.battSOC(1)-battery.MDOD)*Qmax_battbank;
% total.EoutMax_SOCminlimit(1) = (total.battSOC(1)-battery.MDOD)*Qmax_battbank;
%
% total.battEoutMax(1)= min(total.EoutMax_kineticmodellimit(1),total.EoutMax_SOCminlimit(1));
%
%
% total.Ebatt_max_in_lim1(1) = deltaT*k*total.Q_1(1)*e_raised+total.Q_total(1)*k*c*(1-e_raised)/(1-e_raised+c*(k*deltaT-1+e_raised));
% total.Ebatt_max_in_lim2(1) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-total.Q_total(1))/deltaT;
% total.Ebatt_max_in_lim3(1) = deltaT*num_batt*Imax/1;
% %Ebatt_max_in_limSOCmax = (BESS.SOCmax-total.battSOC(1))*Qmax_battbank;
% %Ebatt_max_in_limSOCmax = (1-total.battSOC(1))*Qmax_battbank;
% total.Ebatt_max_in_limSOCmax(1) = (1-total.battSOC(1))*Qmax_battbank;
%
% % energy into the battery during charge (Wh)
% total.battEinMax(1) = min([total.Ebatt_max_in_lim1(1), total.Ebatt_max_in_lim2(1), total.Ebatt_max_in_lim3(1), total.Ebatt_max_in_limSOCmax(1)])/battery.charge_eff;
%

%Critical
critical.Q_1(1) = c*battery.initial_capacity*Qmax_battbank; %these are the starting values of Q_1 and Q_2 for the first iteration - this is a property of the battery bank
critical.Q_2(1) = (1-c)*battery.initial_capacity*Qmax_battbank;

%battEcum(1) = 0;
critical.Q_total(1) = critical.Q_1(1) + critical.Q_2(1);
critical.battSOC(1) = critical.Q_total(1)/(Qmax_battbank);

%removed -1*
% critical.EoutMax_kineticmodellimit(1) = -1*battery.discharge_eff*deltaT*((-k*c*Qmax_battbank + k*critical.Q_1(1)*e_raised + critical.Q_total(1)*k*c*(1-e_raised))/...
%     (1-e_raised+c*(k*deltaT-1+e_raised)));
%
% %EoutMax_SOCminlimit = (critical.battSOC(1)-battery.MDOD)*Qmax_battbank;
% critical.EoutMax_SOCminlimit(1) = (critical.battSOC(1)-battery.MDOD)*Qmax_battbank;
%
% critical.battEoutMax(1)= max(critical.EoutMax_kineticmodellimit(1),critical.EoutMax_SOCminlimit(1));
%
%
% critical.Ebatt_max_in_lim1(1) = deltaT*k*critical.Q_1(1)*e_raised+critical.Q_total(1)*k*c*(1-e_raised)/(1-e_raised+c*(k*deltaT-1+e_raised));
% critical.Ebatt_max_in_lim2(1) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-critical.Q_total(1))/deltaT;
% critical.Ebatt_max_in_lim3(1) = deltaT*num_batt*Imax/1;
% %Ebatt_max_in_limSOCmax = (BESS.SOCmax-critical.battSOC(1))*Qmax_battbank;
% %Ebatt_max_in_limSOCmax = (1-critical.battSOC(1))*Qmax_battbank;
% critical.Ebatt_max_in_limSOCmax(1) = (1-critical.battSOC(1))*Qmax_battbank;
%
% % energy into the battery during charge (Wh)
% critical.battEinMax(1) = min([critical.Ebatt_max_in_lim1(1), critical.Ebatt_max_in_lim2(1), critical.Ebatt_max_in_lim3(1), critical.Ebatt_max_in_limSOCmax(1)])/battery.charge_eff;
%


%total
%initialize vectors
network.total.reliability = zeros(1,1);
network.total.amountserved.cumulative = zeros(length(loadprof.final.total),1);
network.total.nonserved = zeros(length(loadprof.final.total),1);
%network.total.nonserved = [];
%network.total.amountserved.cumulative(1) = 20;
network.total.solar_load_coincidence = zeros(length(loadprof.final.total),1);


SOC.total = zeros(length(loadprof.final.total), 1);
SOC.total_spillage = zeros(length(loadprof.final.total),1);
SOC.total_cycles = zeros(1,1);

%Critical
%initialize vectors
network.critical.reliability = zeros(1,1);
network.critical.amountserved.cumulative = zeros(length(loadprof.final.total),1);
network.critical.nonserved = zeros(length(loadprof.final.total),1);
%network.critical.nonserved = [];
%network.total.amountserved.cumulative(1) = 20;
network.critical.solar_load_coincidence = zeros(length(loadprof.final.total),1);


SOC.critical = zeros(length(loadprof.final.total), 1);
SOC.critical_spillage = zeros(length(loadprof.final.total),1);
SOC.critical_cycles = zeros(1,1);

%count variables
% zero = current; 1 = KIBAM; 2 = SOC
%SOC.charge_total_limit = zeros(length(loadprof.final.total),1);
%SOC.charge_critical_limit = zeros(length(loadprof.final.total),1);
%SOC.discharge_total_limit = zeros(length(loadprof.final.total),1);
%SOC.discharge_critical_limit = zeros(length(loadprof.final.total),1);
SOC.charge_total_limit = [];
SOC.charge_critical_limit = [];
SOC.discharge_total_limit = [];
SOC.discharge_critical_limit = [];



%% Create load stack

total.load_stack = [load.Critical, load.Power/battery.voltage, num.consumers*load.AppNum];
total.load_stack_sort = sortrows(total.load_stack,[1 -3 2]); %sorts by critical
total.non_critical_loads = find(total.load_stack_sort(:,1)<1);
total.non_critical_load_ind = total.non_critical_loads(end);

% 1st column is criticality, 2nd column is current draw, 3rd column is
% number of appliances over network


total.count_phone = zeros(length(loadprof.final.total),1);
total.count_fan = zeros(length(loadprof.final.total),1);
total.count_light = zeros(length(loadprof.final.total),1);
total.curtail = zeros(length(loadprof.final.total),1);

% figure;
% plot(total.load_stack_sort(1:total.non_critical_load_ind,2).*total.load_stack_sort(1:total.non_critical_load_ind,3)*battery.voltage);
% xlabel('Load / Appliance Number');
% ylabel('Total Load, Watts');
% title('Load Stack for N=5 for Non-Critical Loads: Power * Num Appliances');


%% Total (KIBAM) (now Total Load)
%clear j i;
%j=1;
%while j<= num.consumers
%for j=1:num.consumers %column is B customer    
    for i=1
        
        if current.in_out_batt_total(i) < 0     %charging
            
            % need other 1st condition to include solar/load for hour 1
            % total.Q_1(i) = critical.Q_1(i)*e_raised+(critical.Q_total(i)*k*c+current.in_out_batt_total(i))*(1-e_raised)/k+power.in_out_batt_total(i)*c*(k*deltaT-1+e_raised)/k;
            % critical.Q_2(i) = critical.Q_2(i)*e_raised+(critical.Q_total(i)*(1-c))*(1-e_raised)+current.in_out_batt_total(i)*(1-c)*(k*deltaT-1+e_raised)/k;
            % critical.Q_total(i) = critical.Q_1(i) + critical.Q_2(i);
            % critical.battSOC(i) = critical.Q_total(i)/Qmax_battbank;
            
            total.Ebatt_max_in_lim1(i) = (deltaT*((-k*c*Qmax_battbank + k*total.Q_1(i)*e_raised + total.Q_total(i)*k*c*(1-e_raised))/...
                (1-e_raised+c*(k*deltaT-1+e_raised))));
            %total.Ebatt_max_in_lim2(i) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-total.Q_total(i))/deltaT;
            %total.Ebatt_max_in_lim3(i) = deltaT*num_batt*Imax/1;
            total.Ebatt_max_in_limSOCmax(i) = -1*(battery.MLOC-total.battSOC(i))*Qmax_battbank;
            %energy into the battery during charge (Wh)
            %total.battEinMax(i) = min([total.Ebatt_max_in_lim1(i), total.Ebatt_max_in_lim2(i), total.Ebatt_max_in_lim3(i), total.Ebatt_max_in_limSOCmax(i)])/battery.charge_eff;
            total.battEinMax(i) = max([total.Ebatt_max_in_lim1(i),total.Ebatt_max_in_limSOCmax(i),current.in_out_batt_total(i)]);
            
            
            % total spillage (Wh)
            SOC.total_spillage(i) = abs(avail_solar(i)-loadprof.final.total(i,end)-battery.voltage*total.Q_1(i))*ge(avail_solar(i),loadprof.final.total(i,end)+battery.voltage*total.Q_1(i));
            
            
        elseif current.in_out_batt_total(i) > 0 %discharging
            
            % need other 1st condition to include solar/load for hour 1
            
            % critical.Q_1(i) = critical.Q_1(i)*e_raised+(critical.Q_total(i)*k*c+current.in_out_batt_total(i))*(1-e_raised)/k+power.in_out_batt_total(i)*c*(k*deltaT-1+e_raised)/k;
            % critical.Q_2(i) = critical.Q_2(i)*e_raised+(critical.Q_total(i)*(1-c))*(1-e_raised)+current.in_out_batt_total(i)*(1-c)*(k*deltaT-1+e_raised)/k;
            % critical.Q_total(i) = critical.Q_1(i) + critical.Q_2(i);
            % critical.battSOC(i) = critical.Q_total(i)/Qmax_battbank;
            
            total.EoutMax_kineticmodellimit(i) = deltaT*(k*total.Q_1(i)*e_raised+total.Q_total(i)*k*c*(1-e_raised))/(1-e_raised+c*(k*deltaT-1+e_raised));
            
            total.EoutMax_SOCminlimit(i) = (total.battSOC(i)-(1-battery.MDOD))*Qmax_battbank;
            % energy out from the battery during discharge (Wh)
            total.battEoutMax(i)= min([total.EoutMax_kineticmodellimit(i),total.EoutMax_SOCminlimit(i),current.in_out_batt_total(i)]);
            
            % total non-served energy (Wh)
            network.total.nonserved(i) = (loadprof.final.total(i,end) - battery.voltage*total.Q_1(i)-avail_solar(i))*ge(loadprof.final.total(i,end),(avail_solar(i)+battery.voltage*total.Q_1(i)));
 
        end
        
    end
    
    for i=2:length(loadprof.final.total)
        if current.in_out_batt_total(i) < 0       % charging
            
            %limits
            
            total.Ebatt_max_in_lim1(i) = (deltaT*((-k*c*Qmax_battbank + k*total.Q_1(i-1)*e_raised + total.Q_total(i-1)*k*c*(1-e_raised))/...
                (1-e_raised+c*(k*deltaT-1+e_raised))));
            %total.Ebatt_max_in_lim2(i) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-total.Q_total(i-1))/deltaT;
            %total.Ebatt_max_in_lim3(i) = deltaT*num_batt*Imax*battery.voltage/1;
            total.Ebatt_max_in_limSOCmax(i) = -1*(battery.MLOC-total.battSOC(i-1))*Qmax_battbank;
            %total.battEinMax(i) = min([total.Ebatt_max_in_lim1(i), total.Ebatt_max_in_lim2(i), total.Ebatt_max_in_lim3(i), total.Ebatt_max_in_limSOCmax(i)])/battery.charge_eff;
            %total.battEinMax(i) = min([total.Ebatt_max_in_lim1(i), total.Ebatt_max_in_lim3(i), total.Ebatt_max_in_limSOCmax(i)])/battery.charge_eff;
            %changed to min
            total.battEinMax(i) = max([total.Ebatt_max_in_lim1(i), total.Ebatt_max_in_limSOCmax(i),current.in_out_batt_total(i)])*battery.charge_eff;
            
            if total.battEinMax(i) == current.in_out_batt_total(i)
                SOC.charge_total_limit(i) = 1;
                
            elseif total.battEinMax(i) == total.Ebatt_max_in_lim1(i)
                SOC.charge_total_limit(i) = 2;
            elseif total.battEinMax(i) == total.Ebatt_max_in_limSOCmax(i)
                SOC.charge_total_limit(i) = 3;
                
            end
            
            
            total.Q_1(i) = total.Q_1(i-1)*e_raised+(total.Q_total(i-1)*k*c-total.battEinMax(i))*(1-e_raised)/k-total.battEinMax(i)*c*(k*deltaT-1+e_raised)/k;
            total.Q_2(i) = total.Q_2(i-1)*e_raised+(total.Q_total(i-1)*(1-c))*(1-e_raised)-total.battEinMax(i)*(1-c)*(k*deltaT-1+e_raised)/k;
            total.Q_total(i) = total.Q_1(i) + total.Q_2(i);
            total.battSOC(i) = abs(total.Q_total(i)/Qmax_battbank);      % + battSOC(i-1);
            
            
            if total.battSOC(i)+total.battSOC(i-1) < battery.MLOC
                total.battSOC(i) = total.battSOC(i)+total.battSOC(i-1);
            else
                %total.battSOC(i) = total.battSOC(i-1);
                total.battSOC(i) = battery.MLOC;
            end
            
            %spillage: if total.Q_1 is larger than batt_e_in_max
            %set to previous condition
            %limits
            
            SOC.total_spillage(i) = abs(avail_solar(i)-loadprof.final.total(i,end)-battery.voltage*total.Q_1(i))*ge(avail_solar(i),loadprof.final.total(i,end)+battery.voltage*total.Q_1(i));
            
        elseif current.in_out_batt_total(i) == 0
            
            %critical.Q_1(i) = critical.Q_1(i-1);
            %critical.Q_2(i) = critical.Q_2(i-1);
            %critical.Q_total(i) = critical.Q_total(i-1);
            total.battSOC(i) = total.battSOC(i-1);
            
        elseif current.in_out_batt_total(i) > 0       % discharging
            
            %limits
            total.EoutMax_kineticmodellimit(i) = deltaT*(k*total.Q_1(i-1)*e_raised+total.Q_total(i-1)*k*c*(1-e_raised))/(1-e_raised+c*(k*deltaT-1+e_raised));
            total.EoutMax_SOCminlimit(i) = (total.battSOC(i-1)-(1-battery.MDOD))*Qmax_battbank;
            total.battEoutMax(i)= min([total.EoutMax_kineticmodellimit(i),total.EoutMax_SOCminlimit(i),current.in_out_batt_total(i)])*battery.discharge_eff;
            
            
            if total.battEoutMax(i) == current.in_out_batt_total(i)
                SOC.discharge_total_limit(i) = 1;
                
            elseif total.battEoutMax(i) == total.EoutMax_kineticmodellimit(i)
                SOC.discharge_total_limit(i) = 2;
            elseif total.battEoutMax(i) == total.EoutMax_SOCminlimit(i)
                SOC.discharge_total_limit(i) = 3;
                
            end
            

            total.Q_1(i) = total.Q_1(i-1)*e_raised+(total.Q_total(i-1)*k*c - total.battEoutMax(i))*(1-e_raised)/k - total.battEoutMax(i)*c*(k*deltaT-1+e_raised)/k;
            total.Q_2(i) = total.Q_2(i-1)*e_raised+(total.Q_total(i-1)*(1-c))*(1-e_raised) - total.battEoutMax(i)*(1-c)*(k*deltaT-1+e_raised)/k;
            total.Q_total(i) = total.Q_1(i) + total.Q_2(i);
            total.battSOC(i) = abs(total.Q_total(i)/Qmax_battbank);          % + battSOC(i-1);
            
            if total.battSOC(i-1) - total.battSOC(i) > (1- battery.MDOD)
                total.battSOC(i) = total.battSOC(i-1) - total.battSOC(i);
                
            else
                %total.battSOC(i) = total.battSOC(i-1);
                total.battSOC(i) = 1 - battery.MDOD;
            end
            
            
            % non-served energy -- LOAD STACK.
            
%             if total.battEoutMax(i) >= abs(current.in_out_batt_total(i)) %if there is non-served energy
%                 continue;
%             else  %if <=
%                 
%                 total.curtail(i) = current.in_out_batt_total(i) - total.battEoutMax(i);
%                 
%                 while total.curtail(i) >= 0
%                     
%                     if loadprof.cumulative.non_critical_phone(i) ~=0 %phone load in hour
%                         
%                         num_phone = min([loadprof.cumulative.non_critical_phone(i)/(battery.voltage*total.load_stack_sort(1,2)),total.load_stack_sort(1,3)]);
%                         
%                         for o=1:num_phone
%                             if total.curtail(i) > o*loadprof.cumulative.non_critical_phone(i)/battery.voltage
%                                 total.count_phone(i) = o;
%                                 total.curtail(i) = total.curtail(i) - o*loadprof.cumulative.non_critical_phone(i)/battery.voltage;
%                             else
%                                 break;
%                             end
%                         end
%                     end
%                     
%                     if loadprof.cumulative.non_critical_fan(i) ~=0 %fan load in hour
%                         
%                         num_fan = min([loadprof.cumulative.non_critical_fan(i)/(battery.voltage*total.load_stack_sort(2,2)),total.load_stack_sort(2,3)]);
%                         
%                         for o=1:num_fan
%                             if total.curtail(i) > o*loadprof.cumulative.non_critical_fan(i)/battery.voltage
%                                 total.count_fan(i) = o;
%                                 total.curtail(i) = total.curtail(i) - o*loadprof.cumulative.non_critical_fan(i)/battery.voltage;
%                             else
%                                 break;
%                             end
%                         end
%                     end
%                     
%                     break;
%                     
%                     %need to exit while loop
%        
%                 end
                
           % end
            
            %in current (A)
            %battEoutMax is negative
           % if total.count_fan(i)==0 && total.count_phone(i)==0
                
            %    network.total.nonserved(i) = (loadprof.final.total(i) - battery.voltage*total.battEoutMax(i)-avail_solar(i))*ge(loadprof.final.total(i),(avail_solar(i)+battery.voltage*total.Q_1(i)));
                
            %else total.count_fan(i)~=0 || total.count_phone(i)~=0
           % else   
           %     network.total.nonserved(i) = total.count_fan(i)*battery.voltage*total.load_stack_sort(2,2) - total.count_phone(i)*battery.voltage*total.load_stack_sort(1,2);
                %may end up leading to negative values for nonserved energy
                
           % end
            %             nonserved = loadprof.final.total(i) + total.battEoutMax(i) - avail_solar(i);
            %             if nonserved > 0
            %                 network.total.nonserved(i) = nonserved;
            %             end
            
            %if network.total.nonserved(i)<0
           %     network.total.nonserved(i)=0;
           % end
            
        end
        
        
        network.total.nonserved(i) = (loadprof.final.total(i,end) - battery.voltage*total.Q_1(i)-avail_solar(i))*ge(loadprof.final.total(i,end),(avail_solar(i)+battery.voltage*total.Q_1(i)));

        
    end
      
network.total.reliability = 1- (sum(network.total.nonserved)/sum(loadprof.final.total(:,end)));  
    %j=j+1;
    


% % FIGURE

% figure;
% rfmatrix(cycles.total_use);
% set(gca,'fontsize',15);
% xlabel('SOC Amplitude');
% ylabel('SOC Mean Value');
% zlabel('Number of Cycles');
% title(' ');
% 
% figure;
% plot(cycles.total_reshape_2);
% set(gca,'fontsize',15);
% xlabel('Hour of Day');
% ylabel('State of Charge (Percent)');

% figure;
% plot(diff(total.battSOC), 'c');
% hold on;
% plot(total.battSOC);
% set(gca,'fontsize',15);
% xlabel('Hour of Year');
% ylabel('State of Charge (Percent)');
% legend('Differential SOC','Total SOC');



% CALCULATE total RELIABILITY
% CALCULATE total AMOUNT SERVED
reliability_total = network.total.reliability;
total_amount_served = sum(loadprof.final.total(:,end)) - sum(network.total.nonserved);


%% Critical (KIBAM)
%clear j i;
%j=1;
%while j<= num.consumers
%for j=1:num.consumers %column is B customer
    
    for i=1
        
        if current.in_out_batt_critical(i) < 0     %charging
            
            % need other 1st condition to include solar/load for hour 1
            % critical.Q_1(i) = critical.Q_1(i)*e_raised+(critical.Q_total(i)*k*c+current.in_out_batt_critical(i))*(1-e_raised)/k+power.in_out_batt_critical(i)*c*(k*deltaT-1+e_raised)/k;
            % critical.Q_2(i) = critical.Q_2(i)*e_raised+(critical.Q_total(i)*(1-c))*(1-e_raised)+current.in_out_batt_critical(i)*(1-c)*(k*deltaT-1+e_raised)/k;
            % critical.Q_total(i) = critical.Q_1(i) + critical.Q_2(i);
            % critical.battSOC(i) = critical.Q_total(i)/Qmax_battbank;
            
            critical.Ebatt_max_in_lim1(i) = (deltaT*((-k*c*Qmax_battbank + k*critical.Q_1(i)*e_raised + critical.Q_total(i)*k*c*(1-e_raised))/...
                (1-e_raised+c*(k*deltaT-1+e_raised))));
            
            %critical.Ebatt_max_in_lim2(i) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-critical.Q_total(i))/deltaT;
            %critical.Ebatt_max_in_lim3(i) = deltaT*num_batt*Imax/1;
            critical.Ebatt_max_in_limSOCmax(i) = -1*(battery.MLOC-critical.battSOC(i))*Qmax_battbank;
            % energy into the battery during charge (Wh)
            %critical.battEinMax(i) = min([critical.Ebatt_max_in_lim1(i), critical.Ebatt_max_in_lim2(i), critical.Ebatt_max_in_lim3(i), critical.Ebatt_max_in_limSOCmax(i)])/battery.charge_eff;
            critical.battEinMax(i) = max([critical.Ebatt_max_in_lim1(i), critical.Ebatt_max_in_limSOCmax(i), current.in_out_batt_critical(i)]);
            
            
            % critical spillage (Wh)
            %SOC.critical_spillage(i) = abs(avail_solar(i)-loadprof.final.total(i)+battery.voltage*critical.battEinMax(i))*ge(avail_solar(i),loadprof.final.critical(i)-battery.voltage*critical.battEinMax(i));
            SOC.critical_spillage(i) = abs(avail_solar(i)-loadprof.final.critical(i,end)-battery.voltage*critical.Q_1(i))*ge(avail_solar(i),loadprof.final.critical(i,end)+battery.voltage*critical.Q_1(i));
            
            
        elseif current.in_out_batt_critical(i) > 0 %discharging
            
            % need other 1st condition to include solar/load for hour 1
            
            % critical.Q_1(i) = critical.Q_1(i)*e_raised+(critical.Q_total(i)*k*c+current.in_out_batt_critical(i))*(1-e_raised)/k+power.in_out_batt_critical(i)*c*(k*deltaT-1+e_raised)/k;
            % critical.Q_2(i) = critical.Q_2(i)*e_raised+(critical.Q_total(i)*(1-c))*(1-e_raised)+current.in_out_batt_critical(i)*(1-c)*(k*deltaT-1+e_raised)/k;
            % critical.Q_total(i) = critical.Q_1(i) + critical.Q_2(i);
            % critical.battSOC(i) = critical.Q_total(i)/Qmax_battbank;
            
            
            critical.EoutMax_kineticmodellimit(i) = deltaT*(k*total.Q_1(i)*e_raised+total.Q_total(i)*k*c*(1-e_raised))/(1-e_raised+c*(k*deltaT-1+e_raised));
            
            critical.EoutMax_SOCminlimit(i) = (critical.battSOC(i)-(1-battery.MDOD))*Qmax_battbank;
            % energy out from the battery during discharge (Wh)
            critical.battEoutMax(i)= min([critical.EoutMax_kineticmodellimit(i),critical.EoutMax_SOCminlimit(i),current.in_out_batt_critical(i)]);
            
            % critical non-served energy (Wh)
            %network.critical.nonserved(i) = (loadprof.final.critical(i) - battery.voltage*critical.battEoutMax(i)-avail_solar(i))*ge(loadprof.final.critical(i),(avail_solar(i)+battery.voltage*critical.battEoutMax(i)));
            network.critical.nonserved(i) = (loadprof.final.critical(i,end) - battery.voltage*critical.Q_1(i)-avail_solar(i))*ge(loadprof.final.critical(i,end),(avail_solar(i)+battery.voltage*critical.Q_1(i)));
      
        end
        
    end
    
    for i=2:length(loadprof.final.total)
        if current.in_out_batt_critical(i) < 0       % charging
            
            critical.Ebatt_max_in_lim1(i) = (deltaT*((-k*c*Qmax_battbank + k*critical.Q_1(i-1)*e_raised + critical.Q_total(i-1)*k*c*(1-e_raised))/...
                (1-e_raised+c*(k*deltaT-1+e_raised))));
            %critical.Ebatt_max_in_lim2(i) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-critical.Q_total(i))/deltaT;
            %critical.Ebatt_max_in_lim3(i) = deltaT*num_batt*Imax*battery.voltage/1;
            critical.Ebatt_max_in_limSOCmax(i) = -1*(battery.MLOC-critical.battSOC(i-1))*Qmax_battbank;
            %critical.battEinMax(i) = min([critical.Ebatt_max_in_lim1(i), critical.Ebatt_max_in_lim2(i), critical.Ebatt_max_in_lim3(i), critical.Ebatt_max_in_limSOCmax(i)])/battery.charge_eff;
            critical.battEinMax(i) = max([critical.Ebatt_max_in_lim1(i), critical.Ebatt_max_in_limSOCmax(i), current.in_out_batt_critical(i)])*battery.charge_eff;
            

            
            critical.Q_1(i) = critical.Q_1(i-1)*e_raised+(critical.Q_total(i-1)*k*c-critical.battEinMax(i))*(1-e_raised)/k-critical.battEinMax(i)*c*(k*deltaT-1+e_raised)/k;
            critical.Q_2(i) = critical.Q_2(i-1)*e_raised+(critical.Q_total(i-1)*(1-c))*(1-e_raised)-critical.battEinMax(i)*(1-c)*(k*deltaT-1+e_raised)/k;
            critical.Q_total(i) = critical.Q_1(i) + critical.Q_2(i);
            critical.battSOC(i) = abs(critical.Q_total(i)/Qmax_battbank);      % + battSOC(i-1);
            
            
            if critical.battSOC(i)+critical.battSOC(i-1) < battery.MLOC
                critical.battSOC(i) = critical.battSOC(i)+critical.battSOC(i-1);
            else
                %critical.battSOC(i) = critical.battSOC(i-1);
                critical.battSOC(i) = battery.MLOC;
            end
            
            %spillage: if critical.Q_1 is larger than batt_e_in_max
            %set to previous condition
            %limits
            
            SOC.critical_spillage(i) = abs(avail_solar(i)-loadprof.final.critical(i,end)-battery.voltage*critical.Q_1(i))*ge(avail_solar(i),loadprof.final.critical(i,end)+battery.voltage*critical.Q_1(i));
            
        elseif current.in_out_batt_critical(i) == 0
            
            %critical.Q_1(i) = critical.Q_1(i-1);
            %critical.Q_2(i) = critical.Q_2(i-1);
            %critical.Q_total(i) = critical.Q_total(i-1);
            critical.battSOC(i) = critical.battSOC(i-1);
            
        elseif current.in_out_batt_critical(i) > 0       % discharging
            
            %limits
            critical.EoutMax_kineticmodellimit(i) = deltaT*(k*critical.Q_1(i-1)*e_raised+critical.Q_total(i-1)*k*c*(1-e_raised))/(1-e_raised+c*(k*deltaT-1+e_raised));
            critical.EoutMax_SOCminlimit(i) = (critical.battSOC(i-1)-(1-battery.MDOD))*Qmax_battbank;
            critical.battEoutMax(i)= min([critical.EoutMax_kineticmodellimit(i),critical.EoutMax_SOCminlimit(i),current.in_out_batt_critical(i)])*battery.discharge_eff;
            

            critical.Q_1(i) = critical.Q_1(i-1)*e_raised+(critical.Q_total(i-1)*k*c - critical.battEoutMax(i))*(1-e_raised)/k - critical.battEoutMax(i)*c*(k*deltaT-1+e_raised)/k;
            critical.Q_2(i) = critical.Q_2(i-1)*e_raised+(critical.Q_total(i-1)*(1-c))*(1-e_raised) - critical.battEoutMax(i)*(1-c)*(k*deltaT-1+e_raised)/k;
            critical.Q_total(i) = critical.Q_1(i) + critical.Q_2(i);
            critical.battSOC(i) = abs(critical.Q_total(i)/Qmax_battbank);          % + battSOC(i-1);
            
            if critical.battSOC(i-1) - critical.battSOC(i) > (1- battery.MDOD)
                critical.battSOC(i) = critical.battSOC(i-1) - critical.battSOC(i);
                
            else
                %critical.battSOC(i) = critical.battSOC(i-1);
                critical.battSOC(i) = 1 - battery.MDOD;
                % is this right
            end
            
            
            %in current (A)
            %battEoutMax is negative
            
            % BATT EOUT OF Q1 IN TERMS OF NON-SERVED 
            %network.critical.nonserved(i) = current.in_out_batt_critical
            
        end
        
        network.critical.nonserved(i) = (loadprof.final.critical(i,end) - battery.voltage*critical.Q_1(i)-avail_solar(i))*ge(loadprof.final.critical(i,end),(avail_solar(i)+battery.voltage*critical.Q_1(i)));

        
    end
    
    network.critical.reliability = 1 - (sum(network.critical.nonserved)/ sum(loadprof.final.critical(:,end)));
    %j=j+1;


% COUNT CYCLES USING RAINFLOW ALGORITHM 
cycles.critical_reshape = reshape(critical.battSOC,24*365/cycles.num,cycles.num);
cycles.critical_reshape = cycles.critical_reshape';
for k=1:cycles.num
     cycles.critical_mean(k) = mean(cycles.critical_reshape(k,:));  
end

%cycle_total = rainflow(total.cycles_mean);
%hourly cycle
cycles.critical_1 = rainflow(critical.battSOC(:,end));
cycles.critical_1 = cycles.critical_1';

%daily cycle
cycles.critical_2 = rainflow(cycles.critical_mean);
cycles.critical_2 = cycles.critical_2';
%1,: = amplitude, 2: mean value, 3: is it 0.5 or 1 cycle
cycles.final_critical = sum(cycles.critical_2(:,end));


% CALCULATE RELIABILITY
%CALCULATE critical SERVED
reliability_critical = network.critical.reliability;
critical_amount_served = sum(loadprof.final.critical(:,end)) - sum(network.critical.nonserved);


%% RAINFLOW

% COUNT CYCLES USING RAINFLOW ALGORITHM 
cycles.total_reshape = reshape(total.battSOC,[12,365*2]); %reshape from hourly to daily 
cycles.total_reshape_2 = reshape(total.battSOC,[24,365]); %reshape from hourly to daily 

for k=1:365*2
     cycles.total_mean(k) = mean(cycles.total_reshape(:,k));
     %rf_out = rainflow(cycles.total_mean(:,k));
     %cycles.total_use(k) = mean(rf_out(end,:));
end

%cycle_total = rainflow(total.cycles_mean);
%hourly cycle
%cycles.total_1 = rainflow(total.battSOC(:,end)); %8760,1:8760);

%daily cycle
%cycles.total_2 = rainflow(cycles.total_mean);
%cycles.total_2 = cycles.total_2';
%1,: = amplitude, 2: mean value, 3: is it 0.5 or 1 cycle
cycles.total_use = rainflow(cycles.total_mean);
cycles.final_total = sum(cycles.total_use(end,:));

%% Plotting for Thesis 

% figure;
% subplot(2,1,1), plot(critical.EoutMax_kineticmodellimit(:,end));
% hold on;
% plot(critical.EoutMax_SOCminlimit(:,end),'r');
% plot(critical.battEoutMax(:,end),'g');
% plot(max(current.in_out_batt_critical(:,end),0),'c');
% title('Critical Load: Battery Discharge Limits');
% legend('KIBAM Limit', 'SOC Min Limit', 'Battery Discharge','Load/PV Current');
% xlabel('Hour of Year');
% ylabel('Ah');
% hold off;
% 
% subplot(2,1,2), plot(critical.Ebatt_max_in_lim1(:,end));
% hold on;
% plot(critical.Ebatt_max_in_limSOCmax(:,end),'r');
% plot(critical.battEinMax(:,end),'g');
% plot(min(current.in_out_batt_critical(:,end),0),'c');
% title('Critical Load: Battery Charge Limits');
% legend('KIBAM Limit', 'SOC Max Limit', 'Battery Charge','Load/PV Current');
% xlabel('Hour of Year');
% ylabel('Ah');
% hold off;
% set(gca,'fontsize',15);

% figure;
% subplot(2,1,1),plot(available_solar);
% hold on;
% plot(-loadprof.final.total(:,end));
% %plot(total.Q_total);
% legend('120W Solar Output','5 Household Total Load Profile');
% hold off;
% set(gca,'fontsize',15);
% xlabel('Hour of Year');
% ylabel('Watts');
% 
% 
% subplot(2,1,2),plot(total.Q_1);
% hold on;
% plot(total.Q_2);
% %plot(total.Q_total);
% legend('Q1','Q2');
% xlabel('Hour of Year');
% ylabel('Amps');
% hold off;
% set(gca,'fontsize',15);
% 
% 
% figure;
% subplot(2,1,1), plot(total.EoutMax_kineticmodellimit);
% hold on;
% plot(total.EoutMax_SOCminlimit);
% %plot(total.battEoutMax);
% plot(max(current.in_out_batt_total,0));
% %title('Total Load: Battery Discharge Limits');
% legend('KiBaM Discharge Limit', 'MDOD SOC Limit', 'Net Load Profile');
% hold off;
% set(gca,'fontsize',15);
% xlabel('Hour of Year');
% ylabel('Amps');
% 
% subplot(2,1,2), plot(total.Ebatt_max_in_lim1);
% hold on;
% plot(total.Ebatt_max_in_limSOCmax);
% %plot(total.battEinMax);
% plot(min(current.in_out_batt_total,0));
% %title('Total Load: Battery Charge Limits');
% legend('KiBaM Charge Limit', 'MLOC SOC Limit', 'Net Load Profile');
% hold off;
% set(gca,'fontsize',15);
% xlabel('Hour of Year');
% ylabel('Amps');


%% limits from kibam 
% charge_str = {'Solar','KiBaM Limit','MLOC Limit'};
% discharge_str = {'Load','KiBaM Limit', 'MDOD Limit'};
% charge_plot = [sum(SOC.charge_total_limit==1), sum(SOC.charge_total_limit==2),sum(SOC.charge_total_limit==3)];
% discharge_plot = [sum(SOC.discharge_total_limit==1), sum(SOC.discharge_total_limit==2),sum(SOC.discharge_total_limit==3)];
% figure;
% bar(charge_plot);
% set(gca, 'XTickLabel',charge_str, 'XTick',1:numel(charge_str))
% set(gca,'fontsize',15);
% 
% figure;
% bar(discharge_plot);
% set(gca, 'XTickLabel',discharge_str, 'XTick',1:numel(discharge_str))
% set(gca,'fontsize',15);
% 
% figure;
% subplot(2,1,1), hist(total.battSOC);
% set(get(gca,'child'),'FaceColor','b','EdgeColor','c');
% set(gca,'fontsize',15);
% xlabel('State of Charge (%)');
% ylabel('Hourly Counts');
% subplot(2,1,2), plot(total.battSOC);
% set(gca,'fontsize',15);
% xlabel('Hour of Year');
% ylabel('State of Charge (%)');

%% PLOTTING THESIS 


% figure(13);
% hold on;
% plot(loadprof.final.total(:,end)-network.total.amountserved.cumulative(:,end));
% xlabel('Hour of Year');
% ylabel('Watts');
% title('Amount of Non Served Energy Across Network, total + Critical Load')
% plot(loadprof.final.critical(:,end)-network.critical.amountserved.cumulative(:,end));
% legend('total NonServed', 'Critical NonServed');
% hold off;


% figure(9);
% subplot(2,2,3,'align'), plotyy(1:num.consumers,network.total.reliability,1:num.consumers,network.critical.reliability);
% xlabel('Number of B Customers');
% ylabel('Avg Yearly Reliability, Percent');
% title(strcat('Load Served with Decreasing Reliability, A+B Customers, for N=',num.graph))
% legend('total Reliability', 'Critical Reliability');
%
% subplot(2,2,1,'align'),plot(1:length(SOC.total),SOC.total(:,num.consumers)/battery.size);
% xlabel('Hour');
% ylabel('Remaining SOC');
% title(strcat(' total SOC Remaining for N=', num.graph))
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
% ylabel('Amt of total Spillage in Wh');
% title(strcat('Spillage for Serving total Load for N=',num.graph));
%
%
% figure(11);
% subplot(1,2,1,'align')
%
% plot(1:length(network.total.nonserved),network.total.nonserved(:,num.consumers), 'c');
% xlabel('Time');
% ylabel('Amt of total NonServed in Wh');
% title(strcat('total Non Served Energy, Cumulative for N=',num.graph));
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
% legend('total Load', 'Critical Load');
% title('Number of Cycles when SOC was above 70% serving Critical and total Load for Increasing N');
% %% thesis 
% 
% figure;
% plot(loadprof.final.total(:,end),'c');
% hold on;
% plot(network.total.nonserved,'r');
% set(gca,'fontsize',15);
% legend('Total Load','Non-Served Load');
% xlabel('Hour of Year');
% ylabel('Watts');
% 
% figure;
% plot(loadprof.final.critical(:,end),'c');
% hold on;
% plot(network.critical.nonserved,'r');
% set(gca,'fontsize',15);
% legend('Total Load','Non-Served Load');
% xlabel('Hour of Year');
% ylabel('Watts');


% figure;
% subplot(2,1,1),scatter(1:8760,1-network.total.nonserved./loadprof.final.total(:,end),'r');
% ylim([0,1.01]);
% xlabel('Hour of Year');
% ylabel('Reliability (%)');
% set(gca,'fontsize',15);
% 
% subplot(2,1,2),plot(loadprof.final.total(:,end));
% hold on;
% plot(1:8760,network.total.nonserved,'r');
% xlabel('Hour of Year');
% ylabel('Power (Watts)');
% legend('Total Load Profile','Total Non-Served Energy');
% set(gca,'fontsize',15);
% 
%  % reliability plots 
% figure;
% subplot(2,1,1),scatter(1:8760,1-network.critical.nonserved./loadprof.final.critical(:,end),'r');
% ylim([0,1.01]);
% xlabel('Hour of Year');
% ylabel('Reliability (%)');
% set(gca,'fontsize',15);
% 
% subplot(2,1,2),plot(loadprof.final.critical(:,end));
% hold on;
% plot(1:8760,network.critical.nonserved,'r');
% xlabel('Hour of Year');
% ylabel('Power (Watts)');
% legend('Critical Load Profile','Critical Non-Served Energy');
% set(gca,'fontsize',15);

% crit_rel = 1-network.critical.nonserved./loadprof.final.critical(:,end);
% tot_rel = 1-network.total.nonserved./loadprof.final.total(:,end);
% 
% figure;
% hist(crit_rel,100);
% 
% figure;
% hist(tot_rel,100);


% SOC thresholds plot 

% red = NaN(length(total.battSOC),1);
% yellow = NaN(length(total.battSOC),1);
% green = NaN(length(total.battSOC),1);
% 
% for k=1:length(total.battSOC)
%     
%     if total.battSOC(k) <= .55
%         red(k) = total.battSOC(k);
%     elseif total.battSOC(k) >= .55 && total.battSOC(k) <= .75
%         yellow(k) = total.battSOC(k);
%     elseif total.battSOC(k) >= .75
%         green(k) = total.battSOC(k);
%     end
%     
%     
% end


%[red_r red_c red_v] = find(total.battSOC(total.battSOC<.55 & total.battSOC>=.4));
%[yel_r yel_c yel_v]= find(total.battSOC(total.battSOC<.75 & total.battSOC>=.55));
%[gre_r gre_c gre_v] = find(total.battSOC(total.battSOC<=1 & total.battSOC>.75));

% 
% figure;
% plot(1:length(green),green,'g');
% hold on;
% plot(1:length(yellow),yellow,'y');
% plot(1:length(red),red,'r');
% set(gca,'color',[.75 .75 .75]);
% xlabel('Hour of Year');
% ylabel('State-of-Charge (%)');
% set(gca,'fontsize',15);

end




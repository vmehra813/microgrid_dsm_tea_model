function [reliability_critical, critical_amount_served,reliability_total, total_amount_served] = connect_ulink_customer_gen_sizing_kibam_pb_acid(available_solar,battery,loadprof,num,batt)

%% Given remaining SOC from A load, ability to serve another customer.

%right now code can take in multiple a generators but is only using 1 a box
%remaining SOC
%incorporate reliability preference
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
t=battery.t1;
Qmax_batt = (batt*(1-exp(-t*k)*(1-c)+t*c*k))/(t*c*k);
Qmax_battbank = Qmax_batt*num_batt; % Qmax_battbank is max capacity of the battery bank in Ah


e_raised = exp(-k*deltaT);


%Initialize Matrices
total.Q_1 = zeros(size(loadprof.final.total)); %some have an added a row bc of t+1 iteration later on
total.Q_2 = zeros(size(loadprof.final.total));
%battEcum = zeros(size(loadprof.final.total));
total.battEoutMax = zeros(size(loadprof.final.total));
total.battEinMax = zeros(size(loadprof.final.total));
total.Ebatt_max_in_lim1 = zeros(size(loadprof.final.total));
total.Ebatt_max_in_lim2 = zeros(size(loadprof.final.total));
total.Ebatt_max_in_lim3 = zeros(size(loadprof.final.total));
total.Ebatt_max_in_limSOCmax = zeros(size(loadprof.final.total));
total.EoutMax_kineticmodellimit = zeros(size(loadprof.final.total));
total.battSOC = zeros(size(loadprof.final.total));
total.Q_total = zeros(size(loadprof.final.total));


critical.Q_1 = zeros(size(loadprof.final.total)); %some have an added a row bc of t+1 iteration later on
critical.Q_2 = zeros(size(loadprof.final.total));
%battEcum = zeros(size(loadprof.final.total));
critical.battEoutMax = zeros(size(loadprof.final.total));
critical.battEinMax = zeros(size(loadprof.final.total));
critical.Ebatt_max_in_lim1 = zeros(size(loadprof.final.total));
critical.Ebatt_max_in_lim2 = zeros(size(loadprof.final.total));
critical.Ebatt_max_in_lim3 = zeros(size(loadprof.final.total));
critical.Ebatt_max_in_limSOCmax = zeros(size(loadprof.final.total));
critical.EoutMax_kineticmodellimit = zeros(size(loadprof.final.total));
critical.battSOC = zeros(size(loadprof.final.total));
critical.Q_total = zeros(size(loadprof.final.total));


battery.initial_capacity = .8;

%% Defining Battery Charge / Discharge [TOTAL vs NON CRITICAL]
%these are matrices
avail_solar = repmat(available_solar,1,num.consumers);

loadprof.final.non_critical = loadprof.final.total;

%power.in_out_batt_total = avail_solar - loadprof.final.total;
power.in_out_batt_total =  loadprof.final.total - avail_solar; %negative is charging
power.in_out_batt_critical =  loadprof.final.critical - avail_solar;
current.in_out_batt_total = power.in_out_batt_total / battery.voltage;
current.in_out_batt_critical = power.in_out_batt_critical / battery.voltage;


%% Initial Condition
%first hour is the same
%total
total.Q_1(1,:) = c*battery.initial_capacity*Qmax_battbank; %these are the starting values of Q_1 and Q_2 for the first iteration - this is a property of the battery bank
total.Q_2(1,:) = (1-c)*battery.initial_capacity*Qmax_battbank;

%battEcum(1,:) = 0;
total.Q_total(1,:) = total.Q_1(1,:) + total.Q_2(1,:);
total.battSOC(1,:) = total.Q_total(1,:)/(Qmax_battbank);

%removed -1*
% total.EoutMax_kineticmodellimit(1,:) = -1*battery.discharge_eff*deltaT*((-k*c*Qmax_battbank + k*total.Q_1(1)*e_raised + total.Q_total(1)*k*c*(1-e_raised))/...
%     (1-e_raised+c*(k*deltaT-1+e_raised)));

%EoutMax_SOCminlimit = (total.battSOC(1)-battery.MDOD)*Qmax_battbank;
% total.EoutMax_SOCminlimit(1,:) = (total.battSOC(1,:)-battery.MDOD)*Qmax_battbank;
% 
% total.battEoutMax(1,:)= min(total.EoutMax_kineticmodellimit(1,:),total.EoutMax_SOCminlimit(1,:));
% 
% 
% total.Ebatt_max_in_lim1(1,:) = deltaT*k*total.Q_1(1)*e_raised+total.Q_total(1,:)*k*c*(1-e_raised)/(1-e_raised+c*(k*deltaT-1+e_raised));
% total.Ebatt_max_in_lim2(1,:) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-total.Q_total(1,:))/deltaT;
% total.Ebatt_max_in_lim3(1,:) = deltaT*num_batt*Imax/1;
% %Ebatt_max_in_limSOCmax = (BESS.SOCmax-total.battSOC(1))*Qmax_battbank;
% %Ebatt_max_in_limSOCmax = (1-total.battSOC(1))*Qmax_battbank;
% total.Ebatt_max_in_limSOCmax(1,:) = (1-total.battSOC(1,:))*Qmax_battbank;
% 
% % energy into the battery during charge (Wh)
% total.battEinMax(1,:) = min([total.Ebatt_max_in_lim1(1,:), total.Ebatt_max_in_lim2(1,:), total.Ebatt_max_in_lim3(1,:), total.Ebatt_max_in_limSOCmax(1,:)])/battery.charge_eff;
% 

%Critical
critical.Q_1(1,:) = c*battery.initial_capacity*Qmax_battbank; %these are the starting values of Q_1 and Q_2 for the first iteration - this is a property of the battery bank
critical.Q_2(1,:) = (1-c)*battery.initial_capacity*Qmax_battbank;

%battEcum(1,:) = 0;
critical.Q_total(1,:) = critical.Q_1(1,:) + critical.Q_2(1,:);
critical.battSOC(1,:) = critical.Q_total(1,:)/(Qmax_battbank);

%removed -1*
% critical.EoutMax_kineticmodellimit(1,:) = -1*battery.discharge_eff*deltaT*((-k*c*Qmax_battbank + k*critical.Q_1(1)*e_raised + critical.Q_total(1)*k*c*(1-e_raised))/...
%     (1-e_raised+c*(k*deltaT-1+e_raised)));
% 
% %EoutMax_SOCminlimit = (critical.battSOC(1)-battery.MDOD)*Qmax_battbank;
% critical.EoutMax_SOCminlimit(1,:) = (critical.battSOC(1,:)-battery.MDOD)*Qmax_battbank;
% 
% critical.battEoutMax(1,:)= max(critical.EoutMax_kineticmodellimit(1,:),critical.EoutMax_SOCminlimit(1,:));
% 
% 
% critical.Ebatt_max_in_lim1(1,:) = deltaT*k*critical.Q_1(1)*e_raised+critical.Q_total(1,:)*k*c*(1-e_raised)/(1-e_raised+c*(k*deltaT-1+e_raised));
% critical.Ebatt_max_in_lim2(1,:) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-critical.Q_total(1,:))/deltaT;
% critical.Ebatt_max_in_lim3(1,:) = deltaT*num_batt*Imax/1;
% %Ebatt_max_in_limSOCmax = (BESS.SOCmax-critical.battSOC(1))*Qmax_battbank;
% %Ebatt_max_in_limSOCmax = (1-critical.battSOC(1))*Qmax_battbank;
% critical.Ebatt_max_in_limSOCmax(1,:) = (1-critical.battSOC(1,:))*Qmax_battbank;
% 
% % energy into the battery during charge (Wh)
% critical.battEinMax(1,:) = min([critical.Ebatt_max_in_lim1(1,:), critical.Ebatt_max_in_lim2(1,:), critical.Ebatt_max_in_lim3(1,:), critical.Ebatt_max_in_limSOCmax(1,:)])/battery.charge_eff;
% 


%total
%initialize vectors
network.total.reliability = zeros(num.consumers,1);
network.total.amountserved.cumulative = zeros(length(loadprof.final.total),num.consumers);
network.total.nonserved = zeros(length(loadprof.final.total),num.consumers);
%network.total.nonserved = [];
%network.total.amountserved.cumulative(1) = 20;
network.total.solar_load_coincidence = zeros(length(loadprof.final.total),num.consumers);


SOC.total = zeros(length(loadprof.final.total), num.consumers);
SOC.total_spillage = zeros(length(loadprof.final.total),num.consumers);
SOC.total_cycles = zeros(num.consumers,1);

%Critical
%initialize vectors
network.critical.reliability = zeros(num.consumers,1);
network.critical.amountserved.cumulative = zeros(length(loadprof.final.total),num.consumers);
network.critical.nonserved = zeros(length(loadprof.final.total),num.consumers);
%network.critical.nonserved = [];
%network.total.amountserved.cumulative(1) = 20;
network.critical.solar_load_coincidence = zeros(length(loadprof.final.total),num.consumers);


SOC.critical = zeros(length(loadprof.final.total), num.consumers);
SOC.critical_spillage = zeros(length(loadprof.final.total),num.consumers);
SOC.critical_cycles = zeros(num.consumers,1);

%count variables
% zero = current; 1 = KIBAM; 2 = SOC
SOC.charge_total_limit = zeros(length(loadprof.final.total),num.consumers);
SOC.charge_critical_limit = zeros(length(loadprof.final.total),num.consumers);
SOC.discharge_total_limit = zeros(length(loadprof.final.total),num.consumers);
SOC.discharge_critical_limit = zeros(length(loadprof.final.total),num.consumers);

%% total (KIBAM) (now Total Load)
%clear j i;
%j=1;
%while j<= num.consumers
%for j=1:num.consumers %column is B customer
for j=num.consumers
    
    for i=1
        
        if current.in_out_batt_total(i,j) < 0     %charging
            
            % need other 1st condition to include solar/load for hour 1
            % total.Q_1(i) = critical.Q_1(i)*e_raised+(critical.Q_total(i)*k*c+current.in_out_batt_total(i))*(1-e_raised)/k+power.in_out_batt_total(i)*c*(k*deltaT-1+e_raised)/k;
            % critical.Q_2(i) = critical.Q_2(i)*e_raised+(critical.Q_total(i)*(1-c))*(1-e_raised)+current.in_out_batt_total(i)*(1-c)*(k*deltaT-1+e_raised)/k;
            % critical.Q_total(i) = critical.Q_1(i) + critical.Q_2(i);
            % critical.battSOC(i) = critical.Q_total(i)/Qmax_battbank;
            
            total.Ebatt_max_in_lim1(i,j) = (battery.discharge_eff*deltaT*((-k*c*Qmax_battbank + k*total.Q_1(i,j)*e_raised + total.Q_total(i,j)*k*c*(1-e_raised))/...
                (1-e_raised+c*(k*deltaT-1+e_raised))));
            %total.Ebatt_max_in_lim2(i,j) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-total.Q_total(i,j))/deltaT;
            %total.Ebatt_max_in_lim3(i,j) = deltaT*num_batt*Imax/1;
            total.Ebatt_max_in_limSOCmax(i,j) = -1*(battery.MLOC-total.battSOC(i,j))*Qmax_battbank;
            %energy into the battery during charge (Wh)
            %total.battEinMax(i,j) = min([total.Ebatt_max_in_lim1(i,j), total.Ebatt_max_in_lim2(i,j), total.Ebatt_max_in_lim3(i,j), total.Ebatt_max_in_limSOCmax(i,j)])/battery.charge_eff;
            total.battEinMax(i,j) = max([total.Ebatt_max_in_lim1(i,j),total.Ebatt_max_in_limSOCmax(i,j),current.in_out_batt_total(i,j)]);
            
            
            % total spillage (Wh)
            SOC.total_spillage(i,j) = abs(avail_solar(i,j)-loadprof.final.total(i,j)-battery.voltage*total.Q_1(i,j))*ge(avail_solar(i,j),loadprof.final.total(i,j)+battery.voltage*total.Q_1(i,j));
            
            
        elseif current.in_out_batt_total(i,j) > 0 %discharging
            
            % need other 1st condition to include solar/load for hour 1
            
            % critical.Q_1(i,j) = critical.Q_1(i,j)*e_raised+(critical.Q_total(i,j)*k*c+current.in_out_batt_total(i,j))*(1-e_raised)/k+power.in_out_batt_total(i,j)*c*(k*deltaT-1+e_raised)/k;
            % critical.Q_2(i,j) = critical.Q_2(i,j)*e_raised+(critical.Q_total(i,j)*(1-c))*(1-e_raised)+current.in_out_batt_total(i,j)*(1-c)*(k*deltaT-1+e_raised)/k;
            % critical.Q_total(i,j) = critical.Q_1(i,j) + critical.Q_2(i,j);
            % critical.battSOC(i,j) = critical.Q_total(i,j)/Qmax_battbank;
          
            total.EoutMax_kineticmodellimit(i,j) = deltaT*(k*total.Q_1(i,j)*e_raised+total.Q_total(i,j)*k*c*(1-e_raised))/(1-e_raised+c*(k*deltaT-1+e_raised));
             
            total.EoutMax_SOCminlimit(i,j) = (total.battSOC(i,j)-(1-battery.MDOD))*Qmax_battbank;
            % energy out from the battery during discharge (Wh)
            total.battEoutMax(i,j)= min([total.EoutMax_kineticmodellimit(i,j),total.EoutMax_SOCminlimit(i,j),current.in_out_batt_total(i,j)]);
                    
              
            
            % total non-served energy (Wh)
            network.total.nonserved(i,j) = (loadprof.final.total(i,j) - battery.voltage*total.Q_1(i,j)-avail_solar(i,j))*ge(loadprof.final.total(i,j),(avail_solar(i,j)+battery.voltage*total.Q_1(i,j)));
            
                 
            
        end
        
    end
    
    for i=2:length(loadprof.final.total)
        if current.in_out_batt_total(i,j) < 0       % charging
            
            %limits
            
            total.Ebatt_max_in_lim1(i,j) = (battery.discharge_eff*deltaT*((-k*c*Qmax_battbank + k*total.Q_1(i-1,j)*e_raised + total.Q_total(i-1,j)*k*c*(1-e_raised))/...
                (1-e_raised+c*(k*deltaT-1+e_raised))));
            %total.Ebatt_max_in_lim2(i,j) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-total.Q_total(i-1,j))/deltaT;
            %total.Ebatt_max_in_lim3(i,j) = deltaT*num_batt*Imax*battery.voltage/1;
            total.Ebatt_max_in_limSOCmax(i,j) = -1*(battery.MLOC-total.battSOC(i-1,j))*Qmax_battbank;
            %total.battEinMax(i,j) = min([total.Ebatt_max_in_lim1(i,j), total.Ebatt_max_in_lim2(i,j), total.Ebatt_max_in_lim3(i,j), total.Ebatt_max_in_limSOCmax(i,j)])/battery.charge_eff;
            %total.battEinMax(i,j) = min([total.Ebatt_max_in_lim1(i,j), total.Ebatt_max_in_lim3(i,j), total.Ebatt_max_in_limSOCmax(i,j)])/battery.charge_eff;
            %changed to min
            total.battEinMax(i,j) = max([total.Ebatt_max_in_lim1(i,j), total.Ebatt_max_in_limSOCmax(i,j),current.in_out_batt_total(i,j)])/battery.charge_eff;

            if total.battEinMax(i,j) == total.Ebatt_max_in_lim1(i,j)
                SOC.charge_total_limit(i,j) = 1;
                
            elseif total.battEinMax(i,j) == total.Ebatt_max_in_limSOCmax(i,j)
                SOC.charge_total_limit(i,j) = 2;
                
            end
                     
            
            total.Q_1(i,j) = total.Q_1(i-1,j)*e_raised+(total.Q_total(i-1,j)*k*c-total.battEinMax(i,j))*(1-e_raised)/k-total.battEinMax(i,j)*c*(k*deltaT-1+e_raised)/k;
            total.Q_2(i,j) = total.Q_2(i-1,j)*e_raised+(total.Q_total(i-1,j)*(1-c))*(1-e_raised)-total.battEinMax(i,j)*(1-c)*(k*deltaT-1+e_raised)/k;
            total.Q_total(i,j) = total.Q_1(i,j) + total.Q_2(i,j);
            total.battSOC(i,j) = abs(total.Q_total(i,j)/Qmax_battbank);      % + battSOC(i-1);
            
            
            if total.battSOC(i,j)+total.battSOC(i-1,j) < battery.MLOC
                total.battSOC(i,j) = total.battSOC(i,j)+total.battSOC(i-1,j);
            else
                %total.battSOC(i,j) = total.battSOC(i-1);
                total.battSOC(i,j) = battery.MLOC;
            end
            
            %spillage: if total.Q_1 is larger than batt_e_in_max
            %set to previous condition
            %limits
            
            SOC.total_spillage(i,j) = abs(avail_solar(i,j)-loadprof.final.total(i,j)-battery.voltage*total.Q_1(i,j))*ge(avail_solar(i,j),loadprof.final.total(i,j)+battery.voltage*total.Q_1(i,j));
            
        elseif current.in_out_batt_total(i,j) == 0
            
            %critical.Q_1(i,j) = critical.Q_1(i-1);
            %critical.Q_2(i,j) = critical.Q_2(i-1);
            %critical.Q_total(i,j) = critical.Q_total(i-1);
            total.battSOC(i,j) = total.battSOC(i-1,j);
            
        elseif current.in_out_batt_total(i,j) > 0       % discharging
            
             %limits
            total.EoutMax_kineticmodellimit(i,j) = deltaT*(k*total.Q_1(i-1,j)*e_raised+total.Q_total(i-1,j)*k*c*(1-e_raised))/(1-e_raised+c*(k*deltaT-1+e_raised));
            total.EoutMax_SOCminlimit(i,j) = (total.battSOC(i-1,j)-(1-battery.MDOD))*Qmax_battbank;
            total.battEoutMax(i,j)= min([total.EoutMax_kineticmodellimit(i,j),total.EoutMax_SOCminlimit(i,j),current.in_out_batt_total(i,j)]);
            
            
            if total.battEoutMax(i,j) == total.EoutMax_kineticmodellimit(i,j)
                SOC.discharge_total_limit(i,j) = 1;
                
            elseif total.battEoutMax(i,j) == total.EoutMax_SOCminlimit(i,j)
                SOC.discharge_total_limit(i,j) = 2;
                
            end
            
            
            
            total.Q_1(i,j) = total.Q_1(i-1,j)*e_raised+(total.Q_total(i-1,j)*k*c - total.battEoutMax(i,j))*(1-e_raised)/k - total.battEoutMax(i,j)*c*(k*deltaT-1+e_raised)/k;
            total.Q_2(i,j) = total.Q_2(i-1,j)*e_raised+(total.Q_total(i-1,j)*(1-c))*(1-e_raised) - total.battEoutMax(i,j)*(1-c)*(k*deltaT-1+e_raised)/k;
            total.Q_total(i,j) = total.Q_1(i,j) + total.Q_2(i,j);
            total.battSOC(i,j) = abs(total.Q_total(i,j)/Qmax_battbank);          % + battSOC(i-1);
            
            if total.battSOC(i-1,j) - total.battSOC(i,j) > (1- battery.MDOD)
                total.battSOC(i,j) = total.battSOC(i-1,j) - total.battSOC(i,j);
                
            else
                %total.battSOC(i,j) = total.battSOC(i-1);
                total.battSOC(i,j) = 1 - battery.MDOD;
                % is this right
            end
            
            
            %in current (A)
            %battEoutMax is negative
            network.total.nonserved(i,j) = (loadprof.final.total(i,j) - battery.voltage*total.Q_1(i,j)-avail_solar(i,j))*ge(loadprof.final.total(i,j),(avail_solar(i,j)+battery.voltage*total.Q_1(i,j)));
%             nonserved = loadprof.final.total(i,j) + total.battEoutMax(i,j) - avail_solar(i,j);
%             if nonserved > 0
%                 network.total.nonserved(i,j) = nonserved;
%             end
            
            
        end
        
    end
    
    network.total.reliability(j) = 1-(sum(network.total.nonserved(:,j))/ sum(loadprof.final.total(:,j)));
    
    
    
    %j=j+1;
    
    
end

% CALCULATE total RELIABILITY
% CALCULATE total AMOUNT SERVED
reliability_total = network.total.reliability(end);
total_amount_served = sum(loadprof.final.total(:,end)) - sum(network.total.nonserved(:,end));


%% Critical (KIBAM)
%clear j i;
%j=1;
%while j<= num.consumers
%for j=1:num.consumers %column is B customer
for j=num.consumers 

    for i=1
        
        if current.in_out_batt_critical(i,j) < 0     %charging
            
            % need other 1st condition to include solar/load for hour 1
            % critical.Q_1(i) = critical.Q_1(i)*e_raised+(critical.Q_total(i)*k*c+current.in_out_batt_critical(i))*(1-e_raised)/k+power.in_out_batt_critical(i)*c*(k*deltaT-1+e_raised)/k;
            % critical.Q_2(i) = critical.Q_2(i)*e_raised+(critical.Q_total(i)*(1-c))*(1-e_raised)+current.in_out_batt_critical(i)*(1-c)*(k*deltaT-1+e_raised)/k;
            % critical.Q_total(i) = critical.Q_1(i) + critical.Q_2(i);
            % critical.battSOC(i) = critical.Q_total(i)/Qmax_battbank;
            
            critical.Ebatt_max_in_lim1(i,j) = (battery.discharge_eff*deltaT*((-k*c*Qmax_battbank + k*critical.Q_1(i,j)*e_raised + critical.Q_total(i,j)*k*c*(1-e_raised))/...
                (1-e_raised+c*(k*deltaT-1+e_raised))));
            
            %critical.Ebatt_max_in_lim2(i,j) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-critical.Q_total(i,j))/deltaT;
            %critical.Ebatt_max_in_lim3(i,j) = deltaT*num_batt*Imax/1;
            critical.Ebatt_max_in_limSOCmax(i,j) = -1*(battery.MLOC-critical.battSOC(i,j))*Qmax_battbank;
            % energy into the battery during charge (Wh)
            %critical.battEinMax(i,j) = min([critical.Ebatt_max_in_lim1(i,j), critical.Ebatt_max_in_lim2(i,j), critical.Ebatt_max_in_lim3(i,j), critical.Ebatt_max_in_limSOCmax(i,j)])/battery.charge_eff;
            critical.battEinMax(i,j) = max([critical.Ebatt_max_in_lim1(i,j), critical.Ebatt_max_in_limSOCmax(i,j), current.in_out_batt_critical(i,j)])/battery.charge_eff;
           
    
            % critical spillage (Wh)
            %SOC.critical_spillage(i,j) = abs(avail_solar(i,j)-loadprof.final.total(i,j)+battery.voltage*critical.battEinMax(i,j))*ge(avail_solar(i,j),loadprof.final.critical(i,j)-battery.voltage*critical.battEinMax(i,j));
            SOC.critical_spillage(i,j) = abs(avail_solar(i,j)-loadprof.final.total(i,j)-battery.voltage*critical.Q_1(i,j))*ge(avail_solar(i,j),loadprof.final.critical(i,j)+battery.voltage*critical.Q_1(i,j));

            
        elseif current.in_out_batt_critical(i,j) > 0 %discharging
            
            % need other 1st condition to include solar/load for hour 1
            
            % critical.Q_1(i,j) = critical.Q_1(i,j)*e_raised+(critical.Q_total(i,j)*k*c+current.in_out_batt_critical(i,j))*(1-e_raised)/k+power.in_out_batt_critical(i,j)*c*(k*deltaT-1+e_raised)/k;
            % critical.Q_2(i,j) = critical.Q_2(i,j)*e_raised+(critical.Q_total(i,j)*(1-c))*(1-e_raised)+current.in_out_batt_critical(i,j)*(1-c)*(k*deltaT-1+e_raised)/k;
            % critical.Q_total(i,j) = critical.Q_1(i,j) + critical.Q_2(i,j);
            % critical.battSOC(i,j) = critical.Q_total(i,j)/Qmax_battbank;
            
            
            critical.EoutMax_kineticmodellimit(i,j) = deltaT*(k*total.Q_1(i,j)*e_raised+total.Q_total(i,j)*k*c*(1-e_raised))/(1-e_raised+c*(k*deltaT-1+e_raised));
            
            critical.EoutMax_SOCminlimit(i,j) = (critical.battSOC(i,j)-(1-battery.MDOD))*Qmax_battbank;
            % energy out from the battery during discharge (Wh)
            critical.battEoutMax(i,j)= min([critical.EoutMax_kineticmodellimit(i,j),critical.EoutMax_SOCminlimit(i,j),current.in_out_batt_critical(i,j)]);
            
            % critical non-served energy (Wh)
            %network.critical.nonserved(i,j) = (loadprof.final.critical(i,j) - battery.voltage*critical.battEoutMax(i,j)-avail_solar(i,j))*ge(loadprof.final.critical(i,j),(avail_solar(i,j)+battery.voltage*critical.battEoutMax(i,j)));
             network.critical.nonserved(i,j) = (loadprof.final.critical(i,j) - battery.voltage*critical.Q_1(i,j)-avail_solar(i,j))*ge(loadprof.final.critical(i,j),(avail_solar(i,j)+battery.voltage*critical.Q_1(i,j)));

            
            
        end
        
    end
    
    for i=2:length(loadprof.final.total)
        if current.in_out_batt_critical(i,j) < 0       % charging
            
            critical.Ebatt_max_in_lim1(i,j) = (battery.discharge_eff*deltaT*((-k*c*Qmax_battbank + k*critical.Q_1(i-1,j)*e_raised + critical.Q_total(i-1,j)*k*c*(1-e_raised))/...
                (1-e_raised+c*(k*deltaT-1+e_raised))));
            %critical.Ebatt_max_in_lim2(i,j) = deltaT*(1-exp(-alpha_c*deltaT))*(Qmax_battbank-critical.Q_total(i,j))/deltaT;
            %critical.Ebatt_max_in_lim3(i,j) = deltaT*num_batt*Imax*battery.voltage/1;
            critical.Ebatt_max_in_limSOCmax(i,j) = -1*(battery.MLOC-critical.battSOC(i-1,j))*Qmax_battbank;
            %critical.battEinMax(i,j) = min([critical.Ebatt_max_in_lim1(i,j), critical.Ebatt_max_in_lim2(i,j), critical.Ebatt_max_in_lim3(i,j), critical.Ebatt_max_in_limSOCmax(i,j)])/battery.charge_eff;
            critical.battEinMax(i,j) = max([critical.Ebatt_max_in_lim1(i,j), critical.Ebatt_max_in_limSOCmax(i,j), current.in_out_batt_critical(i,j)])/battery.charge_eff;
  
            
             if critical.battEinMax(i,j) == critical.Ebatt_max_in_lim1(i,j)
                SOC.charge_critical_limit(i,j) = 1;
                
            elseif critical.battEinMax(i,j) == critical.Ebatt_max_in_limSOCmax(i,j)
                SOC.charge_critical_limit(i,j) = 2;
                
            end
                  
            critical.Q_1(i,j) = critical.Q_1(i-1,j)*e_raised+(critical.Q_total(i-1,j)*k*c-critical.battEinMax(i,j))*(1-e_raised)/k-critical.battEinMax(i,j)*c*(k*deltaT-1+e_raised)/k;
            critical.Q_2(i,j) = critical.Q_2(i-1,j)*e_raised+(critical.Q_total(i-1,j)*(1-c))*(1-e_raised)-critical.battEinMax(i,j)*(1-c)*(k*deltaT-1+e_raised)/k;
            critical.Q_total(i,j) = critical.Q_1(i,j) + critical.Q_2(i,j);
            critical.battSOC(i,j) = abs(critical.Q_total(i,j)/Qmax_battbank);      % + battSOC(i-1);
            
            
            if critical.battSOC(i,j)+critical.battSOC(i-1,j) < battery.MLOC
                critical.battSOC(i,j) = critical.battSOC(i,j)+critical.battSOC(i-1,j);
            else
                %critical.battSOC(i,j) = critical.battSOC(i-1);
                critical.battSOC(i,j) = battery.MLOC;
            end
            
            %spillage: if critical.Q_1 is larger than batt_e_in_max
            %set to previous condition
            %limits
           
            SOC.critical_spillage(i,j) = abs(avail_solar(i,j)-loadprof.final.total(i,j)-battery.voltage*critical.Q_1(i,j))*ge(avail_solar(i,j),loadprof.final.critical(i,j)+battery.voltage*critical.Q_1(i,j));
            
        elseif current.in_out_batt_critical(i,j) == 0
            
            %critical.Q_1(i,j) = critical.Q_1(i-1);
            %critical.Q_2(i,j) = critical.Q_2(i-1);
            %critical.Q_total(i,j) = critical.Q_total(i-1);
            critical.battSOC(i,j) = critical.battSOC(i-1,j);
            
        elseif current.in_out_batt_critical(i,j) > 0       % discharging
            
             %limits
            critical.EoutMax_kineticmodellimit(i,j) = deltaT*(k*critical.Q_1(i-1,j)*e_raised+critical.Q_total(i-1,j)*k*c*(1-e_raised))/(1-e_raised+c*(k*deltaT-1+e_raised));
            critical.EoutMax_SOCminlimit(i,j) = (critical.battSOC(i-1,j)-(1-battery.MDOD))*Qmax_battbank;
            critical.battEoutMax(i,j)= min([critical.EoutMax_kineticmodellimit(i,j),critical.EoutMax_SOCminlimit(i,j),current.in_out_batt_critical(i,j)]);
            
            if critical.battEoutMax(i,j) == critical.EoutMax_kineticmodellimit(i,j)
                SOC.discharge_critical_limit(i,j) = 1;
                
            elseif critical.battEoutMax(i,j) == critical.EoutMax_SOCminlimit(i,j)
                SOC.discharge_critical_limit(i,j) = 2;
                
            end
               
            
            critical.Q_1(i,j) = critical.Q_1(i-1,j)*e_raised+(critical.Q_total(i-1,j)*k*c - critical.battEoutMax(i,j))*(1-e_raised)/k - critical.battEoutMax(i,j)*c*(k*deltaT-1+e_raised)/k;
            critical.Q_2(i,j) = critical.Q_2(i-1,j)*e_raised+(critical.Q_total(i-1,j)*(1-c))*(1-e_raised) - critical.battEoutMax(i,j)*(1-c)*(k*deltaT-1+e_raised)/k;
            critical.Q_total(i,j) = critical.Q_1(i,j) + critical.Q_2(i,j);
            critical.battSOC(i,j) = abs(critical.Q_total(i,j)/Qmax_battbank);          % + battSOC(i-1);
            
            if critical.battSOC(i-1,j) - critical.battSOC(i,j) > (1- battery.MDOD)
                critical.battSOC(i,j) = critical.battSOC(i-1,j) - critical.battSOC(i,j);
                
            else
                %critical.battSOC(i,j) = critical.battSOC(i-1);
                critical.battSOC(i,j) = 1 - battery.MDOD;
                % is this right
            end
            
           
            %in current (A)
            %battEoutMax is negative
            network.critical.nonserved(i,j) = (loadprof.final.critical(i,j) - battery.voltage*critical.Q_1(i,j)-avail_solar(i,j))*ge(loadprof.final.critical(i,j),(avail_solar(i,j)+battery.voltage*critical.Q_1(i,j)));
            %network.critical.nonserved(i,j) = current.in_out_batt_critical 
            
        end
        
    end
    
    network.critical.reliability(j) = 1-(sum(network.critical.nonserved(:,j))/ sum(loadprof.final.critical(:,j)));
    %j=j+1;
end

% CALCULATE RELIABILITY
%CALCULATE critical SERVED
reliability_critical = network.critical.reliability(end);
critical_amount_served = sum(loadprof.final.critical(:,end)) - sum(network.critical.nonserved(:,end));

%% Plotting 

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
% 
% %count / store limits
% 
% figure;
% subplot(2,1,1), plot(total.EoutMax_kineticmodellimit(:,end));
% hold on;
% plot(total.EoutMax_SOCminlimit(:,end),'r');
% plot(total.battEoutMax(:,end),'g');
% plot(max(current.in_out_batt_total(:,end),0),'c');
% title('Total Load: Battery Discharge Limits');
% legend('KIBAM Limit', 'SOC Limit', 'Battery Discharge','Load/PV Current');
% hold off;
% 
% subplot(2,1,2), plot(total.Ebatt_max_in_lim1(:,end));
% hold on;
% plot(total.Ebatt_max_in_limSOCmax(:,end),'r');
% plot(total.battEinMax(:,end),'g');
% plot(min(current.in_out_batt_total(:,end),0),'c');
% title('Total Load: Battery Charge Limits');
% legend('KIBAM Limit', 'SOC Max Limit', 'Battery Charge', 'Load/PV Current');
% xlabel('Hour of Year');
% ylabel('Ah');
% hold off;


%% OLD CODE (before KIBAM)
%
%
%     for i=1  %special case for i=1 to introduce initial state of charge of battery
%         if available_solar(i,1) >= loadprof.final.total(i,j)
%             SOC.total(i,j) = battery.initial_capacity + (available_solar(i,1)-loadprof.final.total(i,j))*battery.charge_eff;
%
%         elseif available_solar(i,1) < loadprof.final.total(i,j)
%             SOC.total(i,j) = battery.initial_capacity - (loadprof.final.total(i,j)-available_solar(i,1))*battery.discharge_eff;
%         end
%     end
%
%     for i=2:length(loadprof.final.total)  %remaining interations look at previous hour
%         if available_solar(i,1) >= loadprof.final.total(i,j)
%             %charging
%             %put into function w/ kibam
%             SOC.total(i,j) = SOC.total(i-1,j)+(available_solar(i,1)-loadprof.final.total(i,j))*battery.charge_eff;
%             network.total.solar_load_coincidence(i,j) = loadprof.final.total(i,j)/available_solar(i,1); % solar / load coincidence factor
%
%             if SOC.total(i,j) < batt
%                 SOC.total(i,j) = SOC.total(i,j);
%             else
%                 SOC.total(i,j) = SOC.total(i-1,j);
%                 SOC.total_spillage(i,j) = available_solar(i,1) - loadprof.final.total(i,j); %spillage
%
%             end
%
%         elseif available_solar(i,1) < loadprof.final.total(i,j)
%
%             SOC.total(i,j) = SOC.total(i-1,j)-(loadprof.final.total(i,j)-available_solar(i,1))*battery.discharge_eff;
%             network.total.solar_load_coincidence(i,j) = available_solar(i,1) / loadprof.final.total(i,j); % solar load coincidence
%
%             if SOC.total(i,j) > battery.MDOD*batt %MDOD constraint
%                 SOC.total(i,j)= SOC.total(i,j);
%             else SOC.total(i,j) = SOC.total(i-1,j); %nonserved
%                 network.total.nonserved(i,j) = loadprof.final.total(i,j) - available_solar(i,1);
%
%             end
%         end
%
%         network.total.amountserved.cumulative(i,j) = available_solar(i,1) + SOC.total(i-1,j) - SOC.total(i,j);
%
%         if network.total.amountserved.cumulative(i,j) > loadprof.final.total(i,j)
%             network.total.amountserved.cumulative(i,j) = loadprof.final.total(i,j);
%         end
%
%
%     end
%
%     network.total.reliability(j) = (sum(network.total.amountserved.cumulative(:,j)))/ (sum(loadprof.final.total(:,j)));
%
%     SOC.total_cycles(j) = length(findpeaks(SOC.total(:,j)/batt,'MinPeakHeight',0.7,'MinPeakDistance',1));
%
%
% end
%
% reliability_critical = network.total.reliability(end);
% critical_amount_served = sum(network.total.amountserved.cumulative(:,end));
%
%
%
%
% %% Critical Served for B
% %
% network.critical.reliability = zeros(num.consumers,1);
% network.critical.amountserved.cumulative = zeros(length(loadprof.final.critical),num.consumers);
% %network.critical.amountserved.cumulative(1) = 20;
%
% network.critical.percentserved.cumulative = zeros(length(loadprof.final.critical),num.consumers);
% network.critical.percentserved.cumulative(1) = 1;
% network.critical.nonserved = zeros(length(loadprof.final.critical),num.consumers);
% network.critical.solar_load_coincidence = zeros(length(loadprof.final.critical),num.consumers);
%
% SOC.critical = zeros(length(loadprof.final.total), num.consumers);
% SOC.critical_spillage = zeros(length(loadprof.final.critical),num.consumers);
% SOC.critical_cycles = zeros(num.consumers,1);
%
%
% for j=1:num.consumers
%     for i=1  %special case for i=1 to introduce initial state of charge of battery
%         if available_solar(i,1) >= loadprof.final.critical(i,j)
%             SOC.critical(i,j) = battery.initial_capacity + (available_solar(i,1)-loadprof.final.critical(i,j))*battery.charge_eff;
%
%         elseif available_solar(i,1) < loadprof.final.critical(i,j)
%             SOC.critical(i,j) = battery.initial_capacity - (loadprof.final.critical(i,j)-available_solar(i,1))*battery.discharge_eff;
%         end
%     end
%
%     for i=2:length(loadprof.final.critical)  %remaining interations look at previous hour
%         if available_solar(i,1) >= loadprof.final.critical(i,j)
%
%             SOC.critical(i,j) = SOC.critical(i-1,j)+(available_solar(i,1)-loadprof.final.critical(i,j))*battery.charge_eff;
%             network.critical.solar_load_coincidence(i,j) = loadprof.final.critical(i,j)/available_solar(i,1);
%             if SOC.critical(i,j) < batt
%                 SOC.critical(i,j) = SOC.critical(i,j);
%             else
%                 SOC.critical(i,j) = SOC.critical(i-1,j);
%                 SOC.critical_spillage(i,j) = available_solar(i,1) - loadprof.final.critical(i,j); %spillage
%
%             end
%
%         elseif available_solar(i,1) < loadprof.final.critical(i,j)
%
%             SOC.critical(i,j) = SOC.critical(i-1,j)-(loadprof.final.critical(i,j)-available_solar(i,1))*battery.discharge_eff;
%             network.critical.solar_load_coincidence(i,j) = available_solar(i,1) / loadprof.final.critical(i,j);
%             if SOC.critical(i,j) > battery.MDOD*batt
%                 SOC.critical(i,j)= SOC.critical(i,j);
%             else
%                 SOC.critical(i,j) = SOC.critical(i-1,j);
%                 network.critical.nonserved(i,j) = loadprof.final.critical(i,j) - available_solar(i,1); %non served
%
%             end
%
%         end
%
%         network.critical.amountserved.cumulative(i,j) = (available_solar(i,1) + SOC.critical(i-1,j) - SOC.critical(i,j));
%
%         if network.critical.amountserved.cumulative(i,j) > loadprof.final.critical(i,j)
%             network.critical.amountserved.cumulative(i,j) = loadprof.final.critical(i,j);
%         end
%
%
%
%     end
%
%     network.critical.reliability(j) = (sum(network.critical.amountserved.cumulative(:,j)))/ (sum(loadprof.final.critical(:,j)));
%     SOC.critical_cycles(j) = length(findpeaks(SOC.critical(:,j)/batt,'MinPeakHeight',0.7,'MinPeakDistance',1));
%
% end
%
% reliability_critical = network.critical.reliability(end);
% critical_amount_served = sum(network.critical.amountserved.cumulative(:,end));

%% PLOTTING


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





end




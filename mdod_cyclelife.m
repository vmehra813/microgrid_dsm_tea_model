% MDOD / Cycle Life

% universal = [30, 1200;
%             50, 500;
%             100, 200];

%10, 8000;
        
amara = [30, 2250;
         50, 1200;
         70, 1000;
         90, 700];

%20, 5000;
exide = [30, 3500;
         40, 2500;
         50, 2000;
         60, 1500;
         70, 1200;
         80, 1000];

% figure;
% scatter(amara(:,1),amara(:,2));
% hold on;
% scatter(exide(:,1),exide(:,2));


all_mdod = [amara(:,1)/100; exide(:,1)/100];
all_cycle = [amara(:,2); exide(:,2)];

fit1 = fit(all_mdod,all_cycle,'exp1');
fit2 = fit(all_mdod,all_cycle,'exp2');
% figure;
% plot(fit1,all_mdod,all_cycle);
% %plot(fit2,all_mdod,all_cycle);
% xlabel('Maximum Depth of Discharge (%)');
% ylabel('Cycle Life');
% %legend('Data from Exide + Amara Raja + SuKam', 'ExpFit');
% %title('Exponential Fit for Battery Costs + Size');
% set(gca,'fontsize',15);




%% thesis plots w/ results

% [~, ~, raw] = xlsread('/Users/vmehra813/Dropbox (MIT)/[uLink]/Simulation/MatLab/uLink/Modeling Analysis/Generation + Reliability Model/thesis_lcoe_biz.xlsx','Sheet2');
% raw = raw(2:end,:);
% 
% %% Create output variable
% data = reshape([raw{:}],size(raw));
% 
% %% Allocate imported array to column variable names
% MDOD = data(:,1);
% Cycles = data(:,2);
% LCOE = data(:,3);
% CapCost = data(:,4);
% Batt = data(:,5);
% PV = data(:,6);
% 
% %% Clear temporary variables
% clearvars data raw;
% 
% % figure;
% % line(MDOD,LCOE,'Color','r');
% % ax1 = gca;
% % %ax1.YColor = 'r';
% % ax1_pos = ax1.Position; % position of first axes
% % ax2 = axes('Position',ax1_pos,...
% %     'XAxisLocation','top',...
% %     'YAxisLocation','right',...
% %     'Color','none','');
% % 
% % line(Cycles,Batt,'Parent',ax2);
% 
% 
% %line(Cycles,CapCost,'Parent',ax2,'Color','r');
% figure;
% [AX,H1,H2]=plotyy(MDOD,LCOE,MDOD,Batt);
% grid on;
% ylabel(AX(1),'LCOE ($/kWh)') % left y-axis
% ylabel(AX(2),'Battery Capacity (Ah)') % right y-axis
% %set(AX(2),'XAxisLocation','top','xlim',[Cycles(1) Cycles(end)]);
% set(AX,'fontsize',15);
% xlabel('Maximum Depth of Discharge (%)');

% subplot(2,1,2);
% [AX1,H1,H2]=plotyy(MDOD,LCOE,MDOD,CapCost);
% grid on;
% ylabel(AX1(1),'LCOE ($/kWh)') % left y-axis
% ylabel(AX1(2),'Capital Cost ($/W)') % right y-axis
% %set(AX(2),'XAxisLocation','top','xlim',[Cycles(1) Cycles(end)]);
% xlabel('Maximum Depth of Discharge (%)');


%% thesis plot 2 

[~, ~, raw] = xlsread('/Users/vmehra813/Dropbox (MIT)/[ulink]/Simulation/MatLab/uLink/Modeling Analysis/Generation + Reliability Model/thesis_lcoe_biz.xlsx','Sheet4');
raw = raw(2:end,:);

%Create output variable
data = reshape([raw{:}],size(raw));

%Allocate imported array to column variable names
MDOD1 = data(:,1);
Cycles1 = data(:,2);
LCOE1 = data(:,3);
CapCost1 = data(:,4);
Batt1 = data(:,5);
PV1 = data(:,6);
Years = data(:,7);
NumReplace = data(:,8);

%Clear temporary variables
clearvars data raw;

figure;
subplot(2,1,1);
[AX,H1,H2]=plotyy(MDOD1,LCOE1,MDOD1,Batt1);
grid on;
ylabel(AX(1),'LCOE ($/kWh)') % left y-axis
ylabel(AX(2),'Battery Capacity (Ah)') % right y-axis
%set(AX(2),'XAxisLocation','top','xlim',[Cycles(1) Cycles(end)]);
set(AX,'fontsize',15);
xlabel('Maximum Depth of Discharge (%)');

subplot(2,1,2);
[AX1,H1,H2]=plotyy(MDOD1,Years,MDOD1,NumReplace);
grid on;
ylabel(AX1(1),'Calendar Life (Years)') % left y-axis
ylabel(AX1(2),'Number of Replacements') % right y-axis
%set(AX(2),'XAxisLocation','top','xlim',[Cycles(1) Cycles(end)]);
set(AX1,'fontsize',15);
xlabel('Maximum Depth of Discharge (%)');



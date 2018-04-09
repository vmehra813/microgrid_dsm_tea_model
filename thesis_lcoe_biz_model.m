%% plots for thesis on financial sensitivity

%% Import the data
[~, ~, raw] = xlsread('/Users/vmehra813/Dropbox (MIT)/[ulink]/Simulation/MatLab/uLink/Modeling Analysis/Generation + Reliability Model/thesis_lcoe_biz.xlsx','Sheet5');
raw = raw(2:8,1:12);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Allocate imported array to column variable names
Conn = data(:,1);
Conn = Conn(1:5);
MP8 = data(:,2);
MP8 = MP8(1:5);

MP10 = data(:,3);
MP10 = MP10(1:5);

MP12 = data(:,4);
MP12 = MP12(1:5);

MP14 = data(:,5);
MP14 = MP14(1:5);

MP16 = data(:,6);
MP16 = MP16(1:5);

MP18 = data(:,7);
MP18 = MP18(1:5);

MP20 = data(:,8);
MP20 = MP20(1:5);

MP = [MP8, MP10, MP12, MP14, MP16, MP18, MP20];

%VarName9 = data(:,9);
DiscRate = data(:,10);
IRR = data(:,11);
LCOE = data(:,12);

%% Clear temporary variables
%clearvars data raw R;
figure;
a=contourf(Conn,DiscRate,MP');
 clabel(a,'manual','FontSize',15);
 set(gca,'fontsize',15);
 xlabel('Connection Fee ($)');
 ylabel('Discount Rate (%)');
 colorbar;
 legend('Monthly Payment ($)');
 
%  
% figure;
%  plot(DiscRate,100*IRR);
%  hold on;
%  plot(DiscRate,100*LCOE);
%  grid on;
%  hold off;
%  colorbar;
 
%% SECOND PLOT 
 
%  %% Import the data
[~, ~, raw] = xlsread('/Users/vmehra813/Dropbox (MIT)/[ulink]/Simulation/MatLab/uLink/Modeling Analysis/Generation + Reliability Model/thesis_lcoe_biz.xlsx','Sheet6');
raw = raw(1:5,2:16);

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
Payback = data(:,1);
IRR_8 = data(:,2);
IRR_10 = data(:,3);
IRR_12 = data(:,4);
IRR_14 = data(:,5);
IRR_16 = data(:,6);
IRR_18 = data(:,7);
IRR_20 = data(:,8);

IRR_P = 100*[IRR_8, IRR_10, IRR_12, IRR_14, IRR_16, IRR_18, IRR_20];


MP_8_PB = data(:,9);
MP_10_PB = data(:,10);
MP_12_PB = data(:,11);
MP_14_PB = data(:,12);
MP_16_PB = data(:,13);
MP_18_PB = data(:,14);
MP_20_PB = data(:,15);

MP_PB = [MP_8_PB, MP_10_PB, MP_12_PB, MP_14_PB, MP_16_PB, MP_18_PB, MP_20_PB];

%clearvars data raw;
 
figure;  % NO CONNECTION FEE - HOLD CONSTANT 
a=contourf(Payback,DiscRate,MP_PB');
 clabel(a,'manual','FontSize',15);
 set(gca,'fontsize',15);
 xlabel('Payback Period (Years)');
 ylabel('Discount Rate (%)');
 colorbar;
 legend('Monthly Payment ($)');
 
 figure;
a=contourf(Payback,DiscRate,IRR_P');
 clabel(a,'manual','FontSize',15);
 set(gca,'fontsize',15);
 xlabel('Payback Period (Years)');
 ylabel('Discount Rate (%)');
 colorbar;
 legend('IRR (%)');

 


function [available_solar] = singlehome_gen_sizing(config,solar,num,pv)

%% Description: Solar Irradiance Data for Individual SHS simulation (PV generation estimation)                                                              

%% Time Sampling
% parameters of year: sampling for generation simulation
x = 365;                                % number of days sampled (could be 365 to simulate entire year)
pDuration = 24*x;                       % number of hrs
nPeriods = 1;                           % number of time per year to take sample of pDuration

% Select hours of year for simulation                                  
increment = round(8760/nPeriods);
offsetmax = increment - pDuration;
yearFraction = pDuration*nPeriods/8760;
hindex = [];
for i = 1:nPeriods
    offset = round(offsetmax/2);        % place block of hours in middle of chunk
    hindex = [hindex;((i-1)*increment+offset+(1:pDuration)')];
end
%config.hindex = hindex;



%% Solar Output
%Based on Tata Power Solar platinum series, MNRE certified, along with the
%TS 250 Series
%a.available_solar = zeros(length(solar.irradiance_plane),num.A);


PV_eff = config.PVEff;   % 235 watt panels and 180 Ah is (10; 75 Watt and 75 Ah is (7)

PV_area = pv / 1000 * (1/PV_eff);
%PV_area = config.PV_area;     %235 watt panels and 180 Ah is (10); 75 Watt and 75 Ah is (7)
timestep = 1; %solar data is in one hour time increments

PV_output = [];
for i=1:num.A
% 
PV_output(:,i) = solar.irradiance_plane*PV_area*PV_eff*timestep;  % PV output generation, hourly
% 
% %fcn of temperature
end

% PV_eff_temp = zeros(length(solar.irradiance_plane),1);
% for j=1:num.A
%     for i=1:length(PV_eff_temp)
%     PV_eff_temp(i,j) = PV_eff*(1-solar.temp_coeff*(solar.tAmbient(i)-solar.temp_STC));
%     end
% end
% 
% PV_output_temp = zeros(size(a.available_solar));
% 
% % use to compuete output based on temp
% 
% for i=1:num.A
%     PV_output_temp(:,i) = solar.irradiance_plane.*PV_eff_temp*PV_area*timestep;
% end


solar.cc_fcn_gen= polyval(solar.cc_fcn,PV_output);
solar.cc_loss=(PV_output).*(100-solar.cc_fcn_gen)/100;


%available_solar = PV_output - solar.cc_loss; %available solar energy after charge controller in Wh
%a.available_solar_temp = PV_output_temp*CC_efficiency;  % take into account temperature

%a.available_solar_diff = a.available_solar - a.available_solar_temp;

available_solar = PV_output*.95; % 95% charge controller efficiency
% UNSTABLE VALUES ABOVE 300W 





%% Plotting 

% figure(1);
% subplot(2,2,1);
% plot(hindex, a.available_solar, 'g');
% title('Total Average Irradiance Over Patamda, Jharkhand');
% xlabel('Hour of Year');
% ylabel('Watts');
% subplot(2,2,2);
% plot(hindex,a.available_solar_temp,'r');
% title('Total Average Irradiance Over Patamda, Jharkhand (with Temp)');
% xlabel('Hour of Year');
% ylabel('Watts (Temp)');
% subplot(2,2,3);
% plot(hindex,a.available_solar_diff,'c');
% title('Difference With Taking Into Account Temperature');
% xlabel('Hour of Year');
% ylabel('Delta Watts');
% subplot(2,2,4);
% plot(solar.tAmbient);
% title('Temperature Throughout Year');
% xlabel('Hour of Year');
% ylabel('Temperature in C');


end


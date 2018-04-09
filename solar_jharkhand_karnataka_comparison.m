
addpath('/Users/vmehra813/Dropbox (MIT)/[ulink]/Simulation/MatLab/uLink/Modeling Analysis/TPS');


solar.outputs=xlsread('pvwatts_hourly_patamda_jharkhand.xlsx','DATA','A20:K8779');
%solar.irradiance_beam = solar.solaroutputs(:,4);    % plane of array irradiance (W/m^2) (beam = 4)
solar.irradiance_plane = solar.outputs(:,8);    % plane of array irradiance (W/m^2) (beam = 4)
solar.tAmbient=solar.outputs(:,6);             % ambient temperature in C
solar.hourofday=solar.outputs(:,3);
solar.month=solar.outputs(:,1);




solar.solaroutputs1=xlsread('pvwatts_hourly_karnataka.xlsx','DATA','A20:K8779');
%solar.irradiance_beam = solar.solaroutputs(:,4);    % plane of array irradiance (W/m^2) (beam = 4)
solar.irradiance_plane1 = solar.solaroutputs1(:,8);    % plane of array irradiance (W/m^2) (beam = 4)
solar.tAmbient1=solar.solaroutputs1(:,6);             % ambient temperature in C
solar.hourofday1=solar.solaroutputs1(:,3);
solar.month1=solar.solaroutputs1(:,1);


j_tot = sum(solar.irradiance_plane);
k_tot = sum(solar.irradiance_plane1);

j_day = j_tot / 365;
k_day = k_tot / 365;

%% 
s=1;
e=8760;

figure;
subplot(2,1,1),plot(s:e,solar.irradiance_plane1(s:e),'c');
set(gca,'fontsize',15);
xlabel('Hour of Year');
ylabel('Irradiance (W/m2)');
legend('Karnataka (13.25°N, 77.45°E)');
subplot(2,1,2),plot(s:e,solar.irradiance_plane(s:e),'m');
set(gca,'fontsize',15);
xlabel('Hour of Year');
ylabel('Irradiance (W/m2)');
legend('Jharkhand (22.95°N, 86.35°E)');


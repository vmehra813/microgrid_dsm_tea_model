
figure;
subplot(2,1,1), plot(total.Q_1(2:end))
hold on
plot(total.Q_2(2:end))
plot(total.Q_total(2:end))
legend('Q1','Q2','Q Total')
set(gca,'fontsize',14)
xlabel('Hour of Year')
ylabel('Capacity (Amps)')

subplot(2,1,2),
plot(total.Ebatt_max_in_lim1)
hold on
plot(total.Ebatt_max_in_limSOCmax)
plot(total.EoutMax_SOCminlimit)
plot(total.EoutMax_kineticmodellimit)
legend('KiBaM Charge Limit','SOC Charge Limit','KiBaM Discharge Limit','SOC Discharge Limit')
set(gca,'fontsize',14)
xlabel('Hour of Year')
ylabel('Constraint (Amps)')


% 
% figure;
% subplot(2,1,1), plot(critical.Q_1(2:end))
% hold on
% plot(critical.Q_2(2:end))
% plot(critical.Q_total(2:end))
% legend('Q1','Q2','Q Total')
% set(gca,'fontsize',14)
% xlabel('Hour of Year')
% ylabel('Capacity (Amps)')
% 
% subplot(2,1,2),
% plot(critical.Ebatt_max_in_lim1)
% hold on
% plot(critical.Ebatt_max_in_limSOCmax)
% plot(critical.EoutMax_SOCminlimit)
% plot(critical.EoutMax_kineticmodellimit)
% legend('KiBaM Charge Limit','SOC Charge Limit','KiBaM Discharge Limit','SOC Discharge Limit')
% set(gca,'fontsize',14)
% xlabel('Hour of Year')
% ylabel('Constraint (Amps)')

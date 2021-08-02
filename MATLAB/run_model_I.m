close all
clear
clc

t_end = 300;  % simulation end time
tRange = [0 t_end];
Y0 = [0.01; 0.01; 0.01; 0.01];  % initial densities for RS, RL, CS, CL

[tSol, YSol] = ode45(@model_I, tRange, Y0);

RS = YSol(:,1);
RL = YSol(:,2);
CS = YSol(:,3);
CL = YSol(:,4);

figure
subplot(2,1,1)
plot(tSol, RS, '-g', 'LineWidth', 1)
xlim([0 tRange(end)])
hold on
plot(tSol, RL, '-g', 'LineWidth', 2.5)
xlabel({'Time (day)'})
ylabel({'Density (mg/L)'})
legend('RS','RL')
hold off

subplot(2,1,2)
plot(tSol, CS, '-k', 'LineWidth', 1)
xlim([0 tRange(end)])
ylim([0 inf])
hold on
plot(tSol, CL, '-k', 'LineWidth', 2.5)
xlabel({'Time (day)'})
ylabel({'Density (mg/L)'})
legend('CS','CL')
hold off
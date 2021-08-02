close all
clear
clc

t_end = 300;  % simulation end time
tRange = [0 t_end];
Y0 = [0.01; 0.01; 0.01; 0.01; 0.01; 0.1];  % initial densities for RS, RL, JS, AS, JL, AL

[tSol, YSol] = ode45(@model_III, tRange, Y0);

RS = YSol(:,1);
RL = YSol(:,2);
JS = YSol(:,3);
AS = YSol(:,4);
JL = YSol(:,5);
AL = YSol(:,6);

figure
subplot(3,1,1)
plot(tSol, RS, '-g', 'LineWidth', 1)
xlim([0 tRange(end)])
hold on
plot(tSol, RL, '-g', 'LineWidth', 2.5)
xlabel({'Time (day)'})
ylabel({'Algal';'density (mg/L)'})
legend('RS','RL')
hold off

subplot(3,1,2)
plot(tSol, JS, '-k', 'LineWidth', 1)
xlim([0 tRange(end)])
ylim([0 inf])
hold on
plot(tSol, AS, '-k', 'LineWidth', 2.5)
xlabel({'Time (day)'})
ylabel({'Density (mg/L)'})
legend('JS','AS')
hold off

subplot(3,1,3)
plot(tSol, JL, '-k', 'LineWidth', 1)
xlim([0 tRange(end)])
ylim([0 inf])
hold on
plot(tSol, AL, '-k', 'LineWidth', 2.5)
xlabel({'Time (day)'})
ylabel({'Density (mg/L)'})
legend('JL','AL')
hold off
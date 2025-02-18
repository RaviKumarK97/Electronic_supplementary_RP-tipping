clc
clear
close all
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
figure('Color',[1 1 1],'renderer','painters')

% Birhythmicity
subplot(1,2,1)
[t_big,var_big] = ode45(@(t,var)vdP(var,[-0.001 1]),[0 50],[7;0],opts);
[t_small,var_small] = ode45(@(t,var)vdP(var,[-0.001 1]),[0 50],[0.1;0],opts);
[t_unstable,var_unstable] = ode45(@(t,var)vdP(var,[-0.001 1]),50:-0.01:0,[5;0],opts);

hold on
plot([0 0],[0 0],'ko','MarkerSize',7)
plot(var_big(2319:end,1),var_big(2319:end,2),'b','LineWidth',2.5)
plot(var_small(2319:end,1),var_small(2319:end,2),'r','LineWidth',2.5)
plot(var_unstable(end-700:end,1),var_unstable(end-700:end,2),'k--','LineWidth',2.5)

axis([-8 8 -12 12])
ylabel('y','FontSize',13,'FontWeight','bold','Rotation',0)
box on
axis square
set(gca,'FontSize',17,'FontWeight','normal','LineWidth',1.25,'TickDir','in');

%Time series
subplot(1,2,2)
hold on
plot(t_big,var_big(:,2),'b',t_small,var_small(:,2),'r',t_unstable,var_unstable(:,2),'k--','LineWidth',2);
xlim([0,30])
box on 
axis square 
set(gca,'FontSize',17,'FontWeight','normal','LineWidth',1.25,'TickDir','in');

% autonomous birhythmic vdp model
function [dvar] = vdP(var,par)
%parameters
a = 0.09355; b = 0.001942; d = par(1); mu =par(2);
x = var(1);y = var(2);
dx = y;
dy = mu*(1-x^2+a*x^4-b*x^6)*y-x-d*(y-x);
dvar = [dx;dy];
end



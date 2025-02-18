clc
clear
close all
figure('Color',[1 1 1],'renderer','painters')
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
hold on
% BIRHYTHMIC Glycolysis

subplot(1,2,1)
[t_big,var_big] = ode45(@(t,var)glycolysis(var,[0.29 1.15]),[0 1500],[45;5],opts);
[t_small,var_small] = ode45(@(t,var)glycolysis(var,[0.29 1.15]),[0 1500],[55;10],opts);
[t_unstable,var_unstable] = ode45(@(t,var)glycolysis(var,[0.29 1.15]),1500:-0.1:0,[55;10],opts);

v_0 = 0.29; si_0 = 1.15;
%equilibrium_point
syms x y
K=10;L=3600000;sm=10;n=5;q=1;ks=0.06;
eq1 = v_0+((si_0*y^n)/(K^n+y^n))-sm*((x*(1+x)*(1+y)^2)/(L+((1+x)*(1+y))^2)) ;
eq2 = q*sm*((x*(1+x)*(1+y)^2)/(L+((1+x)*(1+y))^2))-ks*y-((q*si_0*y^n)/(K^n+y^n));

eqns = [eq1 == 0, eq2 == 0];
eq_points = solve(eqns, [x, y]);
equilibm = [double(eq_points.x(2)) double(eq_points.y(2))];

hold on
plot(var_big(end-1000:end,1),var_big(end-1000:end,2),'b','LineWidth',2.5)
plot(var_small(end-600:end,1),var_small(end-600:end,2),'r','LineWidth',2.5)
plot(var_unstable(end-3000:end,1),var_unstable(end-3000:end,2),'k--','LineWidth',2.5)
plot([equilibm(1) equilibm(1)],[equilibm(2) equilibm(2)],'ko','MarkerFaceColor','k','MarkerSize',8)

xlim([20 100])
ylim([0 32])
xlabel('x','FontSize',13,'FontWeight','normal')
ylabel('y','FontSize',13,'FontWeight','normal','Rotation',0)
box on
axis square 
set(gca,'FontSize',18,'FontWeight','normal','GridAlpha',0.05,'LineWidth',1.25,'TickDir','none');

%Time series
subplot(1,2,2)
hold on
plot(t_big,var_big(:,2),'b',t_small,var_small(:,2),'r',t_unstable,var_unstable(:,2),'k--','LineWidth',2);
xlim([300,1500])
box on 
axis square 
set(gca,'FontSize',17,'FontWeight','normal','LineWidth',1.25,'TickDir','in');

% autonomous glycolysis model

function [VAR] = glycolysis(var,par)
x=var(1);y=var(2);
%parameters
v=par(1);si=par(2);K=10;L=3600000;sm=10;n=5;q=1;ks=0.06;
phi = (x*(1+x)*(1+y)^2)/(L+((1+x)*(1+y))^2);
dx=v+((si*y^n)/(K^n+y^n))-sm*phi;
dy=q*sm*phi-ks*y-((q*si*y^n)/(K^n+y^n));
VAR = [dx;dy];
end
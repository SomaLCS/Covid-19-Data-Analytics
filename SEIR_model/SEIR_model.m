clear all; clc;

%% Parameters
Pre_infec = 5.2;
f = 1/Pre_infec;

Duration = 7;
r=1/Duration;

R_0 = 2.2; % A single infected person will infect about 2.2 others in a totally susceptible population

N = 145.93e6; % Population of  Russia (2020)
beta = R_0/(N*Duration);


%% Differential equations
tspan = 0:1:365; % We will observe what happens over the next year
y0 = [N-28, 0, 25, 3, 0]; % Coronavirus count of Russia (Mar 12, 2020)

[t,y]=ode45(@(t,y) ode_fun_simple(t,y,beta), tspan, y0);
% ode45 is a built in Numerical Differential Equations solver in MATLAB.
% There are packages for this solver in other programming languages.




%% plot
plot(t,y,'LineWidth', 1.5, 'MarkerSize', 18)
legend('Susceptible','Pre-Infectious','Infectious','Recovered','Death', 'Location', 'Best')
xlabel('Days after March 12, 2020')
ylabel('Population')
title('Predicted Spread of COVID-19 in Russia')
grid on;
grid minor;
set(gca, 'FontSize', 26)
%saveas(gcf, 'corona_simple1.png')
clear all
close all 
clc

%% parameters
disp('Loading parameters')

N = 1000; % population rate
e = 0.5; % tramsmission probability
c = 1; % contact rate
L = 0.025; % lost immunity probability or rate
r = 0.275; % recovery probability rate

PLOTS = 0;
numsteps = 200; %% limit of x-axis


%% Initial conditions 
disp('Setting up initial conditions')

S0=N-1 ;   %%at least one person will be sick/infected at the beginning
I0=N-S0 ;
R0=0 ;

numRuns=50; %%changeable, it stands for sets of experiments done

disp('Running simulation')
for run=1:numRuns 
%% numerical solution
PLOTS = PLOTS + 1;

x(:,1)=[S0;I0;R0];
for k=1:numsteps
   S=x(1,k);
   I=x(2,k);
   R=x(3,k);

StoI=binornd(round(S),(e*c*I/N)); %%flip a coin or randomize the probability of state S to I
ItoR=binornd(I,r);
RtoS=binornd(R,L);

x(:,k+1)=x(:,k)+[-StoI+RtoS; StoI-ItoR; ItoR-RtoS];
end

%% plot results 
PLOTS
plot(0:numsteps,x(1,:),'b','MarkerSize',10,'LineWidth',2)
hold on
plot(0:numsteps,x(2,:),'r','MarkerSize',10,'LineWidth',2)
plot(0:numsteps,x(3,:),'g','MarkerSize',10,'LineWidth',2)

pause(0.1)

end
%% labelling plots
legend('Susceptible','Infected','Recovered')
xlabel('Number of Days')
ylabel('Population in Malaysia')
title('Modelling of Covid-19 Spreading Pattern in Malaysia')

disp('End of simulation')
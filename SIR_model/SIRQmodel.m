clear all
close all 
clc

%% parameters
disp('Loading parameters')

N = 1000; % population rate
e = 0.45; % tramsmission probability
c = 0.35; % contact rate
ft = 0.2; % testing rate
L = 0.025; % lost immunity probability or rate
r = 0.03; % recovery probability rate
v = 0.001; % dealth rate

PLOTS = 0;
numsteps = 400; %% limit of x-axis


%% Initial conditions 
disp('Setting up initial conditions')

S0=N-10 ;   %%at least one person will be sick/infected at the beginning
UQ0=N-S0 ;
Q0=0 ;
R0=0 ;
D0=0 ;

numRuns=10; %%changeable, it stands for sets of experiments done

disp('Running simulation')
for run=1:numRuns 
%% numerical solution
PLOTS = PLOTS + 1;

x(:,1)=[S0;UQ0;Q0;R0;D0];
for k=1:numsteps
   S=x(1,k);
   UQ=x(2,k);
   Q=x(3,k);
   R=x(4,k);
   D=x(5,k);

StoQ=binornd(round(S),(e*c*ft*UQ/N)); 
StoUQ=binornd(round(S),(e*c*(1-ft)*UQ/N)); 
UQtoR=binornd(UQ,r);
QtoR=binornd(Q,r);
UQtoD=binornd(UQ,v);
QtoD=binornd(Q,v);
RtoS=binornd(R,L);

x(:,k+1)=x(:,k)+[-StoQ-StoUQ+RtoS; StoUQ-UQtoR-UQtoD; StoQ-QtoR-QtoD; QtoR+UQtoR-RtoS; UQtoD+QtoD];
end

%% plot results 
PLOTS
plot(0:numsteps,x(1,:),'b','MarkerSize',10,'LineWidth',2)
hold on
plot(0:numsteps,x(2,:),'r','MarkerSize',10,'LineWidth',2)
plot(0:numsteps,x(3,:),'g','MarkerSize',10,'LineWidth',2)
plot(0:numsteps,x(4,:),'c','MarkerSize',10,'LineWidth',2)
plot(0:numsteps,x(5,:),'y','MarkerSize',10,'LineWidth',2)

pause(0.1)

end
%% Plot labelling
legend('Susceptible','Unquarantined','Quarantined','Recovered','Dealth')
xlabel('Number of Days')
ylabel('Population in Malaysia')
title('Modelling of Covid-19 Spreading Pattern in Malaysia')

disp('End of simulation')
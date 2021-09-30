clc; clear all; close all;

beta = 1.29;                  % transmission rate     
gamma = 0.71;                 % recovery rate of H
gamma1 = 0.33;                % recovery rate of I
lambda = 0.62;                % incubation rate
omega = 0.06;                 % lost immunity rate
q = 0.83;                     % hospitalised rate
d = 0.004;                    % dead rate
N = 1000;                     % population size
delta_t = 0.0001;            % time step
t_end = 30;                   % time duration
t_steps = t_end/delta_t;      % total time step
y0 = [N-50 0 50 0 0 0]';      % Initial condition [S0 E0 I0 H0 R0]
runs = 10;                    % number of runs
time = 0:delta_t:t_end;       % time vector

S = zeros(runs,t_steps+1);       % an array of S
E = zeros(runs,t_steps+1);       % an array of E
I = zeros(runs,t_steps+1);       % an array of I
H = zeros(runs,t_steps+1);       % an array of H
R = zeros(runs,t_steps+1);       % an array of R
D = zeros(runs,t_steps+1);       % an array of D

S(:,1) = y0(1);        % Inital number of S
E(:,1) = y0(2);        % Inital number of E
I(:,1) = y0(3);        % Inital number of I
H(:,1) = y0(4);        % Inital number of H
R(:,1) = y0(5);        % Inital number of R
D(:,1) = y0(6);        % Inital number of D

for k = 1:runs
for t = 1:t_steps
    p1 = beta*(I(k,t)+E(k,t))*S(k,t)/N*delta_t;   % probability of a new exposed
    p2 = lambda*E(k,t)*delta_t;                   % probability of a new infection
    p3 = q*I(k,t)*delta_t;                        % probability of a new hospitalised
    p4 = gamma1*I(k,t)*delta_t;                    % probability of a new recovery from I
    p5 = gamma*H(k,t)*delta_t;                    % probability of a new recovery from H
    p6 = omega*R(k,t)*delta_t;                    % probability of a new susceptible 
    p7 = d*I(k,t)*delta_t;                       % probability of a new dead from I
    p8 = d*H(k,t)*delta_t;                       % probability of a new dead from H
    
    r = rand;                          % random number from a uniform distribution of (0,1)
    
    % conditions for a new exposed
    if (r > 0) && (r <= p1)
        S(k,t+1) = S(k,t)-1;
        E(k,t+1) = E(k,t)+1;
        I(k,t+1) = I(k,t);
        H(k,t+1) = H(k,t);
        R(k,t+1) = R(k,t);
        D(k,t+1) = D(k,t);
        
    % conditions for a new infection
    elseif (r > p1) && (r <= p1+p2)
        S(k,t+1) = S(k,t);
        E(k,t+1) = E(k,t)-1;
        I(k,t+1) = I(k,t)+1;
        H(k,t+1) = H(k,t);
        R(k,t+1) = R(k,t);
        D(k,t+1) = D(k,t);
    
    % conditions for a new hospitalised
    elseif (r > p1+p2) && (r <= p1+p2+p3)
        S(k,t+1) = S(k,t);
        E(k,t+1) = E(k,t);
        I(k,t+1) = I(k,t)-1;
        H(k,t+1) = H(k,t)+1;
        R(k,t+1) = R(k,t);
        D(k,t+1) = D(k,t);
    
    % conditions for a new recovery from I
    elseif (r > p1+p2+p3) && (r <= p1+p2+p3+p4)
    	S(k,t+1) = S(k,t);
    	E(k,t+1) = E(k,t);
        I(k,t+1) = I(k,t)-1;
        H(k,t+1) = H(k,t);
        R(k,t+1) = R(k,t)+1;
        D(k,t+1) = D(k,t);
        
    % conditions for a new recovery from H
    elseif (r > p1+p2+p3+p4) && (r <= p1+p2+p3+p4+p5)
    	S(k,t+1) = S(k,t);
    	E(k,t+1) = E(k,t);
        I(k,t+1) = I(k,t);
        H(k,t+1) = H(k,t)-1;
        R(k,t+1) = R(k,t)+1;
        D(k,t+1) = D(k,t);
        
    % conditions for a new susceptible
    elseif (r > p1+p2+p3+p4+p5) && (r <= p1+p2+p3+p4+p5+p6)
    	S(k,t+1) = S(k,t)+1;
    	E(k,t+1) = E(k,t);
        I(k,t+1) = I(k,t);
        H(k,t+1) = H(k,t);
        R(k,t+1) = R(k,t)-1;
        D(k,t+1) = D(k,t);
        
    % conditions for a new dead from I
    elseif (r > p1+p2+p3+p4+p5+p6) && (r <= p1+p2+p3+p4+p5+p6+p7)
    	S(k,t+1) = S(k,t);
    	E(k,t+1) = E(k,t);
        I(k,t+1) = I(k,t)-1;
        H(k,t+1) = H(k,t);
        R(k,t+1) = R(k,t);
        D(k,t+1) = D(k,t)+1;
        
    % conditions for a new dead from H
    elseif (r > p1+p2+p3+p4+p5+p6+p7) && (r <= p1+p2+p3+p4+p5+p6+p7+p8)
    	S(k,t+1) = S(k,t);
    	E(k,t+1) = E(k,t);
        I(k,t+1) = I(k,t);
        H(k,t+1) = H(k,t)-1;
        R(k,t+1) = R(k,t);
        D(k,t+1) = D(k,t)+1;
        
    % conditions for no change
    elseif (r > p1+p2+p3+p4+p5+p6+p7+p8) && (r < 1)
    	S(k,t+1) = S(k,t);
    	E(k,t+1) = E(k,t);
        I(k,t+1) = I(k,t);
        H(k,t+1) = H(k,t);
        R(k,t+1) = R(k,t);
        D(k,t+1) = D(k,t);
    end
end

% stochastic simulation 
figure(1)
plot3(time,k*ones(size(S(k,:))),S(k,:),'Color','#0071BC')
hold on
plot3(time,k*ones(size(E(k,:))),E(k,:),'Color','#D85218')
plot3(time,k*ones(size(I(k,:))),I(k,:),'Color','#ECB01F')
plot3(time,k*ones(size(H(k,:))),H(k,:),'Color','#7D2E8D')
plot3(time,k*ones(size(R(k,:))),R(k,:),'Color','#76AB2F')
plot3(time,k*ones(size(D(k,:))),D(k,:),'Color','#4CBDED')
xlabel('Number of days')
ylabel('Number of trials')
zlabel('Population size')
legend('S(t)', 'E(t)', 'I(t)', 'H(t)', 'R(t)', 'D(t)')

end
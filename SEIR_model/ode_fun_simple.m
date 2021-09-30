function dydt = ode_fun_simple(t, y,beta)

Death = 0.034;  % Death rate of COVID-19 (March, 2020)
Pre_infec = 5.2;
f = 1/Pre_infec;
Duration = 14;
r = 1/Duration;
S = y(1);
E= y(2);
I = y(3);

dS = -beta*I.*S;
dE = beta*I.*S -  f.*E;
dI = f*E - r*I;
dR = r*(1-Death)*I;
dD = (Death)*r*I;

dydt = [dS; dE; dI; dR; dD];
end

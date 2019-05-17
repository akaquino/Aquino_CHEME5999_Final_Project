clear
clc
%global F1 F2 F3 F4 V Eo V_phase1 x1_o x2_o x3_o x1_phase1_ss x2_phase1_ss x3_phase1_ss x4_phase1_ss
F1 = 20; % inlet 1 [uL/hr]
F2 = 20; % inlet 2 [uL/hr] 
F3 = 20; % inlet 3 [uL/hr]
F4 = 60; % outlet [uL/hr] 
V= 40; %uL 
Eo = 10; %uM
V_phase1 = V + 0.1*V;
x1_o = 0.10/Eo; %uM 
x2_o = 0.50/Eo; %uM 
x3_o = 0.20/Eo ; %uM 
x1_EXP2 = 0.485;
x2_EXP2 = 0.2852;
x3_EXP2 = 0.9998;
x4_EXP2 = 0.3148;

n=3.926;
kcat = 2.2741*10^5; % uM/hr
Km1 = 1.02; %uM    
Km2 = 2*Km1; %uM
theta = ((Eo*x3_EXP2)^n)/(50*Eo*x3_EXP2 + (Eo*x3_EXP2)^n); 

%% Use simulated annealing to estimate parameters

%Generate initial conditions
    % Generate random number within bounds a --> b: (b-a).*rand(iterations,1) + a; 
    iterations=1000000;
    Km1_SA_0=(Eo-10^(-3)).*rand(iterations,1) + 10^(-3);
    Km2_SA_0=(Eo-10^(-3)).*rand(iterations,1) + 10^(-3); 
p_SA_0 = [Km1_SA_0,Km2_SA_0];    %Set values for initial parameters

value_initial= (x1_o.*F1)./(kcat.*V_phase1) -  (x1_EXP2.*F4)./(kcat.*V_phase1) - theta.*(x1_EXP2./(Km1_SA_0 + x1_EXP2)).*(x2_EXP2./(Km2_SA_0 + x2_EXP2)) ...
        +(x2_o.*F2)./(kcat.*V_phase1) -  (x2_EXP2.*F4)./(kcat.*V_phase1) - theta.*(x1_EXP2./(Km1_SA_0 + x1_EXP2)).*(x2_EXP2./(Km2_SA_0 + x2_EXP2)) ...
        +(x3_o.*F3)./(kcat.*V_phase1) - (x3_EXP2.*F4)./(kcat.*V_phase1) ...
        +(theta.*(x1_EXP2./(Km1_SA_0 + x1_EXP2)).*(x2_EXP2./(Km2_SA_0 + x2_EXP2)) - (x4_EXP2.*F4)./(kcat.*V_phase1));

    Km1_SA_f =((0.1+0.1).*rand(iterations,1) - 0.1).*Km1_SA_0;
    Km2_SA_f =((0.1+0.1).*rand(iterations,1) - 0.1).*Km2_SA_0;
p_SA_f = [Km1_SA_f,Km2_SA_f];   
value_final=(x1_o.*F1)./(kcat.*V_phase1) -  (x1_EXP2.*F4)./(kcat.*V_phase1) - theta.*(x1_EXP2./(Km1_SA_f + x1_EXP2)).*(x2_EXP2./(Km2_SA_f + x2_EXP2)) ...
        +(x2_o.*F2)./(kcat.*V_phase1) -  (x2_EXP2.*F4)./(kcat.*V_phase1) - theta.*(x1_EXP2./(Km1_SA_f + x1_EXP2)).*(x2_EXP2./(Km2_SA_f + x2_EXP2)) ...
        +(x3_o.*F3)./(kcat.*V_phase1) - (x3_EXP2.*F4)./(kcat.*V_phase1) ...
        +(theta.*(x1_EXP2./(Km1_SA_f + x1_EXP2)).*(x2_EXP2./(Km2_SA_f + x2_EXP2)) - (x4_EXP2.*F4)./(kcat.*V_phase1));

% Simulated annealing
    %T is probability of accepting worse move
T_initial = 1;
T_final = exp(-8);
T=T_initial;
n=0;

while (T>=T_final)
    T=T_initial*0.8^n;
    accept_move=rand;
    
    if abs(value_initial)>= abs(value_final) 
        value_initial = value_final;       
    elseif accept_move>T
        value_initial = value_final;
    end  
    n=n+1;
    %display(T)
    T=T_initial*0.8^n;
    p_SA_1 = ((0.1+0.1).*rand(1,1) - 0.1).*p_SA_0;
    value_final=(x1_o.*F1)./(kcat.*V_phase1) -  (x1_EXP2.*F4)./(kcat.*V_phase1) - theta.*(x1_EXP2./(Km1_SA_f + x1_EXP2)).*(x2_EXP2./(Km2_SA_f + x2_EXP2)) ...
        +(x2_o.*F2)./(kcat.*V_phase1) -  (x2_EXP2.*F4)./(kcat.*V_phase1) - theta.*(x1_EXP2./(Km1_SA_f + x1_EXP2)).*(x2_EXP2./(Km2_SA_f + x2_EXP2)) ...
        +(x3_o.*F3)./(kcat.*V_phase1) - (x3_EXP2.*F4)./(kcat.*V_phase1) ...
        +(theta.*(x1_EXP2./(Km1_SA_f + x1_EXP2)).*(x2_EXP2./(Km2_SA_f + x2_EXP2)) - (x4_EXP2.*F4)./(kcat.*V_phase1));
end

% Choose row with minimum "value_initial" value
T=table([Km1_SA_0,Km2_SA_0, value_initial]);
p_0=[Km1_SA_0, Km2_SA_0, value_initial];
[M,I]=min(abs(p_0));
Km1_SA=p_0(I(3),2);
Km2_SA=p_0(I(3),3);

T=table([Km1;Km1_SA],[Km2;Km2_SA],'RowNames',{'Actual values','simulated annealing values'})
T.Properties.VariableNames = {'Km1','Km2'}


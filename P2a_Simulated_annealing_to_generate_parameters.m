clear
clc
%global F1 F2 F3 F4 V Eo V_phase1 x1_o x2_o x3_o x1_phase1_ss x2_phase1_ss x3_phase1_ss x4_phase1_ss
F1 = 20; % inlet 1 [uL/hr]
F2 = 20; % inlet 2 [uL/hr] 
F3 = 20; % inlet 3 [uL/hr]
F4 = 60; % outlet [uL/hr] 
V= 40; %uL 
Eo = 1.88; %uM
V_phase1 = V + 0.1*V;
x1_o = 10/Eo; %uM 
x2_o = 50/Eo; %uM 
x3_o = 20/Eo ; %uM 
x1_phase1_ss = 5.048291569994046e-04;
x2_phase1_ss = 21.134303315886935;
x3_phase1_ss = 10.604298525209279;
x4_phase1_ss = 5.311844237586983;

a=4;
kcat = 3600*1.5; % uM/hr
Km1 = 1.02; %uM    
Km2 = 2*Km1; %uM
theta = ((Eo*x3_phase1_ss)^a)/(50*Eo*x3_phase1_ss + (Eo*x3_phase1_ss)^a); 

%% Use simulated annealing to estimate parameters

%Generate initial conditions
    % Generate random number within bounds a --> b: (b-a).*rand(iterations,1) + a; 
    iterations=1000000;
    kcat_SA_0=(10^(6)-10^(-3)).*rand(iterations,1) + 10^(-3);                                 
    Km1_SA_0=(Eo-10^(-3)).*rand(iterations,1) + 10^(-3);
    Km2_SA_0=(Eo-10^(-3)).*rand(iterations,1) + 10^(-3); 
    theta_SA_0=rand(iterations,1); 
p_SA_0 = [kcat_SA_0,Km1_SA_0,Km2_SA_0,theta_SA_0];    %Set values for initial parameters

value_initial= (x1_o.*F1)./(kcat_SA_0.*V_phase1) -  (x1_phase1_ss.*F4)./(kcat_SA_0.*V_phase1) - theta_SA_0.*(x1_phase1_ss./(Km1_SA_0 + x1_phase1_ss)).*(x2_phase1_ss./(Km2_SA_0 + x2_phase1_ss)) ...
        +(x2_o.*F2)./(kcat_SA_0.*V_phase1) -  (x2_phase1_ss.*F4)./(kcat_SA_0.*V_phase1) - theta_SA_0.*(x1_phase1_ss./(Km1_SA_0 + x1_phase1_ss)).*(x2_phase1_ss./(Km2_SA_0 + x2_phase1_ss)) ...
        +(x3_o.*F3)./(kcat_SA_0.*V_phase1) - (x3_phase1_ss.*F4)./(kcat_SA_0.*V_phase1) ...
        +(theta_SA_0.*(x1_phase1_ss./(Km1_SA_0 + x1_phase1_ss)).*(x2_phase1_ss./(Km2_SA_0 + x2_phase1_ss)) - (x4_phase1_ss.*F4)./(kcat_SA_0.*V_phase1));

    kcat_SA_f=((0.1+0.1).*rand(iterations,1) - 0.1).*kcat_SA_0;
    Km1_SA_f =((0.1+0.1).*rand(iterations,1) - 0.1).*Km1_SA_0;
    Km2_SA_f =((0.1+0.1).*rand(iterations,1) - 0.1).*Km2_SA_0;
    theta_SA_f =((0.1+0.1).*rand(iterations,1) - 0.1).*theta_SA_0;
p_SA_f = [kcat_SA_f,Km1_SA_f,Km2_SA_f,theta_SA_f];   
value_final=(x1_o.*F1)./(kcat_SA_f.*V_phase1) -  (x1_phase1_ss.*F4)./(kcat_SA_f.*V_phase1) - theta_SA_f.*(x1_phase1_ss./(Km1_SA_f + x1_phase1_ss)).*(x2_phase1_ss./(Km2_SA_f + x2_phase1_ss)) ...
        +(x2_o.*F2)./(kcat_SA_f.*V_phase1) -  (x2_phase1_ss.*F4)./(kcat_SA_f.*V_phase1) - theta_SA_f.*(x1_phase1_ss./(Km1_SA_f + x1_phase1_ss)).*(x2_phase1_ss./(Km2_SA_f + x2_phase1_ss)) ...
        +(x3_o.*F3)./(kcat_SA_f.*V_phase1) - (x3_phase1_ss.*F4)./(kcat_SA_f.*V_phase1) ...
        +(theta_SA_f.*(x1_phase1_ss./(Km1_SA_f + x1_phase1_ss)).*(x2_phase1_ss./(Km2_SA_f + x2_phase1_ss)) - (x4_phase1_ss.*F4)./(kcat_SA_f.*V_phase1));

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
    value_final=(x1_o.*F1)./(kcat_SA_f.*V_phase1) -  (x1_phase1_ss.*F4)./(kcat_SA_f.*V_phase1) - theta_SA_f.*(x1_phase1_ss./(Km1_SA_f + x1_phase1_ss)).*(x2_phase1_ss./(Km2_SA_f + x2_phase1_ss)) ...
        +(x2_o.*F2)./(kcat_SA_f.*V_phase1) -  (x2_phase1_ss.*F4)./(kcat_SA_f.*V_phase1) - theta_SA_f.*(x1_phase1_ss./(Km1_SA_f + x1_phase1_ss)).*(x2_phase1_ss./(Km2_SA_f + x2_phase1_ss)) ...
        +(x3_o.*F3)./(kcat_SA_f.*V_phase1) - (x3_phase1_ss.*F4)./(kcat_SA_f.*V_phase1) ...
        +(theta_SA_f.*(x1_phase1_ss./(Km1_SA_f + x1_phase1_ss)).*(x2_phase1_ss./(Km2_SA_f + x2_phase1_ss)) - (x4_phase1_ss.*F4)./(kcat_SA_f.*V_phase1));
end

% Choose row with minimum "value_initial" value
T=table([kcat_SA_0,Km1_SA_0,Km2_SA_0,theta_SA_0, value_initial]);
p_0=[kcat_SA_0,Km1_SA_0, Km2_SA_0, theta_SA_0,value_initial];
[M,I]=min(abs(p_0));
kcat_SA=p_0(I(5),1);
Km1_SA=p_0(I(5),2);
Km2_SA=p_0(I(5),3);
theta_SA=p_0(I(5),4);

%Solve for parameter, n, using theta_SA
syms a_SA
eqn = ((Eo.*x3_phase1_ss).^a_SA)./(50.*Eo.*x3_phase1_ss + (Eo.*x3_phase1_ss).^a_SA) == theta_SA;
sola_SA = solve(eqn,a_SA); 
n_SA = eval(sola_SA); 
%disp('n='); disp(n1);

T=table([kcat;kcat_SA],[Km1;Km1_SA],[Km2;Km2_SA],[a;n_SA],'RowNames',{'Actual values','simulated annealing values'})
T.Properties.VariableNames = {'kcat','Km1','Km2','n'}


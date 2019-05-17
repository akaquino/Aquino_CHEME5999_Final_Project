% Experiment 2: enzyme is in excess 

global  x1_o x2_o x3_o x4_o Eo kcat Km1 Km2 n theta

Eo = 10; %uM

% x1_o = initial [UG], x2_o = intial [P], x3_o = initial [Mn2+], x4_o = intial [GP]
x1_o = 0.10/Eo; %uM 
x2_o = 0.50/Eo; %uM 
x3_o = 0.20/Eo ; %uM 
x4_o  = 0/Eo; 

% set initial conditions
% x1 = UG (glycan),  x2 = P (target protein/polypeptide), x3 = Mn2+ (cofactor), x4 = GP (glycosylated product)   
x1 = 5/Eo ; x2 = 3/Eo; x3 = 10/Eo; x4 = 3/Eo; 

x0 = [x1 x2 x3 x4] ; 

% Set flowrates and volume 
F1 = 20; % inlet 1 [uL/hr]
F2 = 20; % inlet 2 [uL/hr] 
F3 = 20; % inlet 3 [uL/hr]
F4 = 20; % outlet [uL/hr] 
V = 40; %uL 

% Catalytic parameter set 
n_initial = 4; 
kcat_initial = 3600*1.5; % uM/hr
Km1_initial = 1.02; %uM    
Km2_initial = 2*Km1_initial; %uM

kcat =kcat_initial; 
Km1 = Km1_initial;
Km2 = 2*Km1_initial;
n =n_initial; 

%set time span
t_final = 24;
M = t_final*4;    %measurements taken within time span
t_step = t_final/M; 
tspan = 0:t_step:t_final ;  

% solve system of odes
[t,x] = ode45(@(t,x) dxdt(t,x,V,F1,F2,F3,F4),tspan,x0);

%save concentrations for each species over time
%x_initial_values=[x(1,1) x(1,2) x(1,3) x(1,4);
                  %x(2,1) x(2,2) x(2,3) x(2,4);
                  %x(3,1) x(3,2) x(3,3) x(3,4);
                  %x(end,1) x(end,2) x(end,3) x(end,4)];
                  
x_initial_values=[x(:,1) x(:,2) x(:,3) x(:,4)];

%Plot variation in each overtime

figure
plot(t,x(:,1),'r');
hold on
plot(t,x(:,2),'b');
hold on
plot(t,x(:,3),'g');
hold on
plot(t,x(:,4),'m');
hold on
title ('Experiment 2: Concentration vs. Time') 
xlabel('Dimensionless Time ');
ylabel('Dimensionless Concentration'); 
legend('LLOs (dUG/dt)', 'Protein(dP/dt)','Mn2+(dM/dt)','Glycoprotein (dGP/dt)'); 

 %% Sensitivity analysis S= (dx/x)/(dp/p) --> p=parameters, x=concentration values

%reset parameters and perturb kcat
kcat = kcat_initial*1.1;
Km1 = Km1_initial;
Km2 = 2*Km1_initial;
n = n_initial;

[t,x]= ode45(@(t,x) dxdt(t,x,V,F1,F2,F3,F4),tspan,x0);
x_values = [x(:,1) x(:,2) x(:,3) x(:,4)];
sensitivity_kcat = abs(((x_values-x_initial_values)./x_initial_values)./((kcat-kcat_initial)./kcat));

%reset parameters and perturb Km1
kcat = kcat_initial;
Km1 = Km1_initial*1.1;
Km2 = 2*Km1_initial;
n = n_initial;

[t,x]= ode45(@(t,x) dxdt(t,x,V,F1,F2,F3,F4),tspan,x0);
x_values = [x(:,1) x(:,2) x(:,3) x(:,4)];
sensitivity_Km1 = abs(((x_values-x_initial_values)./x_initial_values)./((Km1-Km1_initial)./Km1));
    
%reset parameters and perturb Km2
kcat = kcat_initial;
Km1 = Km1_initial;
Km2 = 2*Km1_initial*1.1;
n = n_initial;

[t,x]= ode45(@(t,x) dxdt(t,x,V,F1,F2,F3,F4),tspan,x0);
x_values = [x(:,1) x(:,2) x(:,3) x(:,4)];
sensitivity_Km2 = abs(((x_values-x_initial_values)./x_initial_values)./((Km2-Km2_initial)./Km2));

%reset parameters and perturb n
kcat = kcat_initial;
Km1 = Km1_initial;
Km2 = 2*Km1_initial;
n = n_initial*1.1;

[t,x]= ode45(@(t,x) dxdt(t,x,V,F1,F2,F3,F4),tspan,x0);
x_values = [x(:,1) x(:,2) x(:,3) x(:,4)];
sensitivity_n = abs(((x_values-x_initial_values)./x_initial_values)./((n-n_initial)./n));

%% plot sensitivity of glycoprotein concentration to each parameter
    
figure
plot(t,sensitivity_kcat(:,4),'-k');
hold on
plot(t,sensitivity_Km1(:,4),'--r');
hold on
plot(t,sensitivity_Km2(:,4),':g');
hold on
plot(t,sensitivity_n(:,4),'-.b');
hold on
title ('Experiment 2: excess enzyme') 
xlabel('Dimensionless Time ');
ylabel('Sensitivity'); 
legend('kcat', 'Km1','Km2','n'); 

%%
function func = dxdt(t,x,V,F1,F2,F3,F4)
global   x1_o x2_o x3_o n Eo kcat Km1 Km2 theta
% glycan (UG) = x(1) ; target polypeptide (P) = x(2); cofactor (Mn2+) = x(3); product (GP) = x(4) 
% x1_o = initial [UG], x2_o = intial [P], x3_o = initial [Mn2+], x4_o = intial [GP]

% enzyme activity dependency on cofactor w/ hill fxn  
theta = ((Eo*x(3))^n)/(50*Eo*x(3) + (Eo*x(3))^n); 

% define odes 
func1= (x1_o*F1)/(kcat*V) -  (x(1)/Eo*F4)/(kcat*V) - theta*(x(1)/Eo/(Km1 + x(1)/Eo))*(x(2)/Eo/(Km2 + x(2)/Eo)); % dUG/dt
func2= (x2_o*F2)/(kcat*V) -  (x(2)/Eo*F4)/(kcat*V) - theta*(x(1)/Eo/(Km1 + x(1)/Eo))*(x(2)/Eo/(Km2 + x(2)/Eo)); %dP/dt
func3 = (x3_o*F3)/(kcat*V) - (x(3)/Eo*F4)/(kcat*V) ; %dMn2+/dt
func4 = theta*(x(1)/Eo/(Km1 + x(1)/Eo))*(x(2)/Eo/(Km2 + x(2)/Eo)) - (x(4)/Eo*F4)/(kcat*V) ;%dGP/dt 

func = [func1; func2; func3; func4];
end     
     
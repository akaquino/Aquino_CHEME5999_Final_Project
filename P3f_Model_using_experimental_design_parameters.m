%% CHEME 5999 Simulated annealing to estimate parameters

clear all
close all


global  x1_o x2_o x3_o x4_o Eo kcat Km1 Km2 n n_EXP kcat_EXP Km1_EXP Km2_EXP 

Eo = 1.88; %uM

% x1_o = initial [UG], x2_o = intial [P], x3_o = initial [Mn2+], x4_o = intial [GP]
x1_o = 10/Eo; %uM 
x2_o = 50/Eo; %uM 
x3_o = 20/Eo ; %uM 
x4_o  = 0/Eo; 

% set initial conditions
% x1 = UG (glycan),  x2 = P (target protein/polypeptide), x3 = Mn2+ (cofactor), x4 = GP (glycosylated product)   
x1 = 5/Eo ; x2 = 3/Eo; x3 = 10/Eo; x4 = 3/Eo; 

x0 = [x1 x2 x3 x4] ; 

% Set flowrates and volume 
F1 = 20; % inlet 1 [uL/hr]
F2 = 20; % inlet 2 [uL/hr] 
F3 = 20; % inlet 3 [uL/hr]
F4 = 60; % outlet [uL/hr] 
V = 40; %uL 

% Catalytic parameter set 
n = 4; 
kcat = 3600*1.5; % uM/hr
Km1 = 1.02; %uM    
Km2 = 2*Km1; %uM

n_EXP = 3.926; 
kcat_EXP = 227410; % uM/hr
Km1_EXP = 9.6746; %uM    
Km2_EXP = -0.073039; %uM

%set time span
t_final = 24;
M = 10 ;
t_step = t_final/M; 
tspan = 0:t_step:t_final ; 

% solve system of odes
[t,x] = ode45(@(t,x) dxdt(t,x,V,F1,F2,F3,F4),tspan,x0);
[t,y] = ode45(@(t,y) dydt(t,y,V,F1,F2,F3,F4),tspan,x0);

%Plot variation in each overtime
figure 
plot(t,x(:,1),'r');
hold on
plot(t,x(:,2),'b');
hold on
plot(t,x(:,3),'g');
hold on
plot(t,x(:,4),'k');
hold on

plot(t,y(:,1),'--r');
hold on
plot(t,y(:,2),'--b');
hold on
plot(t,y(:,3),'--g');
hold on
plot(t,y(:,4),'--k');
hold on

title ('Model Simulation Parameters vs. Simulated Annealing Parameters') 
xlabel('Dimensionless Time ');
ylabel('Dimensionless Concentration'); 
legend('Glycan (dUG/dt), model', 'Polypeptide(dP/dt), model','Cofactor(dM/dt), model','Glycoprotein (dGP/dt), model','Glycan (dUG/dt), exp design', 'Polypeptide(dP/dt), exp design','Cofactor(dM/dt), exp design','Glycoprotein (dGP/dt), exp design'); 
 
%save steady state values for each species 
x1_orig_ss = x(end,1); x2_orig_ss = x(end,2); x3_orig_ss = x(end,3); x4_orig_ss = x(end,4); 

%%
function func = dxdt(t,x,V,F1,F2,F3,F4)
global   x1_o x2_o x3_o n Eo kcat Km1 Km2 
% glycan (UG) = x(1) ; target polypeptide (P) = x(2); cofactor (Mn2+) = x(3); product (GP) = x(4) 
% x1_o = initial [UG], x2_o = intial [P], x3_o = initial [Mn2+], x4_o = intial [GP]

% enzyme activity dependency on cofactor w/ hill fxn  
theta = ((Eo*x(3))^n)/(50*Eo*x(3) + (Eo*x(3))^n); 

% define odes 
func1= (x1_o*F1)/(kcat*V) -  (x(1)*F4)/(kcat*V) - theta*(x(1)/(Km1 + x(1)))*(x(2)/(Km2 + x(2))); % dUG/dt
func2= (x2_o*F2)/(kcat*V) -  (x(2)*F4)/(kcat*V) - theta*(x(1)/(Km1 + x(1)))*(x(2)/(Km2 + x(2))); %dP/dt
func3 = (x3_o*F3)/(kcat*V) - (x(3)*F4)/(kcat*V) ; %dMn2+/dt
func4 = theta*(x(1)/(Km1 + x(1)))*(x(2)/(Km2 + x(2))) - (x(4)*F4)/(kcat*V) ;%dGP/dt 

func = [func1; func2; func3; func4];
end 

function func = dydt(t,y,V,F1,F2,F3,F4)
global   x1_o x2_o x3_o n_EXP Eo kcat_EXP Km1_EXP Km2_EXP 

% enzyme activity dependency on cofactor w/ hill fxn  
theta = ((Eo*y(3))^n_EXP)/(50*Eo*y(3) + (Eo*y(3))^n_EXP); 

% define odes 
func1= (x1_o*F1)/(kcat_EXP*V) -  (y(1)*F4)/(kcat_EXP*V) - theta*(y(1)/(Km1_EXP + y(1)))*(y(2)/(Km2_EXP + y(2))); % dUG/dt
func2= (x2_o*F2)/(kcat_EXP*V) -  (y(2)*F4)/(kcat_EXP*V) - theta*(y(1)/(Km1_EXP + y(1)))*(y(2)/(Km2_EXP + y(2))); %dP/dt
func3 = (x3_o*F3)/(kcat_EXP*V) - (y(3)*F4)/(kcat_EXP*V) ; %dMn2+/dt
func4 = theta*(y(1)/(Km1_EXP + y(1)))*(y(2)/(Km2_EXP + y(2))) - (y(4)*F4)/(kcat_EXP*V) ;%dGP/dt 

func = [func1; func2; func3; func4];
end 
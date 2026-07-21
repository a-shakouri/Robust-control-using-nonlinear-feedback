clear all;clc

tic

%======Nominal System======%
%=====(continuous-time)====%

% TWIP robot parameters

mB = 45;          % kg (pendulum body mass)
mW = 2;           % kg (wheel mass)
l  = 0.135;       % m  (pendulum length)
r  = 8 * 0.0254;  % m  (wheel radius: 8 in -> meters)
d  = 0.6;         % m  (wheel separation)

I1 = 1.9;         % kg*m^2
I2 = 2.1;         % kg*m^2
I3 = 1.6;         % kg*m^2

K  = 0.04;        % kg*m^2
J  = 0.02;        % kg*m^2

g  = 9.81;        % m/s^2

% Intermediate terms

mu1 = (mB + 2*mW + 2*J/r^2)*(I2 + mB*l^2) - (mB^2)*(l^2);
mu2 = I3 + 2*K + 2*(mW + J/r^2)*d^2;

% Coefficients (reindexed: a1, a2, b1, b2, b3)

a1 = -mB^2 * g * l^2 / mu1;

a2 = ((mB + 2*mW + 2*J/r^2)*mB*g*l) / mu1;

b1 = ((I2 + mB*l^2)/r + mB*l) / mu1;

b2 = -((mB*l)/r + mB + 2*mW + 2*J/r^2) / mu1;

b3 = -(d/r) / mu2;


A0c=[0 1 0  0 0 0
     0 0 a1 0 0 0
     0 0 0  1 0 0
     0 0 a2 0 0 0
     0 0 0  0 0 1
     0 0 0  0 0 0];

B0c=[0   0
     b1  b1
     0   0
     b2  b2
     0   0
     b3 -b3];

[n,m]=size(B0c);%Extract the size
 
dt=0.05; %sample time

%=======Uncertainty=======%

A1c=zeros(n,n);A2c=zeros(n,n);A3c=zeros(n,n);
A4c=zeros(n,n);A5c=zeros(n,n);
B1c=zeros(n,m);B2c=zeros(n,m);B3c=zeros(n,m);
B4c=zeros(n,m);B5c=zeros(n,m);

A1c(2,3)=1;A2c(4,3)=1;B3c(2,1)=1;B3c(2,2)=1;
B4c(4,1)=1;B4c(4,2)=1;B5c(6,1)=1;B5c(6,2)=-1;

ell=5;

A_totalc=[A0c A1c A2c A3c A4c A5c];
B_totalc=[B0c B1c B2c B3c B4c B5c];

A_total=[eye(n) zeros(n,n*ell)]+dt*A_totalc; %Discretization
B_total=dt*B_totalc; %Discretization

thetamax=0.05*abs([a1;a2;b1;b2;b3]);
thetamin=-thetamax;

%=====Covering the set=====%

alpha=0.9; %Rate of decay

[K_total,P_total,~,~,q]=setcover(alpha,A_total,B_total,thetamin,thetamax); %Covering the set

%% %% %% %% %% %% %% %% %% %% %% %% %%
%===Trajectory Simulation===%

T=10; %total time interval

x(:,1)=0.01*[-2;0.5;1;3;0.2;7]; %Initial state
z(:,1)=[0;0]; %Initial state of the controller

[A,B]=sys_eval_box(A_total,B_total,thetamin+(thetamax-thetamin)/10); %choose a system

i=1;

for t=0:dt:T
    
    w=zeros(n,1);
    
    T0(i)=t;
        
    z(:,i+1)=f_fun(alpha,x(:,i),z(:,i),P_total);
    u(:,i)=g_fun(x(:,i),z(:,i),K_total,P_total);
    
    x(:,i+1)=A*x(:,i)+B*u(:,i)+w;
    
    i=i+1;
    
end

X1=x(1,1:end-1)';X2=x(2,1:end-1)';X3=x(3,1:end-1)';
X4=x(4,1:end-1)';X5=x(5,1:end-1)';X6=x(6,1:end-1)';
U1=u(1,1:end)';U2=u(2,1:end)';
Z1=z(1,1:end-1)';
Z2=z(2,1:end-1)';

subplot(3,1,1)
plot(T0,X1,'linewidth',1)
hold on
plot(T0,X2,'linewidth',1)
plot(T0,X3,'linewidth',1)
plot(T0,X4,'linewidth',1)
plot(T0,X5,'linewidth',1)
plot(T0,X6,'linewidth',1)
ylabel('x(t)')
subplot(3,1,2)
plot(T0,U1,'linewidth',1)
hold on
plot(T0,U2,'linewidth',1)
ylabel('u(t)')
subplot(3,1,3)
plot(T0,Z2,'linewidth',1)
ylabel('z_2(t)')
xlabel('t')
hold on

toc
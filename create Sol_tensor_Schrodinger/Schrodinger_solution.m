
clc;clear;

%% Parameters:

K=2^11; 
%dt=1e-2;
T=0.1;

%% Build Potential Matrices:
global V_cos

D = -K/2 : K/2-1;
dx = (2*pi*K^-1); 
x = dx.*D;

V_cos = diag(1-cos(x));
M = diag(2*ones(1,K)) + diag(-1*ones(1,K-1),1) + diag(-1*ones(1,K-1),-1);

%% Initial Data:

% Initial value:
    U0 = orth(rand(K,K));
    S0 = diag(10.^(flip(-K:-1))); 

    V0 = orth(rand(K,K));
    X0=U0*S0*V0';


H = @(Y)   0.5*(M*Y+Y*M)-V_cos*Y*V_cos.'; 


sum=0
t_dis=linspace(0,T,300);
Sol_tensor=zeros(K,K,length(t_dis));
Sol_tensor(:,:,1)=X0;
X_n=X0;
f=@(Y) H(Y);
for j=2:length(t_dis)
    tic
    X_n=matOdeSolver(X_n, f,  t_dis(j-1), t_dis(j));
    Sol_tensor(:,:,j)=X_n;
    j
    toc
end
save('Sol_tensor.mat','Sol_tensor');




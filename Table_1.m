% This file is part of the software library that reproduces 
% numerical experiments from the manuscript
%   Daniel Kressner, Hei Yin Lam: 
%   "Randomized low-rank approximation of parameter-dependent matrices".
%
% Generates Table 1 from the manuscript.
% This script requires Chebfun2, Tensor toolbox.
%
% Since the code of ACA is not public, we did not include the result in
% this sctipt.
%


clear;
clc;
rng(123)

%% create spatial discretization 
n0=70;


x=((kron(ones(n0, 1), [1:n0]')-0.5)/(n0))';
y=((kron([1:n0]', ones(n0, 1))-0.5)/(n0))';


%% obtain affine approximation
sigma=1;
c=@(d,theta) (1./(n0*n0))*sigma.^2.*exp(-d.^2./(2.*theta.^2)); 
rank_fun=18;
f=chebfun2( @(d,theta) (1./(n0*n0))*sigma.^2.*exp(-d.^2./(2.*theta.^2)),[0,sqrt(2),0.1,sqrt(2)],rank_fun);
[C,D,R] = cdr(f);
U=C;
V=D*R'


%% discretize Theta n=300
Nell=300
theta = linspace(0.1,sqrt(2), Nell);
theta_index = 1:Nell;

%% discretize test Theta
Nell_test=300
theta_test = linspace(0.1,sqrt(2), Nell_test);


%% evaulate at grid point
A=zeros(n0^2,n0^2,rank_fun);
for i=1:n0^2
    for j=1:n0^2
        A(i,j,:)=V(norm([x(i)-x(j),y(i)-y(j)]));
    end
    i
end
A_tensor=tensor(A);

%%time test 10 trials

e_GN=[]
time_GN=[]
time_GN_on=[]

e_rsvd=[]
time_rsvd=[]
time_rsvd_on=[]

e_aca=[]
time_aca=[]
time_aca_on=[]

rank=[10,20,30,40,50,60,70];
for r=rank
    
    time_temp=[];
    time_temp_on=[];
    Omega=normrnd(0,1,[n0^2,r]);
    Psi=normrnd(0,1,[n0^2,round(1.2*r)]);
    % GN offline phase
    for count=1:10
        tic
        A_Omega=zeros(n0^2,r,rank_fun);
        A_Psi=zeros(round(1.2*r),n0^2,rank_fun);

        for i=1:rank_fun
            A_Omega(:,:,i)=A(:,:,i)*Omega;
            A_Psi(:,:,i)=Psi'*A(:,:,i);
        end

        A_Omega=tensor(A_Omega);
        A_Psi=tensor(A_Psi);
        time_temp=[time_temp,toc]      
        % GN online phase       
        tic       
        for i= 1:Nell_test     
            GN_factors(ttv(A_Omega,U(theta_test(i))',3).data,ttv(A_Psi,U(theta_test(i))',3).data,Psi);
        end
        time_temp_on=[time_temp_on,toc]

    end
    time_GN=[time_GN,mean(time_temp)]
    time_GN_on=[time_GN_on,mean(time_temp_on)]


    %%HMT offline phase

    time_temp=[];
    time_temp_on=[];
    Omega_1=normrnd(0,1,[n0^2,r]);
    for count=1:10
        tic
        X_HMT=[];
        for i=1:rank_fun
            X_HMT=[X_HMT,A(:,:,i)*Omega_1];         
        end
        [Q,~] = qr(X_HMT,0);
        Y_Omega_1=zeros(r*rank_fun,r,rank_fun);
        Z_Omega_1=zeros(r*rank_fun,n0^2,rank_fun);
        for i=1:rank_fun
           Z_Omega_1(:,:,i)=Q'*A(:,:,i);
           Y_Omega_1(:,:,i)=Z_Omega_1(:,:,i)*Omega_1;
        end
        Z_Omega_1=tensor(Z_Omega_1);
        Y_Omega_1=tensor(Y_Omega_1);
        time_temp=[time_temp,toc]

        %HMT online phase
        tic
        for i= 1:Nell_test      
            RSVD_linear_factors(Q,Y_Omega_1,Z_Omega_1,U(theta_test(i)));
        end
        time_temp_on=[time_temp_on,toc]
    end
    time_rsvd_on=[time_rsvd_on,mean(time_temp_on)]
    time_rsvd=[time_rsvd,mean(time_temp)]
  
end
% This file is part of the software library that reproduces 
% numerical experiments from the manuscript
%   Daniel Kressner, Hei Yin Lam: 
%   "Randomized low-rank approximation of parameter-dependent matrices".
%
% Generates Figure 9 from the manuscript.
% This script requires Chebfun2, Tensor toolbox. We also used Parallel Computing Toolbox
% parfor to speed up the computation.
%
% Since the code of ACA is not public, we did not include the result in
% this sctipt.
%
% Warning: This script take quite long time to run.

clear;
clc;
rng(123)
data = csvread('gas_data_all.csv')
%% create spatial discretization 
n0=70;% 70 orignal

n0_2=n0*n0;
data=normalize(data,1);
x=data(randi([1 13910],1,n0_2),:)

clear data

%% obtain affine approximation
sigma=1;
c=@(d,theta) (1/n0_2).*exp(-d.^2./(2.*theta.^2)); 
rank_fun=18;% 18 orignal
f=chebfun2( @(d,theta) (1/n0_2)*sigma.^2.*exp(-d.^2./(2.*theta.^2)),[0,200,10,120],rank_fun);
[C,D,R] = cdr(f);
U=C;
V=D*R'


%% discretize Theta n=300
Nell=300
theta = linspace(10,120, Nell);
theta_index = 1:Nell;

%% discretize Theta for quardurature
theta_dis=linspace(10,120, 300);

% evaulate at grid point
A=zeros(n0^2,n0^2,rank_fun);
xx=zeros(n0^2,n0^2);
for i=1:n0^2
    fprintf('Creating A:')
    parfor j=1:n0^2
        A(i,j,:)=V(norm(x(i,:)-x(j,:)));
        %xy(i,j)=norm(x(i,:)-x(j,:));
        xx(i,j)=norm(x(i,:)-x(j,:));
    end
    i
end


e=[]
e_rsvd=[]
e_svd=[]
fun_val_svd=[];
rank=[20,50,100,150,200,260];

fun_val_svd=zeros(length(rank),length(theta_dis))
%%error of svd
parfor d=1:length(theta_dis)
    tic
        fprintf('Calculating error of svd')
        cov=cov_operator(theta_dis(d),n0,xx);
        fun_val_svd(:,d)=Error_SVD(cov,rank);
        [d toc]
end

for r=1:length(rank)
    e_svd=[e_svd,sqrt(trapz(theta_dis,fun_val_svd(r,:)))];
end




%% error of HMT and GN
e_rsvd_table=[];
e_table=[];
for count=1:20
        e_rsvd=[];
        e=[];
    for r=rank
 
        Omega=normrnd(0,1,[n0^2,r]);
        Psi=normrnd(0,1,[n0^2,round(1.2*r)]);
        Omega_1=normrnd(0,1,[n0^2,r]);
    
        
        A_Omega=zeros(n0^2,r,rank_fun);
        A_Psi=zeros(round(1.2*r),n0^2,rank_fun);  

        %A_Omega_1=zeros(n0^2,r,rank_fun);

        
        for i=1:rank_fun
            %GN offline
            A_Omega(:,:,i)=A(:,:,i)*Omega;
            A_Psi(:,:,i)=Psi'*A(:,:,i);
        end
        %HMT offline
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
    
        A_Omega=tensor(A_Omega);
        A_Psi=tensor(A_Psi);
        Z_Omega_1=tensor(Z_Omega_1);
        Y_Omega_1=tensor(Y_Omega_1);
    
        fun_val_rsvd=zeros(1,length(theta_dis));
        fun_val_GN=zeros(1,length(theta_dis));

        for d=1:length(theta_dis)
            cov=cov_operator(theta_dis(d),n0,xx);
            fun_val_GN(d)=(norm(cov-GN(ttv(A_Omega,U(theta_dis(d))',3).data,ttv(A_Psi,U(theta_dis(d))',3).data,Psi),"fro"))^2;
            fun_val_rsvd(d)=(norm(cov-RSVD_linear(Q,Y_Omega_1,Z_Omega_1,U(theta_dis(d))),"fro"))^2;                 
            fprintf('randdomized trial: %d, i=%d, rank: %d\n', count,d,r);
            
        end 
        e_rsvd=[e_rsvd,sqrt(trapz(theta_dis,fun_val_rsvd))];
        e=[e,sqrt(trapz(theta_dis,fun_val_GN))];
    end
    e_table=[e_table;e];
    e_rsvd_table=[e_rsvd_table;e_rsvd];
end



%plot
errorbar(rank,mean(e_table),mean(e_table)-min(e_table),max(e_table)-mean(e_table),'LineWidth',2)
hold on
errorbar(rank,mean(e_rsvd_table),mean(e_rsvd_table)-min(e_rsvd_table),max(e_rsvd_table)-mean(e_rsvd_table),'LineWidth',2)
hold on
plot(rank,e_svd,"-.",'LineWidth',2)
set(gca, 'YScale', 'log')
xlabel('rank')
ylabel('L^2 error')
legend('Time-dependent Generlized Nystrom method','Time-dependent HMT','ACA', 'Truncated SVD')



function C=cov_operator(theta,n0,xx)
        C=zeros(n0^2,n0^2);
        for i=1:n0^2
            C(:,i)=(1/(n0*n0)).*exp(-xx(:,i).^2./(2.*theta.^2));
        end
end

function err= Error_SVD(Y,rank)
    [U,S,V]=svd(Y,'econ');
    err=[];
    for i=1:length(rank)
        err=[err;norm(Y-U(:,1:rank(i))*S(1:rank(i),1:rank(i))*V(:,1:rank(i))',"fro").^2];
    end
end


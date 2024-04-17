% This file is part of the software library that reproduces 
% numerical experiments from the manuscript
%   Daniel Kressner, Hei Yin Lam: 
%   "Randomized low-rank approximation of parameter-dependent matrices".
%
% Generates Figure 5 and 6 from the manuscript.

clear;
clc;
rng(123)
%load matrix 
% the solution matrix can be created using the script in create Sol_tensor_Schrodinger folder
% (warning: the solution tensor is large) or can be downloaded in 
% https://drive.google.com/file/d/15LM4_RLhqmsWgP7i3BtZbtM5UxCbUyPs/view?usp=share_link.
load('Sol_tensor_schordinger.mat')
Y=@(t) Sol_tensor(:,:,t);
n_dis=1:300;
T=0.1
t_dis=linspace(0,T,300);
M=size(Sol_tensor,1)
N=size(Sol_tensor,2)

e_table=[];
e_rsvd_table=[];
e_svd=[];
e_svd_fun=[];
tic
 rank=[10,20,30,40,50];
 %rank=[5,15,25,35,45];


%calculate svd error
for j=1:length(n_dis)
    e_svd_fun=[e_svd_fun,Error_SVD(Y(j),rank)];
    j
end
for r=1:length(rank)
    e_svd=[e_svd,sqrt(trapz(t_dis,e_svd_fun(r,:)))];
end

for count=1:20
        e=[];
        e_rsvd=[];
    for r=rank
        tic
        Omega=normrnd(0,1,[N,r]);
        Theta=normrnd(0,1,[M,round(1.2*r)]);
     
        % error function
        Error_GN= @(t)  (norm(Y(t)-GN(Y(t)*Omega,Theta'*Y(t),Theta),"fro"))^2;
        % L2 error
  
        e=[e,sqrt(trapz(t_dis,arrayfun(Error_GN,n_dis)))];

        Omega=normrnd(0,1,[N,r]);
        Error_rSVD= @(t)  (norm(Y(t)-RSVD(Y(t),Y(t)*Omega),"fro"))^2;    
        e_rsvd=[e_rsvd,sqrt(trapz(t_dis,arrayfun(Error_rSVD,n_dis)))];
          
        [r count]
        toc
    end
    e_table=[e_table;e];
    e_rsvd_table=[e_rsvd_table;e_rsvd];
end


%% plot 
toc
errorbar(rank,mean(e_table),mean(e_table)-min(e_table),max(e_table)-mean(e_table),'LineWidth',2)
hold on
errorbar(rank,mean(e_rsvd_table),mean(e_rsvd_table)-min(e_rsvd_table),max(e_rsvd_table)-mean(e_rsvd_table),'LineWidth',2)
hold on
plot(rank,e_svd,'--','LineWidth',2)
set(gca, 'YScale', 'log')
xlabel('rank','Interpreter','latex')
ylabel('$L^2$ error','Interpreter','latex')
legend('Parameter-dependent generalized Nystr\"om','Parameter-dependent HMT', 'Truncated SVD','Interpreter','latex')


figure

subplot(1,3,1)
semilogy(1:50,svds(Y(1),50),'.-','LineWidth',2)
title('$t=0$','Interpreter','latex')
xlim([1 50])
subplot(1,3,2)
semilogy(svds(Y(150),50),'.-','LineWidth',2)
xlim([1 50])
title('$t=0.05$','Interpreter','latex')
subplot(1,3,3)
semilogy(svds(Y(300),50),'.-','LineWidth',2)
xlim([1 50])
title('$t=0.1$','Interpreter','latex')


function err= Error_SVD(Y,rank)
    [U,S,V]=svd(Y,'econ');
    err=[];
    for i=1:length(rank)
        err=[err;norm(Y-U(:,1:rank(i))*S(1:rank(i),1:rank(i))*V(:,1:rank(i))',"fro").^2];
    end
end


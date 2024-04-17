% This file is part of the software library that reproduces 
% numerical experiments from the manuscript
%   Daniel Kressner, Hei Yin Lam: 
%   "Randomized low-rank approximation of parameter-dependent matrices".
%
% Generates Figure 3 and 4 from the manuscript.
clear;
clc;
rng(123)

%%load the matrix, 
% the solution matrix can be created using the script in create Sol_tensor_cookie folder
% or can be downloaded in 
% https://drive.google.com/file/d/15LM4_RLhqmsWgP7i3BtZbtM5UxCbUyPs/view?usp=share_link .
load('Sol_tensor_cookie.mat')
Y=@(t) Sol_tensor(:,:,t);
T=0.9;
n_dis=1:300;
t_dis=linspace(0,T,300);
M=size(Sol_tensor,1)
N=size(Sol_tensor,2)
 rank=[4,8,12,16,20,22,24];
e_table=[]
e_rsvd_table=[]
e_svd=[];
e_svd_fun=[];

%calculate svd error
for j=1:length(n_dis)
    e_svd_fun=[e_svd_fun,Error_SVD(Y(j),rank)];
    j
end
for r=1:length(rank)
    e_svd=[e_svd,sqrt(trapz(t_dis,e_svd_fun(r,:)))];
end


%calculate randomized approximation error
tic
for count=1:20
        e=[];
        e_rsvd=[];       
       
        %rank=[40,50,70,80];
    for r=rank

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
    end
    e_table=[e_table;e];
    e_rsvd_table=[e_rsvd_table;e_rsvd];
end


%% plot 
toc
errorbar(rank,mean(e_table),mean(e_table)-min(e_table),max(e_table)-mean(e_table),'LineWidth',1)
hold on
errorbar(rank,mean(e_rsvd_table),mean(e_rsvd_table)-min(e_rsvd_table),max(e_rsvd_table)-mean(e_rsvd_table),'LineWidth',1)
hold on
plot(rank,e_svd,'--','LineWidth',1)
set(gca, 'YScale', 'log')
xlabel('rank','Interpreter','latex')
ylabel('$L^2$ error','Interpreter','latex')
legend('Parameter-dependent generalized Nystr\"om','Parameter-dependent HMT', 'Truncated SVD','Interpreter','latex')


figure

subplot(1,3,1)
semilogy(1:50,svds(Y(1),50),'.-','LineWidth',1)
title('$t=0$','Interpreter','latex')
xlim([1 50])
subplot(1,3,2)
semilogy(svds(Y(150),50),'.-','LineWidth',1)
xlim([1 50])
title('$t=0.45$','Interpreter','latex')
subplot(1,3,3)
semilogy(svds(Y(300),50),'.-','LineWidth',1)
xlim([1 50])
title('$t=0.9$','Interpreter','latex')

function err= Error_SVD(Y,rank)
    [U,S,V]=svd(Y,'econ');
    err=[];
    for i=1:length(rank)
        err=[err;norm(Y-U(:,1:rank(i))*S(1:rank(i),1:rank(i))*V(:,1:rank(i))',"fro").^2];
    end
end


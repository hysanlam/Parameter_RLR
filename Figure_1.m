% This file is part of the software library that reproduces 
% numerical experiments from the manuscript
%   Daniel Kressner, Hei Yin Lam: 
%   "Randomized low-rank approximation of parameter-dependent matrices".
%
% Generates Figure 1 from the manuscript.

rng(123)
N=100;

%create time dependent matrix
W1=normrnd(0,1,[N,N]);
W1=W1-0.5*(W1+W1');
W2=normrnd(0,1,[N,N]);
W2=W2-0.5*(W2+W2');
D=diag(2.^-(1:N));
Y= @(t) expm(t.*W1)*exp(t)*D*(expm(t.*W2));

%discretize t
t_dis=linspace(0,1,300);


e_table=[]
e_rsvd_table=[]
e_svd_table=[]
tic
for count=1:20
        e=[];
        e_rsvd=[];
        e_svd=[];
        rank=[10,15,20,25,30,35,40];
    for r=rank

        Omega=normrnd(0,1,[N,r]);
        Theta=normrnd(0,1,[N,round(1.2*r)]);
     
        % error function
        Error_GN= @(t)  (norm(Y(t)-GN(Y(t)*Omega,Theta'*Y(t),Theta),"fro"))^2;
        % L2 error
  
        e=[e,sqrt(trapz(t_dis,arrayfun(Error_GN,t_dis)))];

        Omega=normrnd(0,1,[N,r]);
        Error_rSVD= @(t)  (norm(Y(t)-RSVD(Y(t),Y(t)*Omega),"fro"))^2;    
        e_rsvd=[e_rsvd,sqrt(trapz(t_dis,arrayfun(Error_rSVD,t_dis)))];

        e_svd=[e_svd,sqrt(trapz(t_dis,arrayfun(@(t) Error_SVD(Y(t),r),t_dis)))];
        [r,count]
    end
    e_table=[e_table;e];
    e_rsvd_table=[e_rsvd_table;e_rsvd];
    e_svd_table=[e_svd_table;e_svd];
end


%% plot 
toc
errorbar(rank,mean(e_table),mean(e_table)-min(e_table),max(e_table)-mean(e_table),'LineWidth',2)
hold on
errorbar(rank,mean(e_rsvd_table),mean(e_rsvd_table)-min(e_rsvd_table),max(e_rsvd_table)-mean(e_rsvd_table),'LineWidth',2)
hold on
plot(rank,e_svd,'--','LineWidth',2)
set(gca, 'YScale', 'log')
xlabel('rank')
ylabel('L^2 error')
legend('Parameter-dependent Nystrom method','Parameter-dependent HMT', 'Truncated SVD')


function err= Error_SVD(Y,r)
    [U,S,V]=svd(Y,'econ');
    err=norm(Y-U(:,1:r)*S(1:r,1:r)*V(:,1:r)',"fro").^2;

end


% This file is part of the software library that reproduces 
% numerical experiments from the manuscript
%   Daniel Kressner, Hei Yin Lam: 
%   "Randomized low-rank approximation of parameter-dependent matrices".
%
% Generates Figure 7 from the manuscript.
% This script requires Chebfun2 cdr function.


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


xy=zeros(n0^2,n0^2);
for i=1:n0^2
    for j=1:n0^2
        xy(i,j)=norm([x(i)-x(j),y(i)-y(j)]);
    end
end
d_dis=reshape(xy,n0^4,1);
d_dis=unique(d_dis)

%d_dis=linspace(0,sqrt(2), 500);
theta1=[0.1,0.115,0.136,0.166,0.213,0.297,0.491,sqrt(2)];

%% max norm of cdr approximation for fixed theta
error_table=[];
rank_fun=[1:30];
for r_fun=rank_fun
     
    f=chebfun2( @(d,theta) (1./(n0*n0))*sigma.^2.*exp(-d.^2./(2.*theta.^2)),[0,sqrt(2),0.1,sqrt(2)],r_fun);
    [C,D,R] = cdr(f);
    error=[];
    for t=theta1
        fun_val=arrayfun(@(x) c(x,t),d_dis);
        fun_val_cdr=arrayfun(@(x) C(t,:)* D * R(x,:)',d_dis);

        error=[error;max(abs(fun_val-fun_val_cdr))];
        t
    end
    error_table=[error_table,error];
    r_fun
end
semilogy(rank_fun,error_table,'-o','LineWidth',1)
xlabel('$s$','Interpreter','latex')
ylabel('$L^\infty$ error','Interpreter','latex')
legend('$t=0.1$','$t=0.115$','$t=0.136$','$t=0.166$','$t=0.213$','$t=0.297$','$t=0.491$','$t=1.4142$','Interpreter','latex')

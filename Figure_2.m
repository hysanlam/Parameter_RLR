% This file is part of the software library that reproduces 
% numerical experiments from the manuscript
%   Daniel Kressner, Hei Yin Lam: 
%   "Randomized low-rank approximation of parameter-dependent matrices".
%
% Generates Figure 2 from the manuscript.


rng(123)
N=100;

%create time dependent matrix
W1=normrnd(0,1,[N,N]);
W1=W1-0.5*(W1+W1');
W2=normrnd(0,1,[N,N]);
W2=W2-0.5*(W2+W2');
D=diag(2.^-(1:N));
Y= @(t) expm(t.*W1)*exp(t)*D*(expm(t.*W2))';

%discretize t
t_dis=linspace(0,1,300);


e_table=[]
e_constant_table=[];
e_rsvd_table=[]
e_rsvd_constant_table=[]
e_svd_table=[]
tic
for count=1:20
        e=[];
        e_constant=[];
        e_rsvd=[];
        e_rsvd_constant=[];
        e_svd=[];
        rank=[10,20,30,40];
    for r=rank

        Omega=normrnd(0,1,[N,r]);
        Theta=normrnd(0,1,[N,round(1.2*r)]);
        % error of para GN
        Error_GN= @(t)  (norm(Y(t)-GN(Y(t)*Omega,Theta'*Y(t),Theta),"fro"))^2;
        % L2 error
        e=[e,sqrt(trapz(t_dis,arrayfun(Error_GN,t_dis)))];
        % error of constant GN
        e_constant=[e_constant,sqrt(trapz(t_dis,arrayfun(@(t) Error_GN_constant(Y(t),N,r),t_dis)))];
        
        % error of para HMT
        Omega=normrnd(0,1,[N,r]);
        Error_rSVD= @(t)  (norm(Y(t)-RSVD(Y(t),Y(t)*Omega),"fro"))^2;    
        e_rsvd=[e_rsvd,sqrt(trapz(t_dis,arrayfun(Error_rSVD,t_dis)))];
        % error of constant HMT
        e_rsvd_constant=[e_rsvd_constant,sqrt(trapz(t_dis,arrayfun(@(t) Error_rSVD_constant(Y(t),N,r),t_dis)))];


       % e_svd=[e_svd,sqrt(trapz(t_dis,arrayfun(@(t) Error_SVD(Y(t),r),t_dis)))];
        [r,count]
    end
    e_table=[e_table;e];
     e_constant_table=[e_constant_table;e_constant];
    e_rsvd_table=[e_rsvd_table;e_rsvd];
    e_rsvd_constant_table=[e_rsvd_constant_table;e_rsvd_constant];
   
end
toc

%% plot HMT
x1 = {e_rsvd_constant_table,  e_rsvd_table};

ax(1) = subplot(1,2,1);
%%set whisker larger such that it incluldes outliers
boxplotGroup(ax(1), x1,'Colors','bg','GroupType','betweenGroups','Whisker',100)
ax(1) = gca;
ax(1).YAxis.Scale ="log";
groupCenters = @(nGroups,nMembers,interGroupSpace) ...
    nGroups/2+.5 : nGroups+interGroupSpace : (nGroups+interGroupSpace)*nMembers-1;
x1CenterTicks = groupCenters(numel(x1), size(x1{1},2), 1);
set(ax(1),'XTick',x1CenterTicks,'XTickLabels',rank)
%change medium to mean
h = findobj(gca,'Tag','Median');
temp=mean(e_rsvd_constant_table)
temp_1=mean(e_rsvd_table)
set(h,{'YData'},{[temp_1(4),temp_1(4)];[temp_1(3),temp_1(3)];[temp_1(2),temp_1(2)];[temp_1(1),temp_1(1)];...
                    [temp(4),temp(4)];[temp(3),temp(3)];[temp(2),temp(2)];[temp(1),temp(1)]})

xlabel("rank",'Interpreter','latex')
ylabel('$L^2$ error','Interpreter','latex')
title("HMT",'Interpreter','latex')


%plot GN

x1 = {e_constant_table,  e_table};
ax(1) = subplot(1,2,2);
%%set whisker larger such that it incluldes outliers
boxplotGroup(ax(1), x1,'Colors','bg','GroupType','betweenGroups','Whisker',100)
ax(1) = gca;
ax(1).YAxis.Scale ="log";
groupCenters = @(nGroups,nMembers,interGroupSpace) ...
    nGroups/2+.5 : nGroups+interGroupSpace : (nGroups+interGroupSpace)*nMembers-1;
x1CenterTicks = groupCenters(numel(x1), size(x1{1},2), 1);
set(ax(1),'XTick',x1CenterTicks,'XTickLabels',rank)
%change medium to mean
h_2 = findobj(gca,'Tag','Median');
temp=mean(e_constant_table)
temp_1=mean(e_table)
set(h_2,{'YData'},{[temp_1(4),temp_1(4)];[temp_1(3),temp_1(3)];[temp_1(2),temp_1(2)];[temp_1(1),temp_1(1)];...
                    [temp(4),temp(4)];[temp(3),temp(3)];[temp(2),temp(2)];[temp(1),temp(1)]})

xlabel("rank",'Interpreter','latex')
ylabel('$L^2$ error','Interpreter','latex')
title("Generalized Nystrom method",'Interpreter','latex')



function err= Error_rSVD_constant(Y,N,r)
    S=normrnd(0,1,[N,r]);
    err=norm(Y-RSVD(Y,Y*S),"fro").^2;
end

function err= Error_GN_constant(Y,N,r)
    O=normrnd(0,1,[N,r]);
    T=normrnd(0,1,[N,round(1.2*r)]);
    err=norm(Y-GN(Y*O,T'*Y,T),"fro").^2;
end


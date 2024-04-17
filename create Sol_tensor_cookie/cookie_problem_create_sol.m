
rng(121)

X0 = csvread('Y0.csv');
A0 = csvread('A0.csv');
A1 = csvread('A1.csv');
b = csvread('b.csv');
C1 = csvread('C1.csv');

%H = @(t,Y)   -(A0*Y+A1*Y*C1)+b*ones(1,101)

T=0.9
%   reference solution  
X_n=X0;
sum=0
t_dis=linspace(0,T,300);
Sol_tensor=zeros(1580,101,length(t_dis));
Sol_tensor(:,:,1)=X0;

f=@(Y) H(Y,A0,A1,C1,b);

for j=2:length(t_dis)
    X_n=matOdeSolver(X_n, f,  t_dis(j-1), t_dis(j));
    Sol_tensor(:,:,j)=X_n;
    j
end
save('Sol_tensor.mat','Sol_tensor');

function Y = H(Y,A0,A1,C1,b)
    
   Y= -(A0*Y+A1*Y*C1)+b*ones(1,101);
end


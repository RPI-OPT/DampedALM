clear 
close all

%Reading and splitting data to use it in objective and constraints
%mainA and mainb represent total available features and labels
%A and b are data used for minimizing MSE to calculate Lstar
%Ap and Au are data for protected and unprotected groups

mainA = readmatrix("X_train_a9a.csv");
mainb = readmatrix("Y_train_a9a.csv");

[n,d] = size(mainA);
n1    = ceil(2*n/3);
n2    = n - n1;

A     = mainA(1:n1,:);
b     = mainb(1:n1,:);

remA  = mainA(n1+1:n,:);
remb  = mainb(n1+1:n,:);

%Column 72 has the gender feature, we use this to determine protected and
%unprotected groups
indp = find(remA(:,72) == 1);
indu = find(remA(:,72) == 0);

Ap    = remA(indp,:);
Au    = remA(indu,:);
%%
up   = 0.1;
down = -0.1;
lb   = down*ones(d,1);
ub   = up*ones(d,1);
[xmat, fval, exitflag, output] = quadprog((1/n)*(A'*A), -(1/n)*A'*b, [], [], [], [], lb, ub);

fprintf('%f \n ', fval + (0.5/n)*norm(b)^2)

[np,~] = size(Ap);
[nu,~] = size(Au);

sum1 = Ap(1,:)'*Ap(1,:);
sum2 = Au(2,:)'*Au(2,:);

for i = 2:np
    sum1 = sum1 +  Ap(i,:)'*Ap(i,:);
end
for i = 2:nu
    sum2 = sum2 +  Au(i,:)'*Au(i,:);
end

rho = 3/np * norm(sum1) + 3/nu * norm(sum2);



%%
Lstar    = fval + (0.5/n)*norm(b)^2;
maxiter  = 30000;
epsilon  = 1e-2;
L        = 10^2;
x        = ones(d,1);
lambda   = .1;
beta0     = .01; %tune this as needed
const_nu = 'True';
u0       = 1e-1;
u1       = 1e-2;
gammau   = 3;
gammad   = 5;
theta    = mean(A*xmat);
mu       = rho;
kappa    = 2.1*Lstar;
%%

out1      = DPALM(Au, Ap, A,b, maxiter, epsilon, L, x, lambda, beta0, rho, gammau, gammad, theta, mu,...
           u0, const_nu, kappa, Lstar, down, up);

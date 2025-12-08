clear 
close all

%Reading and splitting data to use it in objective and constraints
%mainA and mainb represent total available features and labels
%A and b are data used for minimizing MSE to calculate Lstar
%Ap and Au are data for protected and unprotected groups
compas = readmatrix("COMPAS.csv");
mainA  = compas(:, 2:end);
mainb  = compas(:,1);

[n,d] = size(mainA);
n1    = ceil(2*n/3);
n2    = n - n1;

A     = mainA(1:n1,:);
b     = mainb(1:n1,:);

remA  = mainA(n1+1:n,:);
remb  = mainb(n1+1:n,:);

indp = find(remA(:,5) == 1);
indu = find(remA(:,5) == 0);

Ap    = remA(indp,:);
Au    = remA(indu,:);

up   = .10;
down = -.10;
lb   = down*ones(d,1);
ub   = up*ones(d,1);
[xmat, fval] = quadprog((A'*A), -A'*b, [], [], [], [], lb, ub);

fprintf('%f \n ', (1/n)*fval + (0.5/n)*norm(b)^2)

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
Lstar    = (1/n)*fval + (0.5/n)*norm(b)^2;
maxiter  = 10000;
epsilon  = 1e-2;
L        = 10^2;
% x        = xmat;
% x      = xmat - .1;
x        = ones(d,1);
lambda   = 1;
beta0    = 1;
u0       = 1e-1;
u1       = 1e-3;
const_nu = 'True';
gammau   = 3;
gammad   = 5;
theta    = mean(A*xmat);
mu       = rho;
kappa    = 2*Lstar;
%%
out1     = DPALM(Au, Ap, A,b, maxiter, epsilon, L, x, lambda, beta0, rho, gammau, gammad, theta, mu,...
           u0, const_nu, kappa, Lstar, down, up);
%%
out2     = prox_linear(Au, Ap, A,b, maxiter, epsilon, L, x, rho, gammau, gammad, theta, mu,...
           u0, u1, const_nu, kappa, Lstar, down, up);
%%

% save resultcompas1 out1 out2

%%


colors = [0 0 0; 255 0 0; 0 0 255; 0 255 0; 255 0 255; 192 192 192; 128 128 105; 244 164 96; 8 46 84; 
    210 105 30; 61 89 171; 0 210 87; 0 255 127; 65 105 225;221 160 221; 227 23 13]/255;
lines   = {'-.' '--' ':' '-'};

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
semilogy(out1.ngrad,out1.primal, 'Marker','x','Color',colors(5,:), 'LineWidth',2)
hold on
semilogy(out2.ngrad,out2.primal,'Marker','p', 'Color',colors(8,:), 'LineWidth',2)
hold off
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to PF','FontSize',16)
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
print(fig,'figures/primal_error_compas.pdf','-dpdf')

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
semilogy(out1.ngrad,out1.dual, 'Marker','x','Color',colors(5,:), 'LineWidth',2)
hold on
semilogy(out2.ngrad,out2.dual,'Marker','p', 'Color',colors(8,:), 'LineWidth',2)
hold off
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to DF','FontSize',16)
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
print(fig,'figures/dual_error_compas.pdf','-dpdf')












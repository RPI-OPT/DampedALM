clear; close all

seednum = 100;
m = 10;
n = 1000;

Q = cell(1,m);
c = cell(1,m);
d = zeros(m,1);

rhoweak = 0.1;
for j = 1:m
%     [Q{j},~] = qr(randn(n));
%     Q{j} = Q{j}(:, 1:n-5);
%     v = rand(n-5,1)*5 + 1; % can remove +1 here...
%     Q{j} = Q{j}*diag(v)*Q{j}';
%     Q{j} = (Q{j} + Q{j}') / 2;

    c{j} = randn(n,1)*5;
    c{j} = c{j} / norm(c{j},2);
    d(j) = min(0, randn*2) - 0.1;

    [U,~] = qr(randn(n));
    v = max(0, randn(n,1)*5); % v ~ Gassian.
    Q{j} = U*diag(v)*U' - rhoweak*eye(n);
    Q{j} = (Q{j} + Q{j}') / 2; % symmetrization is important
end
A         = [randn(m, n-m), eye(m)]/n;
b         = rand(m,1)/n; % check if b in the box constraint
opts      = [];
opts.Q    = Q;
opts.c    = c;
opts.down = -5;
opts.up   = 5;
opts.A    = A;
opts.b    = b;
y         = zeros(m,1);
x_l       = -5 * ones(n,1); 
x_u       = 5 * ones(n,1); 
x         = quadprog(eye(n),randn(n,1), [], [], A, b, x_l, x_u);
maxiter   = 10000;
eigval = zeros(m,1);
for i = 1:m
%     eigval(i) = max(eig(Q{i}));
    eigval(i) = norm(Q{i},2);
end

%%
beta      = .010;
epsilon   = 10^(-2);
rho       = sqrt(m)*norm(eigval);
const_nu  = "True";
u0        =  0.001;
out2      = RobustNLLeastSquares_iPALM(opts,x,y,beta,epsilon,m,n,rho,maxiter, const_nu, u0);
%%
rho       = sqrt(m)*norm(eigval);
const_nu  = "True";
epsilon   = 10^(-2);
nu0       = 0.001;
out5      = RobustNLLeastSquares_proxLinear(opts,x,epsilon,m,n,rho,maxiter,const_nu,nu0);
%%
% save result_practice out2 out4 out5
%%
delta = 1;
out4 = out5;

colors = [0 0 0; 255 0 0; 0 0 255; 0 255 0; 255 0 255; 192 192 192; 128 128 105; 244 164 96; 8 46 84; 
    210 105 30; 61 89 171; 0 210 87; 0 255 127; 65 105 225;221 160 221; 227 23 13]/255;
lines   = {'-.' '--' ':' '-'};

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
[~,id] = max(out2.ngrad); 
semilogy(out2.ngrad(1:delta:id),out2.dual(1:delta:id),'Marker','x','Color',colors(5,:), 'LineWidth',2)
hold on
[~,id] = max(out4.ngrad); 
semilogy(out4.ngrad(1:delta:id),out4.dual(1:delta:id),lines{2},'Marker','p','Color',colors(8,:), 'LineWidth',2)
hold off
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to DF','FontSize',16, 'FontWeight','bold')
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
print(fig,'figures/dual_feasibility_nlls.pdf','-dpdf')

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
[~,id] = max(out2.ngrad); 
semilogy(out2.ngrad(1:delta:id),out2.primal(1:delta:id),'Marker','x','Color',colors(5,:), 'LineWidth',2)
hold on
[~,id] = max(out4.ngrad); 
semilogy(out4.ngrad(1:delta:id),out4.primal(1:delta:id),lines{2},'Marker','p','Color',colors(8,:), 'LineWidth',2)
hold off
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to PF','FontSize',16)
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
print(fig,'figures/Primal_feasibility_nlls.pdf','-dpdf')


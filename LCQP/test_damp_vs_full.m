%this is an implementation to check the convergence with a damped and full
%dual stepsize in our DPALM paper
clear; close all;
rng(20220927,"twister");

m            = 10;
n            = 100;
num_trial    = 10;
v0num       = 4;
tol          = 1e-3;
maxouteriter = 1e4;
out         = cell(v0num,num_trial);
totolgradnum = zeros(v0num, num_trial);

x_l = -5 * ones(n,1); % lower box. -5. 0 for aipp to work.
x_u = 5 * ones(n,1); % upper box. 5


for num = 1:num_trial

    [U,~] = qr(randn(n));
    v     = max(0, randn(n,1)*5); % v ~ Gaussian.
    rho   = 1;
    Q     = U*diag(v)*U' - rho*eye(n);
    Q     = (Q + Q') / 2; % symmetrization is important
    A     = [randn(m, n-m), eye(m)];
    b     = randn(m,1) + 0.1;
    I     = eye(n);
    r     = randn(n,1);
    y0    = rand(m,1);
    x0    = quadprog(eye(n),randn(n,1),[],[],A,b,x_l,x_u);
    %%
    QP     = [];
    QP.A   = A;
    QP.b   = b;
    QP.Q   = Q;
    QP.r   = r;
    QP.ell = -5;
    QP.u   = 5;

    beta0 = 0.1;

    mu = rho;


    QP.v0 = 0.1;
    t0 = tic;
    [~,~,out{1, num}] = iPALM_QP(QP,x0, y0,I,...
        beta0,rho,tol,maxouteriter,mu);
    time1 = toc(t0);
    out{1,num}.time = time1;
    out{1,num}.iter = sum(out{1,num}.gradnum~=0);
    totalgradnum(1,num) = max(out{1,num}.gradnum);

    QP.v0 = 1;
    t0 = tic;
    [~,~,out{2, num}] = iPALM_QP(QP,x0, y0,I,...
        beta0,rho,tol,maxouteriter,mu);
    time2 = toc(t0);
    out{2,num}.time = time2;
    out{2,num}.iter = sum(out{2,num}.gradnum~=0);
    totalgradnum(2,num) = max(out{2,num}.gradnum);

    QP.v0 = 10;
    t0 = tic;
    [~,~,out{3, num}] = iPALM_QP(QP,x0, y0,I,...
        beta0,rho,tol,maxouteriter,mu);
    time3 = toc(t0);
    out{3,num}.time = time3;
    out{3,num}.iter = sum(out{3,num}.gradnum~=0);
    totalgradnum(3,num) = max(out{3,num}.gradnum);


    QP.v0 = inf;
    t0 = tic;
    [~,~,out{4, num}] = iPALM_QP(QP,x0, y0,I,...
        beta0,rho,tol,maxouteriter,mu);
    time4 = toc(t0);
    out{4,num}.time = time4;
    out{4,num}.iter = sum(out{4,num}.gradnum~=0);
    totalgradnum(4,num) = max(out{4,num}.gradnum);

end

%%
% save result_damp_vs_full_DPALM out totalgradnum num_trial v0num

%%
fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
plot(totalgradnum(1,:), 'bs', 'MarkerSize', 10);
hold on
plot(totalgradnum(2,:), 'rp', 'MarkerSize', 10);
plot(totalgradnum(3,:), 'c*', 'MarkerSize', 10);
plot(totalgradnum(4,:), 'ko', 'MarkerSize', 10);
legend('v_0=0.1','v_0=1','v_0=10','v_0=+\infty','Location','best');
xlabel('instance');
ylabel('Number of Gradient');
set(gca,'fontsize',16,'fontweight','bold');
print(fig,'figures/damp_vs_full.pdf','-dpdf');

%%
fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
instance_num = 1;
iternum1 = out{1,instance_num}.iter;
semilogy(out{1,instance_num}.alphaval(1:iternum1),'b-.','LineWidth',2);
hold on
iternum = out{2,instance_num}.iter;
semilogy(out{2,instance_num}.alphaval(1:iternum),'r-','LineWidth',2);

semilogy(out{1,instance_num}.betaval(1:iternum1),'k--','LineWidth',2);
legend('\alpha_k when v_0=0.1','\alpha_k when v_0=1','\beta_k','Location','best');
xlabel('Outer iteration number k');
set(gca,'fontsize',16,'fontweight','bold');
print(fig,'figures/alpha_vs_beta_inst1.pdf','-dpdf');

%%
fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
instance_num = 9;
iternum1 = out{1,instance_num}.iter;
semilogy(out{1,instance_num}.alphaval(1:iternum1),'b-.','LineWidth',2);
hold on
iternum = out{2,instance_num}.iter;
semilogy(out{2,instance_num}.alphaval(1:iternum),'r-','LineWidth',2);

semilogy(out{1,instance_num}.betaval(1:iternum1),'k--','LineWidth',2);
legend('\alpha_k when v_0=0.1','\alpha_k when v_0=1','\beta_k','Location','best');
xlabel('Outer iteration number k');
set(gca,'fontsize',16,'fontweight','bold');
print(fig,'figures/alpha_vs_beta_inst9.pdf','-dpdf');






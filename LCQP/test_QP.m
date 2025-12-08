clear; close all;
rng(20220927,"twister");

m            = 10;
n            = 1000;
tol          = 1e-3;
maxouteriter = 1e4;

x_l = -5 * ones(n,1); % lower box. -5. 0
x_u = 5 * ones(n,1); % upper box. 5

rho    = 1; % you can change this weakly-proximal parameter to make the
%problem easier or harder

[U,~] = qr(randn(n));
v     = max(0, randn(n,1)*5); % v ~ Gaussian.
Q     = U*diag(v)*U' - rho*eye(n);
Q     = (Q + Q') / 2; % symmetrization is important
A     = [randn(m, n-m), eye(m)];
b     = randn(m,1) + 0.1;
I     = eye(n);
r     = randn(n,1);
y0    = rand(m,1);
x0    = quadprog(eye(n),randn(n,1),[],[],A,b,x_l,x_u);

QP     = [];
QP.A   = A;
QP.b   = b;
QP.Q   = Q;
QP.r   = r;
QP.ell = -5;
QP.u   = 5;



%you need to tune beta0 if you change rho, default value for 3 values of
%rho are given here
if rho == 10
    beta0 = 5;
elseif rho == 1
    beta0 = 0.0001;
elseif rho == 0.1
    beta0 = 0.001;
end

mu = rho;
t0 = tic;
[~,~,out_lcqp] = iPALM_QP(QP,x0, y0,I,...
              beta0,rho,tol,maxouteriter,mu);
time1 = toc(t0);
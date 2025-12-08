
clear; close all;

m            = 10;
n            = 100;
temp         = 1e-4;

[U,~]        = qr(randn(n));
v            = max(0, randn(n,1)*5); % v ~ Gaussian.
Q            = U*diag(v)*U' + temp*eye(n);
Q            = (Q + Q') / 2;
A         = [randn(m, n-m), eye(m)];
xfeas     = rand(n,1);
I         = eye(n);
b         = A*xfeas;
H         = Q;
c         = randn(n,1);
cl        = -5*ones(n,1);
cu        = 5*ones(n,1);

% Objective = @(x) 0.5 * x(:)' * H * x(:) + c' * x(:);
% 
% sol = ga(Objective,n, A, b, [], [], cl, cu);

[x_mat, fval_mat] = quadprog(H, c, [], [], A, b, cl, cu);

opt = ActiveSet(); %creating an instance of active set method
[x, f, iters] = opt.run(H, -c, A, b, cu, cl, xfeas);

[xproj, ~] = quadprog(2*eye(n),-2*x, [], [], A, b, cl, cu);
x_as = xproj;
realf_as = 0.5*xproj' * Q * xproj + xproj'*c;

%%

QP     = [];
QP.A   = A;
QP.b   = b;
QP.Q   = Q;
QP.r   = c;
QP.ell = -5;
QP.u   = 5;
beta0  = 1e-2;
rho    = abs(-1*temp);
epsilon = 1e-3;
maxouter = 1e3;
mu      = rho;
y       = zeros(size(A(:,1)));
[xinit, ~] = quadprog(eye(n),randn(n,1),[],[],A,b,cl,cu);

% running DPALM and comparing with sqp
[xdp,y,out] = iPALM_QP(QP,xinit, y,I, beta0,rho,epsilon,maxouter,mu);

[xproj, ~] = quadprog(2*eye(n),-2*xdp, [], [], QP.A, QP.b, cl, cu);
realf = 0.5*xproj' * QP.Q * xproj + xproj'*QP.r;
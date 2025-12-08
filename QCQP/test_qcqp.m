clear; close all % Test code for nc_QCQP
%% generate data
rand('seed',20230119); randn('seed',20230119);
n = 1000; % 1000, 100
m = 10; % 10

rho = 1; % you can change this weakly-proximal parameter to make the
%problem easier or harder

Q = cell(1,m+1);
c = cell(1,m+1);
d = zeros(m+1,1);

[U,~] = qr(randn(n));

v = max(0, randn(n,1)*5); % v ~ Gassian.
Q{m+1} = U*diag(v)*U' - rho*eye(n);
Q{m+1} = (Q{m+1} + Q{m+1}') / 2; % symmetrization is important
c{m+1} = randn(n,1);

fprintf('weak convexity constant is %f\n', min(eig(Q{m+1})));

%%
for j = 1:m
    [Q{j},~] = qr(randn(n));
    Q{j} = Q{j}(:, 1:n-5);
    v = rand(n-5,1)*5 + 1; % can remove +1 here...
    Q{j} = Q{j}*diag(v)*Q{j}';
    Q{j} = (Q{j} + Q{j}') / 2;
    c{j} = randn(n,1)*5;
    d(j) = min(0, randn*2) - 0.1;
end

lb = -5*ones(n,1); % -5, -1. can change.
ub = 5*ones(n,1); % 5, 1. can change.

maxsubit = 10000;
tol = 1e-3;

opts1          = [];
opts1.x0       = zeros(n,1);
opts1.z0      = zeros(m,1);
opts1.maxsubit = maxsubit;
opts1.maxit    = 10000;
opts1.beta0    = 0.0001; %default value for rho=0.1, 1, 10; tune for others
opts1.tol      = tol;
opts1.rho      = -2*min(eig(Q{m+1}));
opts1.gammau   = 1.5;
opts1.gammad   = 2;
opts1.L        = norm(Q{m+1});
down           = -5;
up             = 5;

Q0             = Q{m+1};
c0             = c{m+1};
Q1             = Q;
c1             = c;
Q1(end)         = []; 
c1(end)         = [];

t2 = tic;
[out_qcqp] = iPALM_qcqp(Q0, Q1, c0, c1,d,down, up, opts1);
time2 = toc(t2);
fprintf("iPALM time = %f\n", time2)
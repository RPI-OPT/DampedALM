clear; close all

seednum = 100;
m = 10;
n = 1000;

Q = cell(1,m);
c = cell(1,m);
d = zeros(m,1);

rhoweak = 0.1; % you can change this weakly-proximal parameter to make the
%problem easier or harder

for j = 1:m

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
    eigval(i) = norm(Q{i},2);
end

%%
beta      = .01; %tune this for rho that is different than 0.1, 1, 10
epsilon   = 10^(-2);
rho       = sqrt(m)*norm(eigval);
const_nu  = "True"; 
u0        =  0.001;
out2      = RobustNLLeastSquares_iPALM(opts,x,y,beta,epsilon,m,n,rho,maxiter, const_nu, u0);


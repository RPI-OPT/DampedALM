function [out] = iPALM_qcqp(Q0, Q, c0, c,d,down,up,opts)

%# Problem setting
% we are solving a quadratic problem with non convex quadratic objective and a
% quadratic convex constraint, i.e. min_x 0.5*x'Q_0x + c_0'x  s.t
% 0.5*x'Q_jx + c_j'x + d_j <= 0, x \in [down,up]

%Parameters
% x is the primal variable, z is dual variable
% beta0 is the initial parameter for beta = beta0*sqrt(i)
% rho is the proximal parameter
% down and up are the box parameters for x
% L is the initial guess for the Lipschitz constant of smooth part of
% subroblem's objective function


if isfield(opts,'beta0');     beta0    = opts.beta0;       else,  beta0    = 2*1e-4;           end
if isfield(opts,'x0');        x0       = opts.x0;          else,  x0       = randn(n,1);       end
if isfield(opts,'z0');        z0       = opts.z0;          else,  z0       = randn(m,1);       end
if isfield(opts,'tol');       tol      = opts.tol;         else,  tol      = 1e-3;             end
if isfield(opts,'maxit');     maxit    = opts.maxit;       else,  maxit    = 100;              end
if isfield(opts,'maxsubit');  maxsubit = opts.maxsubit;    else,  maxsubit = 10000;            end 
if isfield(opts,'L');         L        = opts.L;           else,  L        = norm(Q0);         end 
if isfield(opts,'rho');       rho      = opts.rho;         else,  rho      = -2*min(eig(Q0));  end
if isfield(opts,'gammau');    gammau   = opts.gammau;      else, gammau    = 3;                end
if isfield(opts,'gammad');    gammad   = opts.gammad;      else,  gammad   = 5;                end

%initializing parameters
x            = x0;
z            = z0;
v0           = 10;

primalviol   = zeros(maxit-1,1);
dualviol     = primalviol;
compslack    = primalviol;
numgrad      = primalviol;
objective    =  primalviol;

mu           = rho;
i            = 1;
kktviola     = 1;
ngrad        = 0;

while i < maxit && kktviola > tol

    %Set the values of beta and inner tolerance and solve the inner subproblem
    beta               = beta0*sqrt(i);
    eps_k              = min([tol/2,tol*sqrt(rho), sqrt(rho),tol/(sqrt(beta*i))]);

    %update the primal variable, x, using inner Nesterov solver
    [x,ngrad,subi,suberr] = inner(x, z, beta , Q0, Q, c0, c,d, rho,down,up,mu,L,eps_k,maxsubit,gammau,gammad,ngrad);

    %set the value of stepsize gamma and update the multiplier z
    g     = zeros(length(z),1);
    ineq  = g;
    for j = 1:length(z)
        g(j)      = 0.5*x'*Q{j}*x + c{j}'*x + d(j);
        ineq(j)   = max(0,0.5*x'*Q{j}*x + c{j}'*x + d(j));
    end
    v       = v0/sqrt(i);
    gamma   = min(beta,v/(norm(ineq)));
    z       = z + gamma*max(-z/beta,g);
    
    %check the primal, dual feasibility and complementary slackness
    [primal, dual, compli] = getouterviol(x,Q0,Q,c0,c,d,z,up,down);
    kktviola               = primal + dual  +compli;
    obj                    = 0.5*x'*Q0*x + c0'*x;
    if beta == gamma
        message = 1;
    else
        message = 0;
    end
    fprintf(['NumIter: %g, objfun: %g, dualres: %.5f, primalres: %.5f, ' ...
        'compslack: %0.7f, TotalSubIter: %g, Subproblem err: %g, beta: %g\n'],...
        i,obj,dual,primal,compli,subi,suberr, beta);
    
    %store the primal, dual infeasibility, complemenatry slacknes, and kkt violation
    primalviol(i) = primal;
    dualviol(i)   = dual;
    compslack(i)  = compli;
    numgrad(i)    = ngrad;
    objective(i)  = obj;
    
    %move on to the next outer iteration
    i      = i + 1;
end

%compute the required oututs
out       = [];
out.primalviol = primalviol;
out.dualviol   = dualviol;
out.compslack  = compslack;
out.objective  = objective;
out.gradnum    = numgrad;%number of gradient computations
end


%inner subproblem solver
function [x,ngrad,i,err] = inner(x, z, beta , Q0, Q, c0, c,d, rho,down,up,mu,L,eps_k,subit,gammau,gammad,ngrad)
        
%inner subproblem setting 

% f(x) is the smooth part of the proximal augmented lagrangian function
% \Psi(x) is the indicator function on the box constraint
%xk is the previous value of x

%initializing parameters
xk      = x;
v       = xk;
gradsum = zeros(length(xk),1);
err     = 1;
i       = 0;
Ak      = 0;

%calculating one iteration of the required variables for the line search
ak      = ((1+mu*Ak)+sqrt((1+mu*Ak)^2 + 2*L*(1+mu*Ak)*Ak))/L; % calaulating a using the given quadratic equation
p       = (Ak*x+ak*v)/(Ak + ak);
gradp   = getgrad(xk,p,c0,Q0, c,d, Q,z,rho,beta);
T       = getT(xk,p,gradp,rho,L,down,up);
gradT   = getgrad(xk,T,c0,Q0, c,d, Q,z,rho,beta);
phiprim = L*(p-T)+gradT-gradp;
ngrad   = ngrad + 2;

%inner iteration loop 
while i < subit && err > eps_k

    %while loop for line search
    while (dot(phiprim, p - T) < (1/L)*norm(phiprim,"fro")^2) 
        L       = L*gammau;
        ak      = ((1+mu*Ak)+sqrt((1+mu*Ak)^2 + 2*L*(1+mu*Ak)*Ak))/L;
        p       = (Ak*x+ak*v)/(Ak + ak);
        gradp   = getgrad(xk,p,c0,Q0, c,d, Q,z,rho,beta);
        T       = getT(xk,p,gradp,rho,L,down,up);
        gradT   = getgrad(xk,T,c0,Q0, c,d, Q,z,rho,beta);
        phiprim = L*(p-T)+gradT-gradp;
        ngrad   = ngrad + 2;
    end
    
    gradx   = getgrad(xk,x,c0,Q0, c,d, Q,z,rho,beta); %calculate gradient at x_k
    gradsum = gradsum + ak * gradx; %use this to calculate v_k
    Ak      = Ak + ak; % Update A using A_k+1 = A_k + a_k+1
    ngrad   = ngrad+1;
    
    %update v, L and x after linesearch is satisfied
    v       = getv(v,Ak,xk,rho,gradsum,down,up);
    L       = L/gammad;
    x       = T;
    
    if rem(i,5) == 0
     err     = getinnerviol(xk,x,Q0,Q,c0,c,d,beta,z,rho,up,down);
    end

    i = i + 1;
end
end

function prox = getprox(xk,down,up,alpha,w,rho)
n = length(xk);
prox = zeros(n,1);
precalx = (w+alpha*rho*xk)/(1+alpha*rho);

id1 = precalx < up & precalx > down;
prox(id1) = precalx(id1);

id2 = -down + w + alpha*rho*(-down+xk) <= 0;
prox(id2) = down;

id3 = -up + w + alpha*rho*(-up+xk) >= 0;
prox(id3) = up;
end

function grad = getgrad(xk,x,c0,Q0, c,d, Q,zk,rho,beta)
sum = 0;
m = length(zk);
for i  = 1:m
    term = Q{i}*x;
    sum  = sum + max(0,zk(i)+beta*(x'*(0.5*term + c{i}) + d(i)))*(term +c{i});
end
grad = Q0*x + c0  + sum + rho*(x-xk);
end

function innerviol = getinnerviol(xk,x,Q0,Q,c0,c,d,beta,z,rho,up,low)
rho = 2*rho;
grad = getgrad(xk,x,c0,Q0, c,d, Q,z,rho,beta);

xi = zeros(length(x),1);
id1 = x > low & x < up;
xi(id1) = 0;

id2 = x == low;
xi(id2) = min(0,-grad(id2));

id3 = x == up;
xi(id3) = max(0,-grad(id3));

innerviol = norm(grad+xi);
end


%calculate the function value at w using the proximal mapping
function T = getT(xk,p,gradp,rho,L,down,up)
alpha = 1/L;
w     = p - alpha*gradp;
T     = getprox(xk,down,up,alpha,w,rho);
end


%calculate v using the proximal mapping and v_k+1 = min_x \psi_k+1(x):=
%\psi_k(x) + a_k+1[f(x_k+1)+< \nabla f(x_{k+1}, x - x_{k+1})> + \Psi(x)]
function v = getv(v,Ak,xk,rho,gradsum,down,up)
alpha = Ak;
grad  = (v-xk)+ gradsum;
w     = v - grad;
v     = getprox(xk,down,up,alpha,w,rho);
end


function [primal, dual, compli] = getouterviol(x,Q0,Q,c0,c,d,z,up,down)

sum1 = 0;
sum3 = sum1;
sum4 = zeros(length(z),1);
for j = 1:length(z)
    sum1    = sum1 + z(j)*(Q{j}*x + c{j});
    sum3    = sum3 + abs(z(j)*(0.5*x'*Q{j}*x+c{j}'*x + d(j)));
    sum4(j) = max(0,0.5*x'*Q{j}*x+c{j}'*x + d(j));
end
grad = Q0*x + c0 + sum1 ;
xi   = zeros(length(x),1);

id1 = x > down & x < up;
xi(id1) = 0;

id2 = x == down;
xi(id2) = min(0,-grad(id2));

id3 = x == up;
xi(id3) = max(0,-grad(id3));

dual    = norm(grad+xi);
primal  = norm(sum4);
compli  = sum3;

end


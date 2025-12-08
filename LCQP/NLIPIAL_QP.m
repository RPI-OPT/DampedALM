function [z,out] = NLIPIAL_QP(opts, z0,p0,I,mf, sigma, beta0, epsilon,maxiter, lambda,tol)
                    
%# Problem setting
% we are solving a quadratic problem with non convex quadratic objective and a
% convex constraint, i.e. min_x 0.5*x'Q_0x + c_0'x  s.t Ax = b, x \in [down,up]

%Parameters
% z is the primal variable, p is the dual variable
% beta0 is the initial parameter for beta = beta0*sqrt(i)
% rho is the proximal parameter
% low and up are the box parameters for x
% L is the Lipschitz constant of f

if isfield(opts,'Q');    Q0   = opts.Q;    else,  Q0   = randn(n,n);  end
if isfield(opts,'r');    c0   = opts.r;    else,  c0   = randn(n,1);  end
if isfield(opts,'A');    A    = opts.A;    else,  A    = randn(m,n);  end
if isfield(opts,'b');    b    = opts.b;    else,  b    = randn(m,1);  end
if isfield(opts,'ell');  down = opts.ell;  else,  down = -5;         end
if isfield(opts,'u');    up   = opts.u;    else,  up   = 5;          end 


%Initialization of variables 
z      = z0;%primal variable
p      = p0;%dual variable 
beta   = beta0;
khat   = 0;
k      = 1;
Csigma = 2*(1+2*sigma)^2/(1-sigma^2);
eps_k  = 10^(-4);
err    = 1;
rho    = 0.5/lambda;%strong convexity constant of the non-smooth part of the subproblem's objective
ngrad  = 0;
mu     = rho;

B      = A'*A;
L1     = norm(Q0  + 2*(rho)*I);
L2     = norm(B);

%initializing zero vectors to store the values later
primalviol = zeros(maxiter-1,1);
dualviol   = primalviol;
numgrad    = primalviol;
objective  =  primalviol;
betaval    = primalviol;
while k < maxiter && err > epsilon

    L    = L1 + beta*L2;
    [z,ngrad,i,suberr, vk] = innerNesterov(z,z,A,B,b,p, beta,rho, Q0,c0,down,up,mu,L,ngrad,eps_k,I);
    pold = p;
    p    = p + beta * (A*z-b); 

    if k == khat + 1
        zkhatplus = z;
        pkhat     = pold;
     end
     %update beta only if inner condition is satisfied
     if k > khat + 1
        
        delta = (1/(k-khat-1))*(getobj(Q0, c0, A, b, zkhatplus, pkhat,beta) ...
                - getobj(Q0, c0, A, b, z, p,beta) - 0.5*norm(p)^2/beta);
        rhs   = 0.5 * lambda * tol^2/Csigma;
        if delta <= rhs && beta < beta0*2^6
            beta = 2 * beta;
            khat = k;
        end

     end
     k = k + 1;

    [what, qhat]  = getouterviol(z,Q0,c0,A,b,p,up,down);
    err           = sqrt(what^2 + qhat^2);
    primalviol(k) = what;
    dualviol(k)   = qhat;
    objective(k)  = 0.5*z'*Q0*z + c0'*z;
    numgrad(k)    = ngrad;
    betaval(k)    = beta;
    
    if err <= epsilon
         fprintf(['NumIter: %d, objfun: %.5f, dualres: %.5f, primalres: %.5f,' ...
             'subiter : %d, suberr = %0.5f, Beta: = %0.5f\n'],k,objective(k),qhat,what,i, suberr, beta);
    end
        
end
%compute the required oututs
out            = [];
out.primalviol = primalviol;
out.dualviol   = dualviol;
out.objective  = objective;
out.gradnum    = numgrad;%number of gradient computations
out.betaval    = betaval;
end

%inner subproblem solver
function [x,ngrad,i,err,xhat] = innerNesterov(x,zk,A,B,b,yk, beta,rho, Q0,c0,down,up,mu,L,ngrad,eps,I)
%inner subproblem setting 

% f(x) is the smooth part of the objective, f(x) = 0.5*x'*Q0*x + c0'*x +
% y'(A*x-b) + 0.5*beta*norm(A*x-b)^2 + 0.5*rho*norm(x-xk)^2
% \Psi(x) = \tau_[-10,10](x) + 0.5*norm(x-xk)^2
%xk is the previous value of x
%a_i and A_i are parameters where A_{i+1} = A_{i}+a_{i+1}
%v_i is argmin \psi_i(x) where psi_i(x) is an estimate function updated as
%\psi_(i+1)(x) = \psi_i(x) + a_{i+1}(x)[f(x_{i+1})+<\nabla f(x_{i+1}),x-x_{i+1}>+\Psi(x)]
% T_L(y) is given by prox_{1/L \Psi(x)}(y - (1/L)*\nabla f(y)) 
    xhat     = x;
    xk       = x;
    alpha    = sqrt(mu/L);
    q        = mu/L;
    i        = 0;
    err      = 1;
    step     = 1/L;
    halfgrad = c0 + A'*yk - beta*A'*b - 2*rho*zk;
    C        = Q0+beta*B+2*rho*I;
    while i < 10000 && err > eps
        i        = i + 1;
        grad     = halfgrad + C*xhat;
        w        = xhat - step*grad;
        xold     = x;
        x        = min(up,max(down,w));
        alphaold = alpha;
        alpha    = 0.5*(q-alpha^2 +sqrt((q-alpha^2)^2 + 4*alpha^2));
        xhat     = x + ((alphaold*(1-alphaold))/(alphaold^2+alpha))*(x-xold);

        if rem(i,10) == 0
            err      = getinnerviol(halfgrad,xk,x,C,rho,up,down);
        end
        ngrad    = ngrad + 1;
    end
end


%calculating the KKT violation of the original problem
function [prim, dual] = getouterviol(x,Q0,c0,A,b,y,up,down)

grad = Q0*x + c0 + A'*y ;
xi = zeros(length(x),1);
id1 = x > down & x < up;
xi(id1) = 0;

id2 = x == down;
xi(id2) = min(0,-grad(id2));

id3 = x == up;
xi(id3) = max(0,-grad(id3));

dual   = norm(xi+grad) ;
prim   = norm(A*x-b);
end
   
function obj = getobj(Q0, c0, A, b, z, p,beta)
obj = 0.5*z'*Q0*z + c0'*z + 0.5*beta*norm(A*z-b)^2 + p'*(A*z-b);
end

% finding the error in the first order optimality condition for the inner
% subproblem
function innerviol = getinnerviol(halfgrad,xk,x,C,rho,up,down)
% grad = Q0*x + c0 + A'*y + beta*A'*(A*x-b) + 2*rho*(x-xk); 
% halfgrad = c0 + A'*yk - beta*A'*b - rho*zk;
grad = halfgrad + C*x;
xi   = zeros(length(x),1);

id1 = x > down & x < up;
xi(id1) = 0;

id2 = x == down;
xi(id2) = min(0,-grad(id2));

id3 = x == up;
xi(id3) = max(0,-grad(id3));

innerviol = norm(grad+xi);
end

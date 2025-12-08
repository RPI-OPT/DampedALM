function [x,y,out] = iPALM_QP(opts,x, y,I, beta0,rho,epsilon,maxouter,mu)

%# Problem setting
% we are solving a quadratic problem with non convex quadratic objective and a
% convex constraint, i.e. min_x 0.5*x'Q_0x + c_0'x  s.t Ax = b, x \in [down,up]

%Parameters
% x is the primal variable, y and z are dual variables
% beta0 is the initial parameter for beta = beta0*sqrt(i)
% rho is the proximal parameter
% low and up are the box parameters for x
% L is the Lipschitz constant of f

if isfield(opts, 'v0')   v0   = opts.v0;   else   v0   = 200;         end

if isfield(opts,'Q');    Q0   = opts.Q;    else,  Q0   = randn(n,n);  end
if isfield(opts,'r');    c0   = opts.r;    else,  c0   = randn(n,1);  end
if isfield(opts,'A');    A    = opts.A;    else,  A    = randn(m,n);  end
if isfield(opts,'b');    b    = opts.b;    else,  b    = randn(m,1);  end
if isfield(opts,'ell');  down = opts.ell;  else,  down = -5;         end
if isfield(opts,'u');    up   = opts.u;    else,  up   = 5;          end 

maxouteriter = maxouter+1;
i            = 1;
kkt          = 1;
ngrad        = 0;
B            = A'*A;
L1           = norm(Q0  + 2*(rho)*I);
L2           = norm(B);
%initializing zero vectors to store the values later
primalviol   = zeros(maxouteriter-1,1);
dualviol     = primalviol;
numgrad      = primalviol;
objective    =  primalviol;
betaval      = primalviol;
alphaval     = primalviol;

hist_id = [];
while i < maxouteriter && kkt > epsilon
    v       = v0/sqrt(i); % we only need v-seq to be summable
    beta    = beta0*sqrt(i);
    L       = L1 + beta*L2;

    eps_k           = min(epsilon/2,sqrt(0.5*rho/beta));

    [x,ngrad,subiter,inner_err] =  innerNesterov(x,x,A,B,b,y, beta,rho, Q0,c0,down,up,mu,L,ngrad,eps_k,I);

    [alpha, id]   = min([beta,v/(norm(A*x-b)+10^(-15))]);
    y       = y + alpha*(A*x-b); % update multiplier y

    [primal, dua] = getouterviol(x,Q0,c0,A,b,y,up,down);
    kkt           = sqrt(primal^2 + dua^2);
    primalviol(i) = primal;
    dualviol(i)   = dua;
    objective(i)  = 0.5*x'*Q0*x + c0'*x;
    numgrad(i)    = ngrad;
    betaval(i)    = beta;
    alphaval(i)   = alpha;
    
    if kkt < epsilon
    fprintf(['NumIter: %d, objfun: %.5f, dualres: %d, primalres: %d, TotalSubIter: %d, ' ...
        'InnerErr: %d, L : %d, Beta: %d, epsk: %d\n'],i,objective(i),dua,primal,...
        subiter,inner_err,L, beta, eps_k);
    end
    i       = i + 1;

end
%compute the required oututs
out            = [];
out.primalviol = primalviol;
out.dualviol   = dualviol;
out.objective  = objective;
out.gradnum    = numgrad;%number of gradient computations
out.id         = hist_id;
out.betaval    = betaval;
out.alphaval   = alphaval;
end


%inner subproblem solver
function [x,ngrad,i,err] = innerNesterov(x,xk,A,B,b,yk, beta,rho, Q0,c0,down,up,mu,L,ngrad,eps,I)
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
    alpha    = sqrt(mu/L);
    q        = mu/L;
    i        = 0;
    err      = 10;
    step     = 1/L;
    halfgrad = c0 + A'*yk - beta*A'*b - 2*rho*xk;
    C        = Q0+beta*B+ 2*rho*I;
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
            err      = getinnerviol(halfgrad,x,C,up,down);
        end
        ngrad    = ngrad + 1;
    end
end


% finding the error in the first order optimality condition for the inner
% subproblem
function innerviol = getinnerviol(halfgrad,x,C,up,down)
% grad = Q0*x + c0 + A'*y + beta*A'*(A*x-b) + 2*rho*(x-xk); halfgrad = c0 + A'*yk - beta*A'*b - rho*zk;
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

%calculating the KKT violation of the original problem
function [prim, dual] = getouterviol(x,Q0,c0,A,b,y,up,down)

grad = Q0*x + c0 + A'*y ;
xi   = zeros(length(x),1);

id1 = x > down & x < up;
xi(id1) = 0;

id2 = x == down;
xi(id2) = min(0,-grad(id2));

id3 = x == up;
xi(id3) = max(0,-grad(id3));

dual   = norm(xi+grad) ;
prim   = norm(A*x-b);
end


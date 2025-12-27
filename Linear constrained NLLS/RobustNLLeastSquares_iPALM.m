function [out] = RobustNLLeastSquares_iPALM(opts,x,y,beta0,epsilon,m,n,rho,maxiter,cons_nu,u0)
%# Problem setting
% we are solving a quadratic problem with non convex quadratic objective and a
% convex constraint, i.e. min_x \|f(x)\|_1  s.t Ax = b, x \in [down,up], where f is a smooth mapping

%Parameters
% x is the primal variable, y and z are dual variables
% beta0 is the initial parameter for beta = beta0*sqrt(i)
% rho is the proximal parameter
% low and up are the box parameters for x
% L is the Lipschitz constant of f

if isfield(opts,'Q');    Q   = opts.Q;    else,  Q   = randn(n,n);     end
if isfield(opts,'c');    c   = opts.c;    else,  c   = randn(n,1);     end
if isfield(opts,'A');    A    = opts.A;    else,  A    = randn(m,n);   end
if isfield(opts,'b');    b    = opts.b;    else,  b    = randn(m,1);   end
if isfield(opts,'down');  down = opts.down;  else,  down = -10;        end
if isfield(opts,'up');    up   = opts.up;    else,  up   = 10;         end 

v0 = 10;
L  = 10^3;
B  = A'*A;
atb = A'*b;
ngrad = 0;
mu    = 2*rho;
gammau = 1.3;
gammad = 1.5;
subit  = 1000;
err    = 1;
i      = 1;
obj    = zeros(1000,1);
dual   = obj;
primal = obj;
gradnum = obj;
while i < maxiter && err > epsilon 
    xold      = x;
    v         = v0/sqrt(i);
    beta      = beta0 * sqrt(i);
%     eps_k     = min([epsilon/2,epsilon*sqrt(rho), sqrt(rho),epsilon/(sqrt(beta*i))]);
    eps_k     = 1e-3;
    if cons_nu == "True"
%         nu = epsilon^2/rho^3;
        u = u0;
    else
        u         = u0/i;
    end
    
    %update primal variable, x, using inner Nesterov solver
    [x,ngrad,subiter, suberr,L,Ak] =  inner(x, y, A, b, B, atb, beta, Q, c, rho,down,up,mu,u,L,eps_k,subit,gammau,gammad,ngrad,m,n);

    alpha     = min(beta,v/(norm(A*x-b)+10^(-15)));
    y         = y + alpha*(A*x-b); % update multiplier y
    
    vectfun   = zeros(m,1);
    for j = 1:m
        vectfun(j) = 0.5*x'*Q{j}*x + c{j}'*x;
    end
    objective    = norm(vectfun,1);

    obj(i)       = objective;
    dual(i)      = rho*norm(x-xold);
    primal(i)    = norm(A*x-b);
    gradnum(i)   = ngrad;

    err          = max(primal(i), dual(i));

    fprintf(['NumIter: %d, Obj value: %d, ' ...
        'suberr: %d, subiteration: %d, norm(x-xold): %d, Lipschitz: %d, primal err: %d, Err : %d\n'],i,objective,suberr,subiter,dual(i),L, norm(A*x-b), err);

    i         = i + 1;

    

end
out = [];
out.primal = primal;
out.dual   = dual;
out.ngrad  = gradnum;
out.obj    = obj;
end

%inner subproblem solver
function [x,ngrad,i,err,L1, Ak] = inner(x, y, A, b, B, atb, beta , Q, c, rho,down,up,mu,nu,L,eps_k,subit,gammau,gammad,ngrad,m,n)
        
%inner subproblem setting 

% f(x) is the smooth part of the Proximal Augmented Lagrangian function 
% \Psi(x) is the indicator function on the box constraint
%xk is value of x at k-th outer iteration

%initializing parameters
xk      = x;
v       = xk;
gradsum = zeros(length(xk),1);
err     = 10;
i       = 0;
Ak      = 0;
L_max   = 1e5;

% calculating the Jacobian at x_k
D       = zeros(m,n);
for j = 1:m
    D(j,:) = Q{j}*xk + c{j};
end

%calculating one iteration of the required variables for the line search
L       = L*gammau;
ak      = ((1+mu*Ak)+sqrt((1+mu*Ak)^2 + 2*L*(1+mu*Ak)*Ak))/L;
p       = (Ak*x+ak*v)/(Ak + ak);
gradp   = getgrad(A,b,D,B,Q,c,xk,p,y,beta,nu,m,n,rho);
T       = getT(xk,p,gradp,rho,L,down,up);
gradT   = getgrad(A,b,D,B,Q,c,xk,T,y,beta,nu,m,n,rho);
phiprim = L*(p-T)+gradT-gradp;
ngrad   = ngrad + 2;

%inner iteration loop 
while i < subit && err > eps_k
    %while loop for line search
    while (dot(phiprim, p - T) < (1/L)*norm(phiprim,"fro")^2) && L < L_max
        L       = L*gammau;
        ak      = ((1+mu*Ak)+sqrt((1+mu*Ak)^2 + 2*L*(1+mu*Ak)*Ak))/L;
        p       = (Ak*x+ak*v)/(Ak + ak);
        gradp   = getgrad(A,b,D,B,Q,c,xk,p,y,beta,nu,m,n,rho);
        T       = getT(xk,p,gradp,rho,L,down,up);
        gradT   = getgrad(A,b,D,B,Q,c,xk,T,y,beta,nu,m,n,rho);
        phiprim = L*(p-T)+gradT-gradp;
        ngrad   = ngrad + 2;
    end
    
    gradx   = getgrad(A,b,D,B,Q,c,xk,x,y,beta,nu,m,n,rho); %calculate gradient at x_k
    gradsum = gradsum + ak * gradx; %use this to calculate v_k
    Ak      = Ak + ak; % Update A using A_k+1 = A_k + a_k+1
    ngrad   = ngrad+1;
    
    L1 = L;
    
    %update v, L and x after linesearch is satisfied
    v       = getv(Ak,xk,rho,gradsum,down,up);
    L       = L/gammad;
    x       = T;
    if rem(i,10) == 0
        err = getinnerviol(A,b,D,B,Q,c,xk,x,y,beta,nu,m,n,rho,up,down);
    end
    i = i + 1;
end
end

% gradient of the smooth part of the subproblem
function grad = getgrad(A,b,jacobianofc,B,Q,c,xk,x,y,beta,v,m,n,rho)
valueofc    = zeros(m,1);
for i = 1:m
    valueofc(i)   = 0.5*x'*(Q{i}*x +2*c{i});
end
yk   = valueofc + jacobianofc*(x-xk);

proxofl = sign(yk).*max(abs(yk)-v,0);
grad    = (1/v)*(jacobianofc'*(yk-proxofl))+ A'*y + beta*A'*(A*x-b);
end
% 
% function grad = getgrad(A,b,jacobianofc,B,atb,Q,c,xk,x,y,beta,u,m,n,rho)
%     
%     valueofc = zeros(m,1);
%     for i = 1:m
%         valueofc(i) = 0.5*xk'*(Q{i}*xk+2*c{i});
%     end
%     yk    = valueofc + jacobianofc*(x - xk); 
%     proxofl = getinnerprox(yk, u);
%     grad  = A'*y + beta*(B*x-atb)...
%             + (1/u)*jacobianofc'*(yk-proxofl);
% %     grad  = A'*y + beta*(B*x-atb)...
% %             + jacobianofc'*gramoreauforabs(yk,u);
% end

function innerprox = getinnerprox(y,v)
innerprox = sign(y).*max(abs(y)-v,0);
end

function f = gramoreauforabs(x,mu)
    absx = abs(x);
    if absx>mu
        f = sign(x) ; 
    else
        f = x/mu;
    end
end

%finding the proximal mapping of \Psi(x) + 0.5*norm(x-xk)^2
function prox = getprox(xk,down,up,alpha,w,rho)
rho     = 2*rho;
n       = length(xk);
prox    = zeros(n,1);
precalx = (w+alpha*rho*xk)/(1+alpha*rho);

id1 = precalx < up & precalx > down;
prox(id1) = precalx(id1);

id2 = -down + w + alpha*rho*(-down+xk) <= 0;
prox(id2) = down;

id3 = -up + w + alpha*rho*(-up+xk) >= 0;
prox(id3) = up;
end

%calculate the function value at w using the proximal mapping
function T = getT(xk,p,gradp,rho,L,down,up)
alpha = 1/L;
w     = p - alpha*gradp;
T     = getprox(xk,down,up,alpha,w,rho);
end

%calculate v using the proximal mapping and v_k+1 = min_x \psi_k+1(x):=
%\psi_k(x) + a_k+1[f(x_k+1)+< \nabla f(x_{k+1}, x - x_{k+1})> + \Psi(x)]
function v = getv(Ak,xk,rho,gradsum,down,up)
alpha = Ak;
w     = xk - gradsum;
v     = getprox(xk,down,up,alpha,w,rho);
end

% finding the error in the first order optimality condition for the inner
% subproblem
function innerviol = getinnerviol(A,b,D,B,Q,c,xk,x,y,beta,u,m,n,rho,up,down)
% grad = Q0*x + c0 + A'*y + beta*A'*(A*x-b) + 2*rho*(x-xk); 
grad = getgrad(A,b,D,B,Q,c,xk,x,y,beta,u,m,n,rho) + 2*rho*(x - xk);
xi   = zeros(length(x),1);

id1 = x > down & x < up;
xi(id1) = 0;

id2 = x == down;
xi(id2) = min(0,-grad(id2));

id3 = x == up;
xi(id3) = max(0,-grad(id3));

innerviol = norm(grad+xi);
end





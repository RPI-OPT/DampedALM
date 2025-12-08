function [out] = IPLA(Q0, Q, c0, c,d,down,up,sigma,opts )

%# Problem setting
% we are solving a quadratic problem with non convex quadratic objective and a
% quadratic convex constraint, i.e. min_x 0.5*x'Q_0x + c_0'x  s.t
% 0.5*x'Q_jx + c_j'x + d_j <= 0, x \in [down,up]

%Parameters
% z is the primal variable, x is dual variable
% beta0 is the penalty parameter
% lambda is the proximal parameter
% down and up are the box parameters for x
% L is the initial guess for the Lipschitz constant of smooth part of
% subroblem's objective function

if isfield(opts,'x0');        x0       = opts.x0;          else,  x0       = randn(n,1);       end
if isfield(opts,'z0');        z0       = opts.z0;          else,  z0       = randn(m,1);       end
if isfield(opts,'tol');       tol      = opts.tol;         else,  tol      = 1e-3;             end
if isfield(opts,'maxit');     maxit    = opts.maxit;       else,  maxit    = 100;              end
if isfield(opts,'maxsubit');  maxsubit = opts.maxsubit;    else,  maxsubit = 10000;            end 
if isfield(opts,'L');         L        = opts.L;           else,  L        = norm(Q0);         end 
if isfield(opts,'lambda');    lambda   = opts.lambda;      else,  lambda   = -2*min(eig(Q0));  end
if isfield(opts,'gammau');    gammau   = opts.gammau;      else, gammau    = 3;                end
if isfield(opts,'gammad');    gammad   = opts.gammad;      else,  gammad   = 5;                end

%Initialization of variables 
z      = x0;%primal variable
p      = z0;%dual variable 
khat   = 0;
k      = 1;
eps_k   = 10^(-4);
err     = 1;
rho     = 0.5/lambda;%strong convexity constant of the non-smooth part of the subproblem's objective
ngrad   = 0;
mu      = rho;
Lf      = L;
Bg = 0;
for i = 1:length(Q)
    Bg = Bg + norm(Q{i}) + norm(c{i});
end

beta    = max(1, Lf/Bg^2);
sigma = sqrt(0.3);
nu    = sqrt(sigma*(lambda*Lf + 1));

tol1 = 1e-5;
primalviol = zeros(maxit,0);
dualviol   = primalviol;
compslack  = primalviol;
gradnum    = primalviol;
obj        = primalviol;
while k < maxit && err > tol
    
    %solve the subproblem using the Nesterov's gradient method
    [z,ngrad,i,suberr] = inner(z, p, beta , Q0, Q, c0, c,d, rho,down,up,mu,L,eps_k,maxsubit,gammau,gammad,ngrad);
    
    %update the dual variable
    pold     = p;
    g        = zeros(length(p),1);
    for j = 1:length(p)
        g(j) = 0.5*z'*Q{j}*z + c{j}'*z + d(j);
    end
    p        = max(0,p + beta*g);
    
    %calculate the kkt violation and objective value of the main problem
    [primal, dual, compsla] = getouterviol(z,Q0,Q,c0,c,d,p,up,down);
    err                     = sqrt(primal^2+dual^2+compsla^2);
    
    %store some useful values
    primalviol(k) = primal;
    dualviol(k)   = dual;
    compslack(k)  = compsla;
    gradnum(k)    = ngrad;
    obj(k)        = 0.5*z'*Q0*z + c0'*z;

   % Print the inner and outer violations, subiteration number and error
    fprintf(['NumIter: %g, objfun: %g, dualres: %.7f, primalres: %g, ' ...
        'compslack: %0.7f, TotalSubIter: %g, Subproblem err: %g, beta: %g\n'],k,obj(k),dual,primal,compsla,i,suberr, beta);
    
    %store the primal and dual values needed for the khat cycle
    if k == khat + 1
        zkhatplus = z;
        pkhat     = pold;
        
    end
    
    % check if the conditions are satisfied to update beta and go to the
    % next cycle
    if k > khat +1
        deltak = (1/(k-khat-1))*(LagrangianValue(zkhatplus,pkhat, Q0, c0, Q, c, d,beta)- ...
                 LagrangianValue(z,p, Q0, c0, Q, c, d,beta)-0.5*norm(p)^2/beta);
        
        %update beta if updating condition satisfied
        if abs(deltak) <= 0.25*lambda * ( 1 - sigma^2) * tol1^2/(1+2*nu)^2
            beta = 2 * beta;
            khat = k;
        end
    end
    k  = k + 1;
end
out = [];
out.primalviol = primalviol;
out.dualviol   = dualviol;
out.compslack  = compslack;
out.gradnum    = gradnum;
out.objvalue   = obj;
end

% function [z,u,eta,j] = acg(z,pk,Q0,c0, Q, c,d, beta,sigmabar,mutilde,L,down,up,lambda)
% zk      = z;
% x0     = z;
% x      = x0;
% A      = 0;
% tau    = 1;
% j      = 0;
% zeta   = 1/(L-mutilde);
% j      = 0;
% 
% a      = 0.5*(zeta*tau+sqrt((zeta*tau)^2+4*tau*A));
% A      = A + a;
% xtilde = (A*z+a*x)/A;
% tau    = tau + mutilde*a;
% alpha  = 1/L;
% gradz  = getgrad(z, zk, pk,Q0, c0, Q, c, d,beta, lambda);
% w      = z- alpha*gradz;
% z      = getprox(zk,down,up,alpha, w, lambda);
% x      = (1/tau)*((a/zeta)*(z-xtilde)+mutilde*a*z+tau*x);
% 
% u     = mutilde*(z-x)+(x0-x)/A;   
% eta   = (0.5/A)*(norm(x0-z)^2-tau*norm(x-z)^2);
% 
% gammau = 3;
% gammad = 10;
% while norm(u)^2 + 2*eta > sigmabar^2*norm(zk-z+u)^2 && j < 10000
% 
%    zstep     = z - (1/L)*(gradz+(0.5/lambda)*(z-zk));
%    lag_zstep = getsubobj(zstep,zk,pk, Q0, c0, Q, c, d,beta,lambda);
%    lag_z      = getsubobj(z,zk,pk, Q0, c0, Q, c, d,beta,lambda);
%     
%    while lag_zstep > lag_z - (1/L)*norm(gradz)^2
%         L         = L*gammau;
%         zstep     = z - (1/L)*(gradz+(0.5/lambda)*(z-zk));
%         lag_zstep = LagrangianValue(zstep,pk, Q0, c0, Q, c, d,beta);
% 
%    end
%     a      = 0.5*(zeta*tau+sqrt((zeta*tau)^2+4*tau*A));
%     A      = A + a;
%     xtilde = (A*z+a*x)/A;
%     tau    = tau + mutilde*a;
%     gradz  = getgrad(z, zk, pk,Q0, c0, Q, c, d,beta, lambda);
%     w      = z - alpha*gradz;
%     z      = getprox(zk,down,up,alpha, w, lambda);
%     x      = (1/tau)*((a/zeta)*(z-xtilde)+mutilde*a*z+tau*x);
%     
%     u     = mutilde*(z-x)+(x0-x)/A;   
%     eta   = (0.5/A)*(norm(x0-z)^2-tau*norm(x-z)^2);
%     L     = L/ gammad;
%     j     = j + 1;
% end
% 
% end

%inner subproblem solver
% function [x,ngrad,i,err] = innerNesterov(x,xk,Q,c,d,zk, beta,rho, Q0,c0,down,up,mu,L,ngrad,eps,I)
% %inner subproblem setting 
% 
% % f(x) is the smooth part of the objective, f(x) = 0.5*x'*Q0*x + c0'*x +
% % y'(A*x-b) + 0.5*beta*norm(A*x-b)^2 + 0.5*rho*norm(x-xk)^2
% % \Psi(x) = \tau_[-10,10](x) + 0.5*norm(x-xk)^2
% %xk is the previous value of x
% %a_i and A_i are parameters where A_{i+1} = A_{i}+a_{i+1}
% %v_i is argmin \psi_i(x) where psi_i(x) is an estimate function updated as
% %\psi_(i+1)(x) = \psi_i(x) + a_{i+1}(x)[f(x_{i+1})+<\nabla f(x_{i+1}),x-x_{i+1}>+\Psi(x)]
% % T_L(y) is given by prox_{1/L \Psi(x)}(y - (1/L)*\nabla f(y)) 
%     lambda   = 0.5/rho;
%     xhat     = x;
%     xk       = x;
%     alpha    = sqrt(mu/L);
%     q        = mu/L;
%     i        = 0;
%     err      = 1;
%     step     = 1/L;
%     grad     = getgrad(x, xk, zk,Q0, c0, Q, c, d,beta, lambda);
%     gammau   = 3;
%     gammad   = 5;
%     while i < 10000 && err > eps
% 
%         i        = i + 1;
%        gradx     = grad + (0.5/lambda)*(x-xk);
%        xstep     = x - (1/L)*(gradx+(0.5/lambda)*(x-xk));
%        lag_xstep = getsubobj(xstep,xk,zk, Q0, c0, Q, c, d,beta,lambda);
%        lag_x     = getsubobj(x,xk,zk, Q0, c0, Q, c, d,beta,lambda);
%         
%        while lag_xstep > lag_x - (1/L)*norm(gradx)^2
%             L         = L*gammau;
%             xstep     = x - (1/L)*(gradx+(0.5/lambda)*(x-xk));
%             lag_xstep = getsubobj(xstep,xk,zk, Q0, c0, Q, c, d,beta,lambda);
%   
%        end
% 
%         grad     = getgrad(x, xk, zk,Q0, c0, Q, c, d,beta, lambda);
%         w        = xhat - step*grad;
%         xold     = x;
%         x        = getprox(xk, down, up, step,w, rho);
%         alphaold = alpha;
%         alpha    = 0.5*(q-alpha^2 +sqrt((q-alpha^2)^2 + 4*alpha^2));
%         xhat     = x + ((alphaold*(1-alphaold))/(alphaold^2+alpha))*(x-xold);
% 
%         if rem(i,20) == 0
%             err      = getinnerviol(xk,x,Q0,Q,c0,c,d,beta,zk,rho,up,down);
%         end
%         ngrad    = ngrad + 1;
%         L        = L/gammad;
%     end
% end
% 
% 
% function value = getsubobj(x,xk,z, Q0, c0, Q, c, d,beta,lambda)
% m = length(z);
% sum = 0;
% for i = 1:m
%     
%     sum = sum + max(0,x'*Q{i}*x+c{i}'*x+ d(i))^2;
% end
% value = 0.5*x'*Q0*x + c0'*x +sum - (0.5/beta)*norm(z)^2 + (0.5/lambda)*norm(x-xk)^2;
% end
% 
% function grad = getgrad(x, xk, zk,Q0, c0, Q, c, d,beta, lambda)
% sum = 0;
% m = length(zk);
% for i  = 1:m
%     term = Q{i}'*x;
%     sum = sum + max(0,zk(i)+beta*(x'*(0.5*term + c{i}) + d(i)))*(term +c{i});
% end
% grad = Q0*x + c0  + sum + (0.5/lambda)*(x-xk);
% end
% 
% 
% %finding the proximal mapping of \Psi(x) = \tau_{[-10,10]}(x) +
% %0.25*gamma*(x-z^k)^2 wrt parameter alpha at w
% function prox = getprox(xk,down,up,alpha,w,lambda)
% n = length(xk);
% prox = zeros(n,1);
% precalx = (2*lambda*w+alpha*xk)/(2*lambda+alpha);
% 
% id1 = precalx < up & precalx > down;
% prox(id1) = precalx(id1);
% 
% id2 = 2*lambda*(-down + w) + alpha*(-down+xk) <= 0;
% prox(id2) = down;
% 
% id3 = 2*lambda*(-up + w) + alpha*(-up+xk) >= 0;
% prox(id3) = up;
% end
% 
% 
% function [primal, dual, compli] = getouterviol(x,Q0,Q,c0,c,d,z,up,down)
% 
% sum1 = 0;
% sum3 = sum1;
% sum4 = zeros(length(z),1);
% for j = 1:length(z)
%     sum1    = sum1 + z(j)*(Q{j}*x + c{j});
%     sum3    = sum3 + abs(z(j)*(0.5*x'*Q{j}*x+c{j}'*x + d(j)));
%     sum4(j) = max(0,0.5*x'*Q{j}*x+c{j}'*x + d(j));
% end
% grad = Q0*x + c0 + sum1 ;
% xi   = zeros(length(x),1);
% 
% id1 = x > down & x < up;
% xi(id1) = 0;
% 
% id2 = x == down;
% xi(id2) = min(0,-grad(id2));
% 
% id3 = x == up;
% xi(id3) = max(0,-grad(id3));
% 
% primal  = norm(grad+xi);
% dual    = norm(sum4);
% compli  = sum3;
% 
% end
% 
% function innerviol = getinnerviol(xk,x,Q0,Q,c0,c,d,beta,zk,rho,up,low)
% rho = 2*rho;
% lambda = 0.5/rho;
% grad = getgrad(x, xk, zk,Q0, c0, Q, c, d,beta, lambda);
% 
% xi = zeros(length(x),1);
% id1 = x > low & x < up;
% xi(id1) = 0;
% 
% id2 = x == low;
% xi(id2) = min(0,-grad(id2));
% 
% id3 = x == up;
% xi(id3) = max(0,-grad(id3));
% 
% innerviol = norm(grad+xi);
% end
% 
%inner subproblem solver
function [x,ngrad,i,err] = inner(x, z, beta , Q0, Q, c0, c,d, rho,down,up,mu,L,eps_k,subit,gammau,gammad,ngrad)
        
%inner subproblem setting 

% f(x) is the smooth part of the objective, f(x) = 0.5*x'*Q0*x + c0'*x +
% y'(A*x-b) + 0.5*beta*norm(A*x-b)^2 + 0.5*rho*norm(x-xk)^2
% \Psi(x) = \tau_[-10,10](x) + 0.5*rho*norm(x-xk)^2
%xk is the previous value of x
%a_i and A_i are parameters where A_{i+1} = A_{i}+a_{i+1}
%v_i is argmin \psi_i(x) where psi_i(x) is an estimate function updated as
%\psi_(i+1)(x) = \psi_i(x) + a_{i+1}(x)[f(x_{i+1})+<\nabla f(x_{i+1}),x-x_{i+1}>+\Psi(x)]
% T_L(y) is given by prox_{1/L \Psi(x)}(y - (1/L)*\nabla f(y)) 

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
%     if L < 5
%         keyboard;
%     end
    if rem(i,5) == 0
     err     = getinnerviol(xk,x,Q0,Q,c0,c,d,beta,z,rho,up,down);
    end

    i = i + 1;
end
% disp(err)
% disp(i)
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

%calculate the first order optimality violation in the subproblem 
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

%Calculate the Lagrangian value at given x,z and beta
function value = LagrangianValue(x,z, Q0, c0, Q, c, d,beta)
m = length(z);
sum = zeros(m,1);
for i = 1:m   
    sum(i) =  max(0,x'*Q{i}*x+c{i}'*x+ d(i)+ z(i)/beta);
end
value = 0.5*x'*Q0*x + c0'*x +(0.5/beta)* norm(sum)^2 - (0.5/beta)*norm(z)^2;
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


%calculate the kkt violation of the original problem
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


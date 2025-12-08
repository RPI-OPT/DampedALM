function out = DPALM(Au, Ap, A,b, maxiter, epsilon, L, x,...
               lambda, beta0, rho, gammau, gammad, theta, mu, u0, const_nu, kappa, Lstar, down, up )

% initializing variables, Ap and Au represent data for "protected" and
% "unprotected" group, 
[np,~]  = size(Ap);
[nu,~]  = size(Au);
n       = np + nu;
i       = 1;
err     = 1;
v0      = 1000;
obj     = zeros(maxiter,1);
dual    = obj;
primal  = obj;
gradnum = obj;
multipl = obj;
subit   = 1000;
ngrad   = 0;
while i < maxiter && err > epsilon 
    xk      = x;
    lamold    = lambda;
    v         = v0/sqrt(i);
    beta      = beta0 * sqrt(i);
    eps_k     = min(epsilon,1/(beta*i));
    if const_nu == "True"
        u = u0;
    else
        u = u0/i;
    end
    [x,ngrad,subiter, suberr,L,Ak] =  inner(x, lambda, Au, Ap, A, b, theta, beta, np, nu,...
                                      rho, mu, u, L, eps_k,subit,gammau,gammad,ngrad, kappa, Lstar, down, up);

    alpha      = min(beta, v/(max(0, (0.5/n)*norm(A*x-b)^2 - Lstar - kappa)));
    lambda     = lambda + beta* max(-lambda/beta, (0.5/n)*norm(A*x-b)^2 - Lstar - kappa); % update multiplier lambda
    
    objective  = (0.5/n)*norm(A*x-b)^2;
    obj(i)     = objective;
    dual(i)    = rho*norm(x-xk);
%     primal(i)  =  lambda-lamold;
    primal(i)  =  max(0, (0.5/n)*norm(A*x-b)^2 - Lstar - kappa);
    gradnum(i) = ngrad;
    multipl(i) = lambda;

    err        = sqrt(primal(i)^2 + dual(i)^2);

    fprintf(['NumIter: %d, Obj value: %d, ' ...
                 'suberr: %d, subiteration: %d, norm(x-xold): %d, Lipschitz: %d, primal err: %d, Err : %d\n'],...
                 i,objective,suberr,subiter,dual(i),L, primal(i), err);

    i         = i + 1;

end
out = [];
out.x      = x;
out.primal = primal;
out.dual   = dual;
out.ngrad  = gradnum;
out.obj    = obj;
out.mult   = multipl;
end


function [x,ngrad,i,err,L1, Ak] = inner(x, lambda, Au, Ap, A, b, theta, beta, np, nu,...
                                    rho,mu, u, L,eps_k,subit,gammau,gammad,ngrad, kappa, Lstar, down, up)
        
%inner subproblem setting 

% f(x) is the smooth part of the objective, f(x) = 0.5*x'*Q0*x + c0'*x +
% y'(A*x-b) + 0.5*beta*norm(A*x-b)^2 + 0.5*rho*norm(x-xk)^2
% \Psi(x) = \tau_[-10,10](x) + 0.5*rho*norm(x-xk)^2
%oldx is the previous value of x
%a_i and A_i are parameters where A_{i+1} = A_{i}+a_{i+1}
%v_i is argmin \psi_i(x) where psi_i(x) is an estimate function updated as
%\psi_(i+1)(x) = \psi_i(x) + a_{i+1}(x)[f(x_{i+1})+<\nabla f(x_{i+1}),x-x_{i+1}>+\Psi(x)]
% T_L(y) is given by prox_{1/L \Psi(x)}(y - (1/L)*\nabla f(y)) 

%initializing parameters
xk      = x;
v       = xk;
gradsum = zeros(length(xk),1);
err     = 10;
i       = 0;
Ak      = 0;

% calculating the Jacobian at x_k to calculate the gradient of the smooth
% part
[jacobofc, valueofc] = getvalandgrad(Au, Ap, theta, xk, np, nu);

%calculating one iteration of the required variables for the line search
L       = L*gammau;
ak      = ((1+mu*Ak)+sqrt((1+mu*Ak)^2 + 2*L*(1+mu*Ak)*Ak))/L;
p       = (Ak*x+ak*v)/(Ak + ak);
gradp   = getgrad(jacobofc, valueofc, A, b, p, xk, u, rho, beta, lambda, kappa, Lstar);
T       = getT(xk,p,gradp, rho, L, down, up);
gradT   = getgrad(jacobofc, valueofc, A, b, T, xk, u, rho, beta, lambda, kappa, Lstar);
phiprim = L*(p-T)+gradT-gradp;
ngrad   = ngrad + 2;

%inner iteration loop 
while i < subit && err > eps_k
    %while loop for line search
    while (dot(phiprim, p - T) < (1/L)*norm(phiprim,"fro")^2)
        L       = L*gammau;
        ak      = ((1+mu*Ak)+sqrt((1+mu*Ak)^2 + 2*L*(1+mu*Ak)*Ak))/L;
        p       = (Ak*x+ak*v)/(Ak + ak);
        gradp   = getgrad(jacobofc, valueofc, A, b, p, xk, u, rho, beta, lambda, kappa, Lstar);
        T       = getT(xk,p,gradp, rho, L, down, up);
        gradT   = getgrad(jacobofc, valueofc, A, b, T, xk, u, rho, beta, lambda, kappa, Lstar);
        phiprim = L*(p-T)+gradT-gradp;
        ngrad   = ngrad + 2;
    end
    
    gradx   = getgrad(jacobofc, valueofc, A, b, x, xk, u, rho, beta, lambda, kappa, Lstar); %calculate gradient at x_k
    gradsum = gradsum + ak * gradx; %use this to calculate v_k
    Ak      = Ak + ak; % Update A using A_k+1 = A_k + a_k+1
    ngrad   = ngrad+1;
   
    %update v, L and x after linesearch is satisfied
    v       = getv(Ak,xk,rho,gradsum,down,up);
    L       = L/gammad;
    x       = T;
    if rem(i,1) == 0
        err = getinnerviol(jacobofc, valueofc, A, b, x, xk, u, rho, beta, lambda, kappa, Lstar, down, up);
%         fprintf('L: %f, Error: %f \n ', L1, err)
    end
    i = i + 1;
end
L1 = L;
end

function grad = getgrad(jacobofc, valueofc, A, b, x, xk, u, rho, beta, lambda, kappa, Lstar)
    [n,~] = size(A);
    yk  = valueofc + jacobofc'*(x - xk);

    Axminusb = A*x-b;
    grad     = (1/u)*(jacobofc * (yk - sign(yk)*max(abs(yk)-u, 0)))...
               + (beta/n) * A'*(Axminusb) * (max((0.5/n)*norm(Axminusb)^2 - Lstar - kappa + lambda/beta, 0))...
               + rho*(x - xk);
end

function [jacobofc, valueofc] = getvalandgrad(Au, Ap, theta, xk, np, nu)
    grad1 = zeros(size(xk));
    grad2 = grad1;
    % defining sigmoid and gradient of sigmoid
    sigmoid      = @(y) exp(y)/ (1 + exp(y));
    gradsigmoid  = @(y) sigmoid(y)*(1-sigmoid(y));
    
    % Calculating the gradient of inner function, c(x) at oldx
    for i = 1:np
        grad1 = grad1 + gradsigmoid(Ap(i,:)* xk - theta)*Ap(i, :)';
    end
    for i = 1:nu
        grad2 = grad2 + gradsigmoid(Au(i,:)* xk - theta)*Au(i, :)';
    end
    jacobofc = (1/np)*grad1 - (1/nu)*grad2;

    % Calculating the value of inner function, c(x) at oldx
    cxk1 = 0;
    cxk2 = 0;
    for i = 1:np
        cxk1 = cxk1 + sigmoid(Ap(i,:)*xk- theta);
    end
    for i = 1:nu
        cxk2 = cxk2 + sigmoid(Au(i,:)*xk- theta);
    end
    valueofc = (1/np)*cxk1 - (1/nu)*cxk2;
end

%finding the proximal mapping of \Psi(x) + 0.5*norm(x-xk)^2
function prox = getprox(xk,down,up,alpha,w,rho)
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
function innerviol = getinnerviol(jacobofc, valueofc, A, b, x, xk, u, rho, beta, lambda, kappa, Lstar, down, up)

grad = getgrad(jacobofc, valueofc, A, b, x, xk, u, rho, beta, lambda, kappa, Lstar) + rho*(x - xk);
xi   = zeros(length(x),1);

id1 = x > down & x < up;
xi(id1) = 0;

id2 = x == down;
xi(id2) = min(0,-grad(id2));

id3 = x == up;
xi(id3) = max(0,-grad(id3));

innerviol = norm(grad+xi);
end










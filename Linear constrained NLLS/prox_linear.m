function out = prox_linear(Au, Ap, A,b, maxiter, epsilon, L, x,...
                rho, gammau, gammad, theta, mu, u0, u1, const_nu, kappa, Lstar, down, up )

% initializing variables, Ap and Au represent data for "protected" and
% "unprotected" group, 
[np,~]  = size(Ap);
[nu,~]  = size(Au);
n       = np + nu;
i       = 1;
err     = 1;
obj     = zeros(maxiter,1);
dual    = obj;
primal  = obj;
gradnum = obj;
subit   = 20000;
ngrad   = 0;
while i < maxiter && err > epsilon 
    xk      = x;
    eps_k     = min(epsilon, 1/i^(1));
    if const_nu == "True"
        u = u0;
    else
        u = u0/i;
    end
%     u1 = u;
    [x,ngrad,subiter, suberr,L,Ak] =  inner(x, Au, Ap, A, b, theta, np, nu,...
                                      rho, mu, u, u1, L,eps_k,subit,gammau,gammad,ngrad, kappa, Lstar, down, up);

    
    dual(i)    = rho*norm(x-xk);
    primal(i)  = max(0, (0.5/n)*norm(A*x-b)^2 - Lstar - kappa);
    gradnum(i) = ngrad;
   

    err        = max(primal(i), dual(i));

    fprintf(['NumIter: %d,' ...
                 'suberr: %d, subiteration: %d, norm(x-xold): %d, Lipschitz: %d, primal err: %d, Err : %d\n'],...
                 i,suberr,subiter,dual(i),L, primal(i), err);

    i         = i + 1;
end
out        = [];
out.x      = x;
out.primal = primal;
out.dual   = dual;
out.ngrad  = gradnum;
out.obj    = obj;
end


function [x,ngrad,i,err,L1, Ak] = inner(x, Au, Ap, A, b, theta, np, nu,...
                                    rho,mu, u, u1, L,eps_k,subit,gammau,gammad,ngrad, kappa, Lstar, down, up)
        
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
err     = 10;
i       = 0;
Ak      = 0;

% calculating the Jacobian and value of c at x_k to calculate the gradient of the smooth
% part
[jacobofc, valueofc] = getvalandgrad(Au, Ap, A, b, theta, xk, np, nu, kappa, Lstar);

%calculating one iteration of the required variables for the line search
ak      = ((1+mu*Ak)+sqrt((1+mu*Ak)^2 + 2*L*(1+mu*Ak)*Ak))/L;
p       = (Ak*x+ak*v)/(Ak + ak);
gradp   = getgrad(jacobofc, valueofc, p, xk, u, u1, rho);
T       = getT(xk,p,gradp, rho, L, down, up);
gradT   = getgrad(jacobofc, valueofc, T, xk, u, u1, rho);
phiprim = L*(p-T)+gradT-gradp;
ngrad   = ngrad + 2;

%inner iteration loop 
while i < subit && err > eps_k
    %while loop for line search
    while (dot(phiprim, p - T) < (1/L)*norm(phiprim,"fro")^2)
        L       = L*gammau;
        ak      = ((1+mu*Ak)+sqrt((1+mu*Ak)^2 + 2*L*(1+mu*Ak)*Ak))/L;
        p       = (Ak*x+ak*v)/(Ak + ak);
        gradp   = getgrad(jacobofc, valueofc, p, xk, u, u1, rho);
        T       = getT(xk,p,gradp, rho, L, down, up);
        gradT   = getgrad(jacobofc, valueofc, T, xk, u, u1, rho);
        phiprim = L*(p-T)+gradT-gradp;
        ngrad   = ngrad + 2;
    end
    
    gradx   = getgrad(jacobofc, valueofc, x, xk, u, u1, rho); %calculate gradient at x_k
    gradsum = gradsum + ak * gradx; %use this to calculate v_k
    Ak      = Ak + ak; % Update A using A_k+1 = A_k + a_k+1
    ngrad   = ngrad+1;
    
    %update v, L and x after linesearch is satisfied
    v       = getv(Ak,xk,rho,gradsum,down,up);
    L       = L/gammad;
    x       = T;

    if rem(i,1) == 0
        err = getinnerviol(jacobofc, valueofc, x, xk, u, u1, rho, down, up);
    end
    i = i + 1;
end
L1 = L;
end

function grad = getgrad(jacobofc, valueofc, x, xk, u, u1, rho)    
    yk          = valueofc + jacobofc*(x - xk);
    proxofl     = getinnerprox(yk,u);
%     grad        = (1/u)*jacobofc'*(yk-proxofl) ...
%                   + rho*(x - xk) ;
    jacobofc(1) = jacobofc(1)/u;
    jacobofc(2) = jacobofc(2)/u1;
    grad        = jacobofc'*(yk-proxofl) ...
                  + rho*(x - xk) ;
end

function [jacobofc, valueofc] = getvalandgrad(Au, Ap, A, b, theta, xk, np, nu, kappa, Lstar)
    n        = np + nu;
    valueofc = zeros(2,1);
    jacobofc = zeros(2,length(xk));
    grad1    = zeros(size(xk));
    grad2    = grad1;
    
    % defining sigmoid and gradient of sigmoid
    sigmoid      = @(y) exp(y)/ (1 + exp(y));
    gradsigmoid  = @(y) sigmoid(y)*(1-sigmoid(y));

    % Calculating the gradient of inner function, c(x) at xk
    for i = 1:np
        grad1 = grad1 + gradsigmoid(Ap(i,:)* xk - theta)*Ap(i, :)';
    end
    for i = 1:nu
        grad2 = grad2 + gradsigmoid(Au(i,:)* xk - theta)*Au(i, :)';
    end
    gradcxk = (1/np)*grad1 - (1/nu)*grad2;
    jacobofc(1,:) = gradcxk;
    jacobofc(2,:) = (1/n)*A'*(A*xk-b);

    % Calculating the value of inner function, c(x) at xk
    cxk1 = 0;
    cxk2 = 0;
    for i = 1:np
        cxk1 = cxk1 + sigmoid(Ap(i,:)*xk- theta);
    end
    for i = 1:nu
        cxk2 = cxk2 + sigmoid(Au(i,:)*xk- theta);
    end
    valueofc(1) = (1/np)*cxk1 - (1/nu)*cxk2;
    valueofc(2) = (0.5/n)*norm(A*xk - b)^2 - Lstar - kappa;

end

function innerprox = getinnerprox(y,u)
    innerprox    = zeros(size(y));
    innerprox(1) = sign(y(1)).*max(abs(y(1)) - u, 0);
    innerprox(2) = min(0, y(2));
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
function innerviol = getinnerviol(jacobofc, valueofc, x, xk, u, u1, rho, down, up)

grad = getgrad(jacobofc, valueofc, x, xk, u, u1, rho) + rho*(x - xk);
xi   = zeros(length(x),1);

id1 = x > down & x < up;
xi(id1) = 0;

id2 = x == down;
xi(id2) = min(0, -grad(id2));

id3 = x == up;
xi(id3) = max(0, -grad(id3));

innerviol = norm(grad+xi);
end










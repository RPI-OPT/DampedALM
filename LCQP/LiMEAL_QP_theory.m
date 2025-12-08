% LiMEAL for quadratic programming
% Considering the following quadratic programming:
% min_x 1/2 x^TQx+r^Tx  s.t. Ax=b, ell_i<=xi<=ui for all i
% Ref: Jinshan Zeng, Wotao Yin and Ding-Xuan Zhou, Moreau envelope
% augmented Lagrangian method for noncovnex optimization with linear
% constraints, JSC, 2021.
function [x,z,lamb,out]=LiMEAL_QP_theory(QP,x0,z0,lamb0,beta,gamma,eta,NumIter,epsilon,I)
% Input:
% QP: the settings of the problem of quadratic programming
% QP.Q: matrix Q; QP.r: vector r; QP.A: matrix A; QP.b: vector b;
% QP.l: lower bound ell; QP.u: upper bound u
% lamb0 -- initialization of Lagrangian multipliers
% x0,z0 -- initialization
% gamma -- proximal parameter (0,1/||Q||) in this case
% eta -- z-step size (0,2), eta=1 then LiMEAL reduces to proximal ALM
% beta -- penalty parameter (sufficiently large), used as dual stepsize
% NumIter -- number of iterations

% Output:
% x,z,lamb -- the outputs of algorithmic iterative sequences
% objfun -- the trend of objective function
% dualres -- dual residual defined as: dualres = ||Ax-b||
% primalres -- primal residual defined as: primalres = ||x(k+1)-z(k)||
% staterr -- stationary error defined as the norm of gradient of Moreau Envelope
objfun = zeros(NumIter,1);
dualres = zeros(NumIter,1);
primalres = zeros(NumIter,1);
staterr = zeros(NumIter,1);

% calculating an inverse of matrix beta*A'A+1/gamma*I
% invAA = (beta*(QP.A)'*QP.A+eye(size(QP.A,2))/gamma)\eye(size(QP.A,2));
invAA = inv(beta*(QP.A)'*QP.A+eye(size(QP.A,2))/gamma);
i     = 1;
% calculating the vector beta*A'b-r
Abr = beta*(QP.A)'*(QP.b)-QP.r;
err = 1;
dualviol1 = [];
% run iterations

L          = norm(beta*(QP.A'*QP.A) + eye(size(QP.A,2)) /gamma);
mu         = 1/gamma;
ngrad = 0;    
B     = QP.A'*QP.A;
%initializing zero vectors to store the values later
primalviol   = zeros(NumIter-1,1);
dualviol     = primalviol;
numgrad      = primalviol;
objective    =  primalviol;
while i < NumIter && err > epsilon 
%     tempx = invAA*(z0/gamma-QP.Q*x0-(QP.A)'*lamb0+Abr);
%     x = (tempx>=QP.u).*QP.u+(tempx<=QP.ell).*QP.ell+(QP.ell<tempx&tempx<QP.u).*tempx;% projection onto [ell,u]
    eps_k   = 1/(beta*i);

%     eps_k = max(1e-2, eps_k);

    eps_k = min(epsilon^(-4), eps_k);  % added by YX, 
                                    % solve subproblem to sufficient accuracy sometimes is very important
    [x, inner_err,ngrad,subiter] = innerNesterov(x0,z0,QP.A,B,QP.b,lamb0, beta, QP.Q,QP.r,gamma,QP.ell,QP.u,mu,L,ngrad,eps_k,I);
    
    z    = z0-eta*(z0-x);
    lamb = lamb0+beta*(QP.A*x-QP.b);
    
%     objfun(i)      = 1/2*x'*QP.Q*x+(QP.r)'*x; % objective function
%     dualres(i)     = norm(QP.A*x-QP.b); % dual residual
%     primalres(i)   = norm(x-z0); % primal residual
    err1           = sqrt(norm(gamma^(-1)*(z0-x)+QP.Q*(x-x0))^2+beta^(-2)*norm(lamb-lamb0)^2);
    staterr(i)     = err1;
    [primal, dual] = getouterviol(x,QP.Q,QP.r,QP.A,QP.b,lamb,QP.u,QP.ell);
%     temp_dualviol = [primal, dual, norm(x - x0), norm(x - z0)/gamma];
%     dualviol1     = [dualviol1; temp_dualviol];

    primalviol(i) = primal;
    dualviol(i)   = dual;
    objective(i)  = 0.5*x'*QP.Q*x + QP.r'*x;
    numgrad(i)    = ngrad;
    % print the main evaluation metrics
    err = sqrt(dual^2+primal^2);
    if err <= epsilon
        fprintf(['NumIter: %g, objfun: %g, dualres: %g, primalres: %g, staterr: %g, ' ...
         'subproblem err: %g, subiter: %d\n'], ...
         i, objective(i), dual, primal, staterr(i), inner_err, subiter);
    end
    i = i+1;

    x0 = x;
    z0 = z;
    lamb0 = lamb;

    
end
%compute the required oututs
out            = [];
out.primalviol = primalviol;
out.dualviol   = dualviol;
out.objective  = objective;
out.gradnum    = numgrad;%number of gradient computations
clear invAA; clear Abr; clear tempx;
end

function [prim,dual] = getouterviol(x,Q0,c0,A,b,y,up,down)

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

%outerviol = sqrt(prim^2+dual^2);
outerviol = [prim, dual];
end

function [x, err,ngrad,i] = innerNesterov(x,zk,A,B,b,yk, beta, Q0,c0,gamma,down,up,mu,L,ngrad,eps_k,I)
    xhat  = x;
    xk    = x;
    alpha = sqrt(mu/L);
    q     = mu/L;
    i     = 0;
    err   = 1;
    step     = 1/L;
    halfgrad1 = Q0*xk + c0 + A'*(yk - beta*b) - zk/gamma;
    halfgrad2 = beta * B + (I./gamma);
    while i < 10000 && err > eps_k
        i        = i + 1;
        grad     = halfgrad1 + halfgrad2*xhat;
        w        = xhat - step*grad;
        xold     = x;
        x        = getprox(down, up, step, w,gamma, zk);
        alphaold = alpha;
        alpha    = 0.5*(q-alpha^2 +sqrt((q-alpha^2)^2 + 4*alpha^2));
        xhat     = x + ((alphaold*(1-alphaold))/(alphaold^2+alpha))*(x-xold);
        if rem(i,10) == 0
            err      = getinnerviol(Q0,c0,x,xk,A,yk,beta,b,zk,gamma,up,down);
        end
%         disp(err)
%         disp(i);
        ngrad    = ngrad + 1;
    end
    
%     if i >= 1000
%         keyboard;
%     end
end


function innerviol = getinnerviol(Q0,c0,x,xk,A,yk,beta,b,zk,gamma,up,down)

grad = Q0*xk +c0 + A'*yk + beta* A'*(A*x-b) + (1/gamma)*(x-zk);

xi = zeros(length(x),1);
id1 = x > down & x < up;
xi(id1) = 0;

id2 = x == down;
xi(id2) = min(0,-grad(id2));

id3 = x == up;
xi(id3) = max(0,-grad(id3));

innerviol = norm(grad+xi);
end

function prox = getprox(down,up,alpha,w,gamma,zk)
% n = length(zk);
% prox = zeros(n,1);
% precalx = (2*gamma*w+alpha*zk)/(2*gamma+alpha);
% 
% id1 = precalx < up & precalx > down;
% prox(id1) = precalx(id1);
% 
% id2 = 2*gamma*(-down + w) + alpha*(-down+zk) <= 0;
% prox(id2) = down;
% 
% id3 = 2*gamma*(-up + w) + alpha*(-up+zk) >= 0;
% prox(id3) = up;
prox = min(up,max(down,w));
end




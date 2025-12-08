function [x,y,out] = HiAPeM_qp(Q, c, A, b, x_l, x_u, opts)

% ====================================================================
%
% Hybrid iALM and Penalty method for solving box-constrained QP
% min_x 0.5*(x'*Q*x)+c'*x s.t. Ax == b; x_l <= x <= x_u.
% Input:
%       model data: Q,c,A,b
%       lower and upper bound: x_l, x_u
%       algorithm parameters:
%           in_iter,
%           tol,
%           betamax
% Output:
%       x and y: primal-dual solution
%

% =====================================================================

At    = A';
AtA   = At*A;
up    = x_u(1);
down  = x_l(1);

nobj  = 0; % numObj
ngrad = 0; % numGrad

[m, n]=size(A);

if isfield(opts,'K0')        K0 = opts.K0;             else  K0 = 100;         end
if isfield(opts,'K1')        K1 = opts.K1;             else  K1 = 2;           end
if isfield(opts,'sig')       sig = opts.sig;           else  sig = 3;          end
if isfield(opts,'gamma')     gamma = opts.gamma;       else  gamma = 1.1;      end
if isfield(opts,'beta0')     beta0 = opts.beta0;       else  beta0 = 1e-2;     end
if isfield(opts,'x0')        x0 = opts.x0;             else  x0 = randn(n,1);  end
if isfield(opts,'tol')       tol = opts.tol;           else  tol = 1e-3;       end
if isfield(opts,'maxit')     maxit = opts.maxit;       else  maxit = 100;      end
if isfield(opts,'APGmaxit')  APGmaxit = opts.APGmaxit; else  APGmaxit = 10000; end
if isfield(opts,'adp')       adp = opts.adp;           else  adp = 1;          end
if isfield(opts,'frc')       frc = opts.frc;           else  frc = 0.5;        end

x = x0;
y = zeros(m,1);

beta = beta0;

rho = min(min(eig(Q)),0); % rho-weakly convex

rho = max(abs(rho),tol); % arbitrary...

s = 1;

K_total = K0;

hist_numPG = [];

%initializing zero vectors to store the values later
primalviol   = zeros(maxit-1,1);
dualviol     = primalviol;
numgrad      = primalviol;
objective    =  primalviol;
% start the initial stage by calling iALM
for k = 1:K0
    xk = x;
    
    [x,y,beta,total_PG_k] = iALM_str(x, xk, tol*frc, rho, beta0);
    
    hist_numPG = [hist_numPG; total_PG_k];
    
    gradL0 = Q*x + c + At*y;
    ngrad = ngrad + 1;

    [prim, dual]  = getouterviol(x,Q,c,A,b,y,up,down);
    primalviol(k) = prim;
    dualviol(k)   = dual;
    objective(k)  = 0.5*x'*Q*x + c'*x;
    numgrad(k)    = ngrad;
    I1 = (x == x_l);
    I2 = (x == x_u);
    I3 = (x > x_l & x < x_u);
    
    dres = norm( min(0,gradL0(I1)) )^2 + norm( max(0,gradL0(I2)) )^2 + norm(gradL0(I3))^2;
    dres = sqrt(dres);
    
    if dres <= tol % only check the dual_res because primal_res is below tol/2
        fprintf('succeed during the initial stage!\n');
        break;
    end
end

if dres > tol
    K_total = K0 + K1;
    s = s + 1;
    y_fix = y; beta_fix = beta;
    
    tol_pen = tol*frc; 
    
end


while dres > tol % only check dual_res since primal_res must below tol
    
    k = k + 1; beta_fix = beta;
    
    if k < K_total
        % call the penalty method
        xk = x;
        [x, y, beta,total_PG_k] = PenMM(xk, y_fix, rho, tol_pen, beta_fix);
        
        hist_numPG = [hist_numPG; total_PG_k];
        
        if norm(x - xk) <= (1-frc)*tol/2/rho
            fprintf('final subproblem solved by PenMM\n');
            break;
        end
        
    else
        % call the iALM method
        xk = x;
        
        [x,y,beta,total_PG_k] = iALM_str(x, xk, tol*frc, rho, beta0);
        
        hist_numPG = [hist_numPG; total_PG_k];
        
        if norm(x - xk) <= (1-frc)*tol/2/rho
            fprintf('final subproblem solved by iALM\n');
            break;
        end
        
        y_fix = y;
        
    end
    
    gradL0 = Q*x + c + At*y;
    
    ngrad = ngrad + 1;
    
    I1 = (x == x_l);
    I2 = (x == x_u);
    I3 = (x > x_l & x < x_u);
    
    dres = norm( min(0,gradL0(I1)) )^2 + norm( max(0,gradL0(I2)) )^2 + norm(gradL0(I3))^2;
    dres = sqrt(dres);
    
    if k == K_total & dres > tol
        Ks = ceil(gamma^s*K1);
        s = s + 1;
        K_total = K_total + Ks;
    end
    
end

if k == K_total | s == 1
    fprintf('final subproblem solved by iALM\n');
else
    fprintf('final subproblem solved by PenMM\n');
end

if norm(x - xk) <= (1-frc)*tol/2/rho
    % compute dres in this case
    % otherwise dres already computed
    gradL0 = Q*x + c + At*y;
    
    ngrad = ngrad + 1;
    
    I1 = (x == x_l);
    I2 = (x == x_u);
    I3 = (x > x_l & x < x_u);
    
    dres = norm( min(0,gradL0(I1)) )^2 + norm( max(0,gradL0(I2)) )^2 + norm(gradL0(I3))^2;
    dres = sqrt(dres);
end

% out.nobj = nobj; % = 0 for non-adaptive!
% out.ngrad = ngrad;
% out.obj = 1/2*x'*Q*x+c'*x;
% out.dres = dres;
% out.pres = norm(A*x-b);
% out.numStage = s;
% out.totalProb = k;
% out.numPG = hist_numPG;

%compute the required oututs
out            = [];
out.primalviol = primalviol;
out.dualviol   = dualviol;
out.objective  = objective;
out.gradnum    = numgrad;%number of gradient computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [x,y,beta,total_PG] = iALM_str(x0, xk, tol, rho, beta0)
        
        % =========================================
        % inexact ALM for solving k-subprob:
        %
        %       min f0(x) + rho|x-xk|^2, s.t. Ax == b.
        %
        % where f0(x) = 0.5x'Qx+c'x
        %
        % AL: f0(x) +  rho|x-xk|^2 + <y,Ax-b> + beta/2*|Ax-b|^2
        %
        % the output (x,y) will satisfy tol - KKT conditions
        %
        % ============================================
        x = x0;        y = zeros(m,1);    beta = beta0;
        
        tolk = 0.999 * min( tol, sqrt((sig-1)/sig*tol*rho)/2 );
        
        if adp
            [x, numPG] = APG_adp(x0, xk, y, rho, beta, tolk);  
        else
            [x, numPG] = APG(x0, xk, y, rho, beta, tolk);
        end
        Res = A*x - b;
        y = y + beta*Res;
        
        pres = norm(Res);
        
        total_PG = numPG;
        
        call_APG = 1;
        
        while pres > tol % only check primal_res since dual_res automatically below tol
            beta = beta*sig;
            if adp
                [x, numPG] = APG_adp(x0, xk, y, rho, beta, tolk);  
            else
                [x, numPG] = APG(x, xk, y, rho, beta, tolk);
            end
            total_PG = total_PG + numPG;
            Res = A*x - b;
            pres = norm(Res);
            y = y + beta*Res;
            call_APG = call_APG + 1;
        end
        
        fprintf('%d-th subproblem by iALM succeed with %d APG calls and beta = %5.4f,\n',k, call_APG, beta);
    end % of iALM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [x,y,beta,total_PG] = PenMM(xk, yb, rho, tol, beta0)
        
        % =========================================
        % penalty method with estimated multiplier for solving k-subprob:
        %
        %       min f0(x) + rho|x-xk|^2, s.t. Ax == b.
        %
        % where f0(x) = 0.5x'Qx+c'x
        %
        % AL: f0(x) +  rho|x-xk|^2 + <y,Ax-b> + beta/2*|Ax-b|^2
        %
        % the output (x,y) will satisfy
        %           tol/2 - primal feasibility condition
        %           tol/2*min(1,sqrt(rho)) - dual feasibility condition
        %
        % ============================================
        
        x = xk; tolk = (1-frc)*tol/2*min(1, 1/sqrt(rho))*min(1,sqrt(rho));
        beta = beta0;
        
        if adp
            [x, numPG] = APG_adp(x, xk, yb, rho, beta, tolk);
        else
            [x, numPG] = APG(x, xk, yb, rho, beta, tolk);
        end
        Res = A*x - b;
        
        pres = norm(Res);
        
        total_PG = numPG;
        
        call_APG = 1;
        
        while pres > tol % only check primal_res since dual_res automatically below tolk
            beta = beta*sig;
            if adp
                [x, numPG] = APG_adp(x, xk, yb, rho, beta, tolk);
            else
                [x, numPG] = APG(x, xk, yb, rho, beta, tolk);
            end
            total_PG = total_PG + numPG;
            Res = A*x - b;
            pres = norm(Res);
            call_APG = call_APG + 1;
        end
        y = yb + beta*Res;
        
        fprintf('%d-th subproblem by PenMM succeed with %d APG calls and beta = %5.4f,\n',k, call_APG, beta);
        
    end % of PenMM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [x, numPG] = APG(x0, xk, y, rho, beta, tol)
        
        % APG to solve
        %     min f_0(x) + rho|x-xk|^2 + <y,Ax-b> + beta/2*|Ax-b|^2
        %
        % the output x will satisfy
        %       dist( 0, sub_diff(x) ) <= tol
        
        % compuate and use constant Lipschitz constant
        
        HessMat = Q + beta*AtA + 2*rho*eye(n);
        LinVec = c - 2*rho*xk + At*y - beta*At*b;
        LipG = norm(HessMat);
        alpha = sqrt(rho/LipG);       
        extrap_weight = (1 - alpha) / (1 + alpha);
        
        % initialization
        x = x0; xhat = x0;
        
        grad0 = HessMat*x0 + LinVec;
        ngrad = ngrad+1;
        
        gradhat = grad0;
        
        I1 = (x == x_l);
        I2 = (x == x_u);
        I3 = (x > x_l & x < x_u);
        
        grad_res = norm( min(0,grad0(I1)) )^2 + norm( max(0,grad0(I2)) )^2 + norm(grad0(I3))^2;
        grad_res = sqrt(grad_res);
        
        hist_grad_res = grad_res;
        
        numPG = 0;
        
        while grad_res > tol & numPG < APGmaxit
            numPG = numPG + 1;
            x = xhat - gradhat / LipG;
            x = min(x_u,max(x_l,x));
            grad = HessMat*x + LinVec;
            ngrad = ngrad + 1;
            xhat = x + extrap_weight*(x - x0);
            gradhat = grad + extrap_weight*(grad - grad0);
            x0 = x; grad0 = grad;
            I1 = (x == x_l);
            I2 = (x == x_u);
            I3 = (x > x_l & x < x_u);
            
            grad_res = norm( min(0,grad(I1)) )^2 + norm( max(0,grad(I2)) )^2 + norm(grad(I3))^2;
            grad_res = sqrt(grad_res);
            
        end
        
        if numPG == APGmaxit
            fprintf('APG maxiter reached!\n');
        end
        
    end % end of APG.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [x, numPG] = APG_adp(x0, xk, y, rho, beta, tol)
        
        % APG to solve
        %     min f_0(x) + rho|x-xk|^2 + <y,Ax-b> + beta/2*|Ax-b|^2
        %
        % the output x will satisfy
        %       dist( 0, sub_diff(x) ) <= tol
        
        % compuate and use constant Lipschitz constant
        
        HessMat = Q + beta*AtA + 2*rho*eye(n);
        LinVec = c - 2*rho*xk + At*y - beta*At*b;
%        LipG = norm(HessMat); % expensive!!
        
        inc = 2; % tune. 2.
        dec = 2.5; % 2.5. >2 <3
        
        Lip0 = 20 *beta; % tune. 20 ~ 40 good.  >10
        %Lip = Lip0/dec; % *deleted
        
        % initialization
        x = x0; xhat = x0;
        grad0 = HessMat*x0 + LinVec;
        ngrad = ngrad+1;
        gradhat = grad0;
         
        I1 = (x == x_l); % *don't use 'find'!
        I2 = (x == x_u);
        I3 = (x > x_l & x < x_u);
        
        grad_res = norm( min(0,grad0(I1)) )^2 + norm( max(0,grad0(I2)) )^2 + norm(grad0(I3))^2;
        grad_res = sqrt(grad_res);
        
        hist_grad_res = grad_res;
        
        numPG = 0;
        
        
        while grad_res > tol && numPG < APGmaxit
            numPG = numPG + 1;
            objhat = 0.5*(gradhat + LinVec)'*xhat;
            nobj = nobj + 1;
            
            Lip = Lip0/dec;
            obj = inf;
            while obj > objhat + gradhat'*(x - xhat) + 0.5*Lip*norm(x - xhat)^2 + 1e-10 %& Lip < LipG
                %Lip = min(Lip*inc, LipG); % *deleted
                Lip = Lip*inc;
                x = xhat - gradhat / Lip;
                x = min(x_u,max(x_l,x));
                grad = HessMat*x + LinVec;
                ngrad = ngrad + 1;
                numPG = numPG + 1;
                obj = 0.5*(grad + LinVec)'*x;
                nobj = nobj + 1;
            end
            
            alpha = sqrt(rho/Lip);
            extrap_weight = (1 - alpha) / (1 + alpha);
            xhat = x + extrap_weight*(x - x0);
            gradhat = grad + extrap_weight*(grad - grad0); % for QP!
            x0 = x; grad0 = grad; Lip0 = Lip;
            I1 = (x == x_l); I2 = (x == x_u); I3 = (x > x_l & x < x_u);
            grad_res = norm( min(0,grad(I1)) )^2 + norm( max(0,grad(I2)) )^2 + norm(grad(I3))^2;
            grad_res = sqrt(grad_res);
        end
        
        if numPG == APGmaxit
            fprintf('APG_adp maxiter reached!\n');
        end
    end % end of APG_adp.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
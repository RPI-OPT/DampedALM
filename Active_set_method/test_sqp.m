% % min 1/2*x'Qx+c'x s.t. Ax=b; 0 <= x_i <= x_u. % added box for nonconvex!
% Data: c; A,b ~ U[0,1]; Q ~ random PSD matrix with mean = I, covariance = I.
% nonconvex Q case!
clear; close all;
rand('seed',20220927); randn('seed',20220927); % same seed for debug!

m         = 10;
n         = 1000;
num_trial = 10;
rhonum    = 3;
out1      = cell(rhonum,num_trial); 
out2      = cell(rhonum,num_trial);
out_mat   = cell(rhonum,num_trial);
for hardness = 1:rhonum %1:rhonum
    if hardness == 1
        temp = - 0.1;
    end

    if hardness == 2
        temp = - 1;
    end

    if hardness == 3
        temp = - 10;
    end


    for num = 1:num_trial
       
        [U,~]        = qr(randn(n));
        v            = max(0, randn(n,1)*5); % v ~ Gaussian.
        Q            = U*diag(v)*U' + temp*eye(n);
        Q            = (Q + Q') / 2;
        A         = [randn(m, n-m), eye(m)];
        xfeas     = rand(n,1);
        I         = eye(n);
        b         = A*xfeas;
        H         = Q;
        c         = randn(n,1);
        cl        = -5*ones(n,1);
        cu        = 5*ones(n,1);
        lb        = zeros(size(A(1,:)'));
        
        % [x_mat, fval_mat] = quadprog(H, -c, [], [], A, b, cl, cu);
        
        t0 = tic;
        opt = ActiveSet(); %creating an instance of active set method
        [x_sqp, f, iters] = opt.run(H, c, A, b, cu, cl, xfeas);
        t1 = toc(t0);
        
        [xproj, ~] = quadprog(2*eye(n),-2*x_sqp, [], [], A, b, cl, cu);
        realf = 0.5*xproj' * Q * xproj - xproj'*c;

        out1{hardness, num}.xproj = xproj;
        out1{hardness, num}.realf = realf;
        out1{hardness, num}.fval = f;
        out1{hardness, num}.x    = x_sqp;
        out1{hardness, num}.iter = iters;
        out1{hardness, num}.time = t1;

        %% Compare to iPALM
        QP     = [];
        QP.A   = A;
        QP.b   = b;
        QP.Q   = Q;
        QP.r   = -c;
        QP.ell = -5;
        QP.u   = 5;
        
        if hardness == 3
            beta0 = 5;
        elseif hardness == 2
            beta0 = 0.001;
        elseif hardness == 1
            beta0 = 0.0001;
        end

        rho    = -1*temp;
        epsilon = 1e-3;
        maxouter = 1e3;
        mu      = rho;
        y       = zeros(size(A(:,1)));
        
        % running DPALM and comparing with sqp
        t0 = tic;
        [xdp,~,out2{hardness, num}] = iPALM_QP(QP,xfeas, y,I, beta0,rho,epsilon,maxouter,mu);
        time2 = toc(t0);
             
        [xproj, ~] = quadprog(2*eye(n),-2*xdp, [], [], QP.A, QP.b, cl, cu);
        realf = 0.5*xproj' * QP.Q * xproj + xproj'*QP.r;

        out2{hardness,num}.time = time2;
        out2{hardness, num}.xproj = xproj;
        out2{hardness, num}.realf = realf;
        
        t0 = tic;
        [x_mat, fval_mat] = quadprog(H, -c, [], [], A, b, cl, cu);
        time_mat = toc(t0);
        out_mat{hardness, num}.time = time_mat;
        out_mat{hardness, num}.fval = fval_mat;
        out_mat{hardness, num}.x = x_mat;


    end
    
    save result_sqp_dpalm out1 out2 out_mat m n num_trial
end
%%

num = num_trial;
varianceval = cell(3,1);
meanvalues = cell(3,1);
for i = 1:3
    
    objavg1 = zeros(1,num);
    timeavg1 = zeros(1,num);

    objavg2 = zeros(1,num);
    timeavg2 = zeros(1,num);

    objavg3 = zeros(1,num);
    timeavg3 = zeros(1,num);

    for j = 1:num
        [objavg1(j),k]   = max(out1{i,j}.realf);
        timeavg1(j)      = out1{i,j}.time;
     
        [objavg2(j),k]   = max(out2{i,j}.realf);
        timeavg2(j)      = out2{i,j}.time;
        
        [objavg3(j),k]   = max(out_mat{i,j}.fval);
        timeavg3(j)      = out_mat{i,j}.time;
        
    end
    [varianceval{i}.obj_sqp,meanvalues{i}.obj_sqp] = var(objavg1);
    [varianceval{i}.time_sqp, meanvalues{i}.time_sqp] = var(timeavg1);
    
    [varianceval{i}.obj_dpalm,meanvalues{i}.obj_dpalm] = var(objavg2);
    [varianceval{i}.time_dpalm, meanvalues{i}.time_dpalm] = var(timeavg2);

    [varianceval{i}.obj_mat,meanvalues{i}.obj_mat] = var(objavg3);
    [varianceval{i}.time_mat, meanvalues{i}.time_mat] = var(timeavg3);

end

%%

count = zeros(3,1);
for i = 1:num_trial
    for j = 1:3
        if out1{j,i}.realf > out2{j,i}.realf
            count(j) = count(j) + 1;
        end
    end
end

%%
clear; close all
rand('seed',20220927); randn('seed',20220927); % same seed for debug!

m         = 10;
n         = 1000;
num_trial = 10;
rhonum    = 1;
out3      = cell(rhonum,num_trial); 
out4      = cell(rhonum,num_trial);
out_mat_cvx = cell(rhonum, num_trial);

for hardness = 1:rhonum %1:rhonum
    if hardness == 1
        temp = 1e-3;
    end

    for num = 1:num_trial
       
        [U,~]        = qr(randn(n));
        v            = max(0, randn(n,1)*5); % v ~ Gaussian.
        Q            = U*diag(v)*U' + temp*eye(n);
        Q            = (Q + Q') / 2;
        A         = [randn(m, n-m), eye(m)];
        xfeas     = rand(n,1);
        I         = eye(n);
        b         = A*xfeas;
        H         = Q;
        c         = randn(n,1);
        cl        = -5*ones(n,1);
        cu        = 5*ones(n,1);
        lb        = zeros(size(A(1,:)'));
        
        t0 = tic;
        opt = ActiveSet(); %creating an instance of active set method
        [x_sqp, f, iters] = opt.run(H, c, A, b, cu, cl, xfeas);
        t1 = toc(t0);
        
        [xproj, ~] = quadprog(2*eye(n),-2*x_sqp, [], [], A, b, cl, cu);
        realf = 0.5*xproj' * Q * xproj - xproj'*c;

        out3{hardness, num}.xproj = xproj;
        out3{hardness, num}.realf = realf;
        out3{hardness, num}.fval = f;
        out3{hardness, num}.x    = x_sqp;
        out3{hardness, num}.iter = iters;
        out3{hardness, num}.time = t1;

        %% Compare to iPALM
        QP     = [];
        QP.A   = A;
        QP.b   = b;
        QP.Q   = Q;
        QP.r   = -c;
        QP.ell = -5;
        QP.u   = 5;
        
        
        beta0 = 0.01;

        rho    = abs(-1*temp);
        epsilon = 1e-3;
        maxouter = 1e3;
        mu      = rho;
        y       = zeros(size(A(:,1)));
        
        % running DPALM and comparing with sqp
        t0 = tic;
        [xdp,~,out4{hardness, num}] = iPALM_QP(QP,xfeas, y,I, beta0,rho,epsilon,maxouter,mu);
        time2 = toc(t0);
             
        [xproj, ~] = quadprog(2*eye(n),-2*xdp, [], [], QP.A, QP.b, cl, cu);
        realf = 0.5*xproj' * QP.Q * xproj + xproj'*QP.r;

        out4{hardness,num}.time = time2;
        out4{hardness, num}.xproj = xproj;
        out4{hardness, num}.realf = realf;

        t0 = tic;
        [x_mat, fval_mat] = quadprog(H, -c, [], [], A, b, cl, cu);
        time_mat = toc(t0);
        out_mat_cvx{hardness, num}.time = time_mat;
        out_mat_cvx{hardness, num}.fval = fval_mat;
        out_mat_cvx{hardness, num}.x = x_mat;
       

    end
    
    save result_sqp_convex out3 out4 out_mat_cvx m n num_trial
end

%%

num = num_trial;
varianceval = cell(1,1);
meanvalues = cell(1,1);

objavg1 = zeros(1,num);
timeavg1 = zeros(1,num);

objavg2 = zeros(1,num);
timeavg2      = zeros(1,num);

objavg3 = zeros(1, num);
timeavg3 = zeros(1, num);

for j = 1:num
    [objavg1(j),k]   = max(out3{1,j}.realf);
    timeavg1(j)      = out3{1,j}.time;
 
    [objavg2(j),k]   = max(out4{1,j}.realf);
    timeavg2(j)      = out4{1,j}.time;

    [objavg3(j),k]   = max(out_mat_cvx{1,j}.fval);
    timeavg3(j)      = out_mat_cvx{1,j}.time;
    
end
[varianceval{1}.obj_sqp,meanvalues{1}.obj_sqp] = var(objavg1);
[varianceval{1}.time_sqp, meanvalues{1}.time_sqp] = var(timeavg1);

[varianceval{1}.obj_dpalm,meanvalues{1}.obj_dpalm] = var(objavg2);
[varianceval{1}.time_dpalm, meanvalues{1}.time_dpalm] = var(timeavg2);

[varianceval{1}.obj_mat,meanvalues{1}.obj_mat] = var(objavg3);
[varianceval{1}.time_mat, meanvalues{1}.time_mat] = var(timeavg3);
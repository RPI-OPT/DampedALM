%%
clear; close all;
rng(20220927,"twister");

m            = 10;
n            = 1000;
num_trial    = 10;
rhonum       = 3;
tol          = 1e-3;
maxouteriter = 1e4;
out1         = cell(rhonum,num_trial);
out2         = cell(rhonum,num_trial);
out3         = cell(rhonum,num_trial);
out4         = cell(rhonum,num_trial);
out5         = cell(rhonum,num_trial);

x_l = -5 * ones(n,1); % lower box. -5. 0 for aipp to work.
x_u = 5 * ones(n,1); % upper box. 5

for hardness = 1:rhonum
   
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
        [U,~] = qr(randn(n));
        v     = max(0, randn(n,1)*5); % v ~ Gaussian.
        Q     = U*diag(v)*U' + temp*eye(n);
        Q     = (Q + Q') / 2; % symmetrization is important
        A     = [randn(m, n-m), eye(m)];
        b     = randn(m,1) + 0.1;
        I     = eye(n);
        r     = randn(n,1);
        y0    = rand(m,1);
        x0    = quadprog(eye(n),randn(n,1),[],[],A,b,x_l,x_u);
%%
        QP     = [];
        QP.A   = A;
        QP.b   = b;
        QP.Q   = Q;
        QP.r   = r;
        QP.ell = -5;
        QP.u   = 5;

        rho    = -1*temp;

        if hardness == 3
            beta0 = 5;
        elseif hardness == 2
            beta0 = 0.0001;
        elseif hardness == 1
            beta0 = 0.001;
        end

        mu = rho;
        t0 = tic;
        [~,~,out1{hardness, num}] = iPALM_QP(QP,x0, y0,I,...
                      beta0,rho,tol,maxouteriter,mu);
        time1 = toc(t0);
        out1{hardness,num}.time = time1;
        
        %% Compare to NLNL-IAPIAL

        inp     = [];
        inp.A   = A;
        inp.b   = b;
        inp.Q   = Q;
        inp.r   = r;
        inp.ell = -5;
        inp.u   = 5;
        mf      = -1*min(eig(Q));
        sigma   = 1/sqrt(2);
        lambda  = 0.5/(-1*temp);
        if hardness == 3
            tol1 = 1e-5;
            beta0 = 100;
        elseif hardness == 2
            tol1 = 1;
            beta0 = 0.0001;
        elseif hardness == 1
            tol1 = 1;
            beta0 = 0.001;

        end
        t0     = tic;
        [~,out2{hardness,num}] = NLIPIAL_QP(inp,x0,y0,I,mf,...
                    sigma, beta0, tol,maxouteriter,lambda, tol1);
        time2  = toc(t0);
        
        out2{hardness,num}.time = time2;
        
        if hardness == 3
            tol1 = 5e2;
        end
        %% compare to IPLA

        t0              = tic;
        [~,out3{hardness,num}] = IPLA_QP(inp,x0,y0,I,mf,...
                    sigma, beta0, tol,maxouteriter,lambda, tol1);
        time3           = toc(t0);

        out3{hardness,num}.time = time3;

        %Testing LiMEAL 
        NormQ  = norm(QP.Q);
        z0     = x0;
        lamb0  = zeros(m,1);
        gamma  = 0.5/NormQ;
        if hardness == 3
            beta   = 100;
        else
            beta = -1*temp;
        end

        NumIter = 1e4;
        epsilon = tol;
        eta     = 1.5;

        t0 = tic;
        [~,~,~,out5{hardness,num}] = LiMEAL_QP_theory(QP,x0,z0,...
                   lamb0,beta,gamma,eta,NumIter,epsilon,I);
        time4 = toc(t0);

        out5{hardness,num}.time = time4;

        %%many iALM, with constant Lip constant
        opts          = [];
        opts.tol      = tol;
        opts.maxit    = maxouteriter;
        opts.K0       = 1e4;
        opts.K1       = 2;
        opts.sig      = 3;
        opts.gamma    = 1.1;
        opts.x0       = x0;
        opts.beta0    = 0.01;
        opts.APGmaxit = 1000000;
        opts.adp      = 0;
        opts.frc      = 0.9;
        
        % Before calling HiAPeM_qp, save current globals
        saved_Q  = Q;
        saved_r  = r;
        saved_A  = A;
        saved_b  = b;
        saved_xl = x_l;
        saved_xu = x_u;
        t0    = tic;
        [~,~,out4{hardness,num}] = HiAPeM_qp(Q,r,A,b,x_l,x_u,opts);
        time1 = toc(t0);

        out4{hardness,num}.time = time1;
        fprintf('\n\n');        

        % Restore globals to original values
        Q   = saved_Q;
        r   = saved_r;
        A   = saved_A;
        b   = saved_b;
        x_l = saved_xl;
        x_u = saved_xu;

    end
end

%%
save result_ncLCQP out1 out2 out3 out4 out5 num_trial rhonum
%%

gradnum1 = zeros(3,1);
pres1    = zeros(3,1);
dres1    = zeros(3,1);
time1    = zeros(3,1);

varianceval = cell(3,1);
meanvalues  = cell(3,1);
for i = 1:rhonum
    gradnumavg1   = zeros(1,num_trial);
    presavg1      = zeros(1,num_trial);
    dualresavg1   = zeros(1,num_trial);
    timeavg1      = zeros(1,num_trial);

    gradnumavg2   = zeros(1,num_trial);
    presavg2      = zeros(1,num_trial);
    dualresavg2   = zeros(1,num_trial);
    timeavg2      = zeros(1,num_trial);

    gradnumavg3   = zeros(1,num_trial);
    presavg3      = zeros(1,num_trial);
    dualresavg3   = zeros(1,num_trial);
    timeavg3      = zeros(1,num_trial);

    gradnumavg4   = zeros(1,num_trial);
    presavg4      = zeros(1,num_trial);
    dualresavg4   = zeros(1,num_trial);
    timeavg4      = zeros(1,num_trial);

    gradnumavg5   = zeros(1,num_trial);
    presavg5      = zeros(1,num_trial);
    dualresavg5   = zeros(1,num_trial);
    timeavg5      = zeros(1,num_trial);

    for j = 1:num_trial
        [gradnumavg1(j),k]   = max(out1{i,j}.gradnum);
        presavg1(j)          = out1{i,j}.primalviol(k);
        dualresavg1(j)       = out1{i,j}.dualviol(k);
        timeavg1(j)          = out1{i,j}.time;

        [gradnumavg2(j),k]   = max(out2{i,j}.gradnum);
        presavg2(j)          = out2{i,j}.primalviol(k);
        dualresavg2(j)       = out2{i,j}.dualviol(k);
        timeavg2(j)          = out2{i,j}.time;

        [gradnumavg3(j),k]   = max(out3{i,j}.gradnum);
        presavg3(j)          = out3{i,j}.primalviol(k);
        dualresavg3(j)       = out3{i,j}.dualviol(k);
        timeavg3(j)          = out3{i,j}.time;

        [gradnumavg4(j),k]   = max(out4{i,j}.gradnum);
        presavg4(j)          = out4{i,j}.primalviol(k);
        dualresavg4(j)       = out4{i,j}.dualviol(k);
        timeavg4(j)          = out4{i,j}.time;
        
        [gradnumavg5(j),k]   = max(out5{i,j}.gradnum);
        presavg5(j)          = out5{i,j}.primalviol(k);
        dualresavg5(j)       = out5{i,j}.dualviol(k);
        timeavg5(j)          = out5{i,j}.time;

    end
    [varianceval{i}.gradnum1,meanvalues{i}.gradnum1] = var(gradnumavg1);
    [varianceval{i}.pres1, meanvalues{i}.pres1]      = var(presavg1);
    [varianceval{i}.dres1, meanvalues{i}.dres1]      = var(dualresavg1);
    [varianceval{i}.time1, meanvalues{i}.time1]      = var(timeavg1);

    [varianceval{i}.gradnum2,meanvalues{i}.gradnum2] = var(gradnumavg2);
    [varianceval{i}.pres2, meanvalues{i}.pres2]      = var(presavg2);
    [varianceval{i}.dres2, meanvalues{i}.dres2]      = var(dualresavg2);
    [varianceval{i}.time2, meanvalues{i}.time2]      = var(timeavg2);

    [varianceval{i}.gradnum3,meanvalues{i}.gradnum3] = var(gradnumavg3);
    [varianceval{i}.pres3, meanvalues{i}.pres3]      = var(presavg3);
    [varianceval{i}.dres3, meanvalues{i}.dres3]      = var(dualresavg3);
    [varianceval{i}.time3, meanvalues{i}.time3]      = var(timeavg3);
    
    [varianceval{i}.gradnum4,meanvalues{i}.gradnum4] = var(gradnumavg4);
    [varianceval{i}.pres4, meanvalues{i}.pres4]      = var(presavg4);
    [varianceval{i}.dres4, meanvalues{i}.dres4]      = var(dualresavg4);
    [varianceval{i}.time4, meanvalues{i}.time4]      = var(timeavg4);

    [varianceval{i}.gradnum5,meanvalues{i}.gradnum5] = var(gradnumavg5);
    [varianceval{i}.pres5, meanvalues{i}.pres5]      = var(presavg5);
    [varianceval{i}.dres5, meanvalues{i}.dres5]      = var(dualresavg5);
    [varianceval{i}.time5, meanvalues{i}.time5]      = var(timeavg5);

end
%%

rhos = [0.1, 1, 10];
% Open the file for writing
fileID = fopen('Tables/mean_and_variance_qp.txt', 'w');

fprintf(fileID, '%s\n', '                   DPALM_pres     DPALM_dres   DPALM_time   DPALM_gradnum    NLIAPIAL_pres   NLIAPIAL_dres   NLIAPIAL_time   NLIAPIAL_gradnum      IPLA_pres      IPLA_dres      IPLA_time       IPLA_gradnum     HiAPeM_pres     HiAPeM_dres      HiAPeM_time   HiAPeM_gradnum   LiMEAL_pres   LiMEAL_dres   LiMEAL_time   LiMEAL_gradnum');
fprintf(fileID, '%s\n', ' ');
fprintf(fileID, '%s\n', ['mean for rho = .1   ',...
        num2str(meanvalues{1}.pres1,3), '       ', num2str(meanvalues{1}.dres1,3), '       ', num2str(meanvalues{1}.time1,3), '          ', num2str(meanvalues{1}.gradnum1,3), '        '...
        num2str(meanvalues{1}.pres2,3), '         ' num2str(meanvalues{1}.dres2,3), '          ', num2str(meanvalues{1}.time2,3), '           ', num2str(meanvalues{1}.gradnum2,3), '           '...
        num2str(meanvalues{1}.pres3,3), '       ' num2str(meanvalues{1}.dres3,3), '         ', num2str(meanvalues{1}.time3,3), '           ', num2str(meanvalues{1}.gradnum3,3), '         '...
        num2str(meanvalues{1}.pres4,3), '        ' num2str(meanvalues{1}.dres4,3), '            ', num2str(meanvalues{1}.time4,3), '        ', num2str(meanvalues{1}.gradnum4,3), '         '...
         num2str(meanvalues{1}.pres5,3), '      ' num2str(meanvalues{1}.dres5,3), '         ', num2str(meanvalues{1}.time5,3), '           ', num2str(meanvalues{1}.gradnum5,3), '       ']);

fprintf(fileID, '%s\n', ['var for rho = .1    ',...
        num2str(varianceval{1}.pres1,3), '       ', num2str(varianceval{1}.dres1,3), '       ', num2str(varianceval{1}.time1,3), '          ', num2str(varianceval{1}.gradnum1,3), '        '...
        num2str(varianceval{1}.pres2,3), '        ' num2str(varianceval{1}.dres2,3), '          ', num2str(varianceval{1}.time2,3), '           ', num2str(varianceval{1}.gradnum2,3), '           '...
        num2str(varianceval{1}.pres3,3), '       ' num2str(varianceval{1}.dres3,3), '         ', num2str(varianceval{1}.time3,3), '           ', num2str(varianceval{1}.gradnum3,3), '         '...
        num2str(varianceval{1}.pres4,3), '        ' num2str(varianceval{1}.dres4,3), '         ', num2str(varianceval{1}.time4,3), '     ', num2str(varianceval{1}.gradnum4,3), '         '...
        num2str(varianceval{1}.pres5,3), '      ' num2str(varianceval{1}.dres5,3), '        ', num2str(varianceval{1}.time5,3), '           ', num2str(varianceval{1}.gradnum5,3), '       ']);

fprintf(fileID, '%s\n', ['mean for rho = 1    ',...
        num2str(meanvalues{2}.pres1,3), '       ', num2str(meanvalues{2}.dres1,3), '       ', num2str(meanvalues{2}.time1,3), '          ', num2str(meanvalues{2}.gradnum1,3), '        '...
        num2str(meanvalues{2}.pres2,3), '        ' num2str(meanvalues{2}.dres2,3), '          ', num2str(meanvalues{2}.time2,3), '           ', num2str(meanvalues{2}.gradnum2,3), '           '...
        num2str(meanvalues{2}.pres3,3), '       ' num2str(meanvalues{2}.dres3,3), '         ', num2str(meanvalues{2}.time3,3), '           ', num2str(meanvalues{2}.gradnum3,3), '         '...
        num2str(meanvalues{2}.pres4,3), '        ' num2str(meanvalues{2}.dres4,3), '           ', num2str(meanvalues{2}.time4,3), '        ', num2str(meanvalues{2}.gradnum4,3), '         '...
         num2str(meanvalues{2}.pres5,3), '      ' num2str(meanvalues{2}.dres5,3), '       ', num2str(meanvalues{2}.time5,3), '            ', num2str(meanvalues{2}.gradnum5,3), '       ']);

fprintf(fileID, '%s\n', ['var for rho = 1     ',...
        num2str(varianceval{2}.pres1,3), '       ', num2str(varianceval{2}.dres1,3), '       ', num2str(varianceval{2}.time1,3), '          ', num2str(varianceval{2}.gradnum1,3), '        '...
        num2str(varianceval{2}.pres2,3), '        ' num2str(varianceval{2}.dres2,3), '          ', num2str(varianceval{2}.time2,3), '           ', num2str(varianceval{2}.gradnum2,3), '           '...
        num2str(varianceval{2}.pres3,3), '        ' num2str(varianceval{2}.dres3,3), '         ', num2str(varianceval{2}.time3,3), '            ', num2str(varianceval{2}.gradnum3,3), '         '...
        num2str(varianceval{2}.pres4,3), '        ' num2str(varianceval{2}.dres4,3), '          ', num2str(varianceval{2}.time4,3), '     ', num2str(varianceval{2}.gradnum4,3), '          '...
        num2str(varianceval{2}.pres5,3), '      ' num2str(varianceval{2}.dres5,3), '       ', num2str(varianceval{2}.time5,3), '       ', num2str(varianceval{2}.gradnum5,3), '      ']);

fprintf(fileID, '%s\n', ['mean for rho = 10   ',...
        num2str(meanvalues{3}.pres1,3), '       ', num2str(meanvalues{3}.dres1,3), '        ', num2str(meanvalues{3}.time1,3), '          ', num2str(meanvalues{3}.gradnum1,3), '        '...
        num2str(meanvalues{3}.pres2,3), '        ' num2str(meanvalues{3}.dres2,3), '          ', num2str(meanvalues{3}.time2,3), '           ', num2str(meanvalues{3}.gradnum2,3), '           '...
        num2str(meanvalues{3}.pres3,3), '       ' num2str(meanvalues{3}.dres3,3), '         ', num2str(meanvalues{3}.time3,3), '           ', num2str(meanvalues{3}.gradnum3,3), '         '...
        num2str(meanvalues{3}.pres4,3), '        ' num2str(meanvalues{3}.dres4,3), '          ', num2str(meanvalues{3}.time4,3), '        ', num2str(meanvalues{3}.gradnum4,3), '         '...
        num2str(meanvalues{3}.pres5,3), '       ' num2str(meanvalues{3}.dres5,3), '        ', num2str(meanvalues{3}.time5,3), '           ', num2str(meanvalues{3}.gradnum5,3), '       ']);

fprintf(fileID, '%s\n', ['var for rho = 1     ',...
        num2str(varianceval{3}.pres1,3), '       ', num2str(varianceval{3}.dres1,3), '       ', num2str(varianceval{3}.time1,3), '          ', num2str(varianceval{3}.gradnum1,3), '        '...
        num2str(varianceval{3}.pres2,3), '        ' num2str(varianceval{3}.dres2,3), '          ', num2str(varianceval{3}.time2,3), '           ', num2str(varianceval{3}.gradnum2,3), '           '...
        num2str(varianceval{3}.pres3,3), '       ' num2str(varianceval{3}.dres3,3), '          ', num2str(varianceval{3}.time3,3), '          ', num2str(varianceval{3}.gradnum3,3), '         '...
        num2str(varianceval{3}.pres4,3), '        ' num2str(varianceval{3}.dres4,3), '          ', num2str(varianceval{3}.time4,3), '         ', num2str(varianceval{3}.gradnum4,3), '         '...
        num2str(varianceval{3}.pres5,3), '      ' num2str(varianceval{3}.dres5,3), '       ', num2str(varianceval{3}.time5,3), '           ', num2str(varianceval{3}.gradnum5,3), '       ']);
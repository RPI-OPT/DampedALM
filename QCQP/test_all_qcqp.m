clear; close all % Test code for nc_QCQP
%% generate data
rand('seed',20230119); randn('seed',20230119);
n = 1000; % 1000, 100
m = 10; % 10

num = 10;
rho_val = [.1,1,10];
out1 = cell(length(rho_val),num);
out2 = cell(length(rho_val),num);
out3 = cell(length(rho_val),num);
out4 = cell(length(rho_val),num);

% load result_ncQCQP_v1
for rho_num = 1:length(rho_val)
    rho = rho_val(rho_num);
    
    for ii = 1:num
        
        Q = cell(1,m+1);
        c = cell(1,m+1);
        d = zeros(m+1,1);
        
        [U,~] = qr(randn(n));
        
        v = max(0, randn(n,1)*5); % v ~ Gassian.
        Q{m+1} = U*diag(v)*U' - rho*eye(n);
        Q{m+1} = (Q{m+1} + Q{m+1}') / 2; % symmetrization is important
        c{m+1} = randn(n,1);
        
        fprintf('weak convexity constant is %f\n', min(eig(Q{m+1})));
        
        %%
        for j = 1:m
            [Q{j},~] = qr(randn(n));
            Q{j} = Q{j}(:, 1:n-5);
            v = rand(n-5,1)*5 + 1; % can remove +1 here...
            Q{j} = Q{j}*diag(v)*Q{j}';
            Q{j} = (Q{j} + Q{j}') / 2;
            c{j} = randn(n,1)*5;
            d(j) = min(0, randn*2) - 0.1;
        end
        
        lb = -5*ones(n,1); % -5, -1. can change.
        ub = 5*ones(n,1); % 5, 1. can change.
        
        maxsubit = 10000;
        tol = 1e-3;
        
        %%
        opts = [];
        opts.x0 = zeros(n,1);
        opts.maxsubit = maxsubit;
        opts.maxit = 10000;
        opts.sig = 3;
        opts.beta0 = 0.01;
        opts.tol = tol;
        opts.ver = 3; % opts.ver = 3; %% 1 for const beta
        opts.inc = 2;
        opts.dec= 2.5;
        opts.K0 = 10000;    
        
        t1 = tic;
        [x1,z,out1{rho_num,ii}] = HiAPeM_qcqp(Q,c,d,m,lb,ub,opts);
        time1 = toc(t1); out1{rho_num,ii}.time = time1;
        fprintf('HiAPeM time = %f\n', time1);

        %%
        opts1          = [];
        opts1.x0       = zeros(n,1);
        opts1.z0      = zeros(m,1);
        opts1.maxsubit = maxsubit;
        opts1.maxit    = 10000;
        opts1.beta0    = 0.0001;
        opts1.tol      = tol;
        opts1.rho      = -2*min(eig(Q{m+1}));
        opts1.gammau   = 1.5;
        opts1.gammad   = 2;
        opts1.L        = norm(Q{m+1});
        down           = -5;
        up             = 5;

        Q0             = Q{m+1};
        c0             = c{m+1};
        Q1             = Q;
        c1             = c;
        Q1(end)         = []; 
        c1(end)         = [];

        t2 = tic;
        [out2{rho_num, ii}] = iPALM_qcqp(Q0, Q1, c0, c1,d,down, up, opts1);
        time2 = toc(t2);
        out2{rho_num, ii}.time  = time2;
        fprintf("iPALM time = %f\n", time2)
%%      
        opts1.lambda   = 0.5/(-1.0*min(eig(Q{m+1})));
        sigma          = 0.5/sqrt(2);
        if rho_num == 3
            opts1.beta = 0.001;
        else
            opts1.beta     = .0002;
        end
        t3 = tic;
        [out3{rho_num, ii}] = IAPIAL_qcqp(Q0, Q, c0, c,d,down, up,sigma, opts1);
        time3 = toc(t3);
        out3{rho_num,ii}.time  = time3;
        fprintf("IAPIAL time = %f\n", time3)

        [out4{rho_num, ii}] = IPLA_qcqp(Q0, Q, c0, c,d,down, up,sigma, opts1);
        time4 = toc(t3);
        out4{rho_num,ii}.time  = time4;
        fprintf("IPLA time = %f\n", time4)

    end
end
%%
save result_ncQCQP out1 out2 out3 out4 m n
%%
num=10;
gradnum1 = zeros(4,1);
pres1    = zeros(4,1);
dres1    = zeros(4,1);
compslack1 = zeros(4,1);
time1      = zeros(4,1);

varianceval = cell(3,1);
meanvalues = cell(3,1);
for i = 1:3
    gradnumavg1   = zeros(1,num);
    presavg1      = zeros(1,num);
    dualresavg1   = zeros(1,num);
    compslackavg1 = zeros(1,num);
    timeavg1      = zeros(1,num);

    gradnumavg2   = zeros(1,num);
    presavg2      = zeros(1,num);
    dualresavg2   = zeros(1,num);
    compslackavg2 = zeros(1,num);
    timeavg2      = zeros(1,num);

    gradnumavg3   = zeros(1,num);
    presavg3      = zeros(1,num);
    dualresavg3   = zeros(1,num);
    compslackavg3 = zeros(1,num);
    timeavg3      = zeros(1,num);
    
    gradnumavg4   = zeros(1,num);
    presavg4      = zeros(1,num);
    dualresavg4   = zeros(1,num);
    compslackavg4 = zeros(1,num);
    timeavg4      = zeros(1,num);
    for j = 1:num
        [gradnumavg1(j),k]   = max(out1{i,j}.gradnum);
        presavg1(j)          = out1{i,j}.primalviol(k);
        dualresavg1(j)       = out1{i,j}.dualviol(k);
        compslackavg1(j)     = out1{i,j}.compslack(k);
        timeavg1(j)          = out1{i,j}.time;

        [gradnumavg2(j),k]   = max(out2{i,j}.gradnum);
        presavg2(j)          = out2{i,j}.primalviol(k);
        dualresavg2(j)       = out2{i,j}.dualviol(k);
        compslackavg2(j)     = out2{i,j}.compslack(k);
        timeavg2(j)          = out2{i,j}.time;

        [gradnumavg3(j),k]   = max(out3{i,j}.gradnum);
        presavg3(j)          = out3{i,j}.primalviol(k);
        dualresavg3(j)       = out3{i,j}.dualviol(k);
        compslackavg3(j)     = out3{i,j}.compslack(k);
        timeavg3(j)          = out3{i,j}.time;

        [gradnumavg4(j),k]   = max(out4{i,j}.gradnum);
        presavg4(j)          = out4{i,j}.primalviol(k);
        dualresavg4(j)       = out4{i,j}.dualviol(k);
        compslackavg4(j)     = out4{i,j}.compslack(k);
        timeavg4(j)          = out4{i,j}.time;
    end
    [varianceval{i}.gradnum1,meanvalues{i}.gradnum1] = var(gradnumavg1);
    [varianceval{i}.pres1, meanvalues{i}.pres1] = var(presavg1);
    [varianceval{i}.dres1, meanvalues{i}.dres1] = var(dualresavg1);
    [varianceval{i}.compslack1, meanvalues{i}.compslack1] = var(compslackavg1);
    [varianceval{i}.time1, meanvalues{i}.time1] = var(timeavg1);

    [varianceval{i}.gradnum2,meanvalues{i}.gradnum2] = var(gradnumavg2);
    [varianceval{i}.pres2, meanvalues{i}.pres2] = var(presavg2);
    [varianceval{i}.dres2, meanvalues{i}.dres2] = var(dualresavg2);
    [varianceval{i}.compslack2, meanvalues{i}.compslack2] = var(compslackavg2);
    [varianceval{i}.time2, meanvalues{i}.time2] = var(timeavg2);

    [varianceval{i}.gradnum3,meanvalues{i}.gradnum3] = var(gradnumavg3);
    [varianceval{i}.pres3, meanvalues{i}.pres3] = var(presavg3);
    [varianceval{i}.dres3, meanvalues{i}.dres3] = var(dualresavg3);
    [varianceval{i}.compslack3, meanvalues{i}.compslack3] = var(compslackavg3);
    [varianceval{i}.time3, meanvalues{i}.time3] = var(timeavg3);

    [varianceval{i}.gradnum4,meanvalues{i}.gradnum4] = var(gradnumavg4);
    [varianceval{i}.pres4, meanvalues{i}.pres4] = var(presavg4);
    [varianceval{i}.dres4, meanvalues{i}.dres4] = var(dualresavg4);
    [varianceval{i}.compslack4, meanvalues{i}.compslack4] = var(compslackavg4);
    [varianceval{i}.time4, meanvalues{i}.time4] = var(timeavg4);
end

%%

rhos = [0.1, 1, 10];
% Open the file for writing
fileID = fopen('Tables/mean_and_variance_qcqp.txt', 'w');

fprintf(fileID, '%s\n', '                     HiAPeM_pres   HiAPeM_dres   HiAPeM_compslack   HiAPeM_time   HiAPeM_gradnum   NLIAPIAL_pres   NLIAPIAL_dres NLIAPIAL_compslack   NLIAPIAL_time   NLIAPIAL_gradnum   DPALM_pres       DPALM_dres  DPALM_compslack DPALM_time   DPALM_gradnum      IPLA_pres       IPLA_dres    IPLA_compslack      IPLA_time   IPLA_gradnum');
fprintf(fileID, '%s\n', ' ');
 
fprintf(fileID, '%s\n', ['mean for rho = .1     ', num2str(meanvalues{1,1}.pres1,3), '       ', num2str(meanvalues{1,1}.dres1,3), '       ', num2str(meanvalues{1,1}.compslack1,3), '           ', num2str(meanvalues{1,1}.time1,3), '          ', num2str(meanvalues{1,1}.gradnum1,3), '        '...
        num2str(meanvalues{1,1}.pres3,3), '        ', num2str(meanvalues{1,1}.dres3,3), '          ', num2str(meanvalues{1,1}.compslack3,3), '            ', num2str(meanvalues{1,1}.time3,3), '           ', num2str(meanvalues{1,1}.gradnum3,3), '         ',...
        num2str(meanvalues{1,1}.pres2,3), '        ', num2str(meanvalues{1,1}.dres2,3), '       ', num2str(meanvalues{1,1}.compslack2,3), '       ', num2str(meanvalues{1,1}.time2,3), '         ', num2str(meanvalues{1,1}.gradnum2,3), '         '...
        num2str(meanvalues{1,1}.pres4,3), '        ', num2str(meanvalues{1,1}.dres4,3), '       ', num2str(meanvalues{1,1}.compslack4,3), '             ', num2str(meanvalues{1,1}.time4,3), '         ', num2str(meanvalues{1,1}.gradnum4,3), '         ']);
fprintf(fileID, '%s\n', ['mean for rho = .1     ', num2str(varianceval{1,1}.pres1,3), '       ', num2str(varianceval{1,1}.dres1,3), '       ', num2str(varianceval{1,1}.compslack1,3), '           ', num2str(varianceval{1,1}.time1,3), '          ', num2str(varianceval{1,1}.gradnum1,3), '        '...
        num2str(varianceval{1,1}.pres3,3), '         ', num2str(varianceval{1,1}.dres3,3), '          ', num2str(varianceval{1,1}.compslack3,3), '            ', num2str(varianceval{1,1}.time3,3), '           ', num2str(varianceval{1,1}.gradnum3,3), '         '...
        num2str(varianceval{1,1}.pres2,3), '        ', num2str(varianceval{1,1}.dres2,3), '       ', num2str(varianceval{1,1}.compslack2,3), '       ', num2str(varianceval{1,1}.time2,3), '         ', num2str(varianceval{1,1}.gradnum2,3), '         ',...
        num2str(varianceval{1,1}.pres4,3), '        ', num2str(varianceval{1,1}.dres4,3), '       ', num2str(varianceval{1,1}.compslack4,3), '             ', num2str(varianceval{1,1}.time4,3), '        ', num2str(varianceval{1,1}.gradnum4,3), '         ']);

fprintf(fileID, '%s\n', ['mean for rho = 1      ', num2str(meanvalues{2,1}.pres1,3), '       ', num2str(meanvalues{2,1}.dres1,3), '         ', num2str(meanvalues{2,1}.compslack1,3), '           ', num2str(meanvalues{2,1}.time1,3), '           ', num2str(meanvalues{2,1}.gradnum1,3), '        '...
        num2str(meanvalues{2,1}.pres3,3), '        ', num2str(meanvalues{2,1}.dres3,3), '          ', num2str(meanvalues{2,1}.compslack3,3), '            ', num2str(meanvalues{2,1}.time3,3), '             ', num2str(meanvalues{2,1}.gradnum3,3),'         '...
        num2str(meanvalues{2,1}.pres2,3), '        ', num2str(meanvalues{2,1}.dres2,3), '       ', num2str(meanvalues{2,1}.compslack2,3), '       ', num2str(meanvalues{2,1}.time2,3), '         ', num2str(meanvalues{2,1}.gradnum2,3), '         ',...
        num2str(meanvalues{2,1}.pres4,3), '        ', num2str(meanvalues{2,1}.dres4,3), '        ', num2str(meanvalues{2,1}.compslack4,3), '                ', num2str(meanvalues{2,1}.time4,3), '         ', num2str(meanvalues{2,1}.gradnum4,3),]);
fprintf(fileID, '%s\n', ['mean for rho = 1      ', num2str(varianceval{2,1}.pres1,3), '       ', num2str(varianceval{2,1}.dres1,3), '       ', num2str(varianceval{2,1}.compslack1,3), '           ', num2str(varianceval{2,1}.time1,3), '           ', num2str(varianceval{2,1}.gradnum1,3), '        '...
        num2str(varianceval{2,1}.pres3,3), '        ', num2str(varianceval{2,1}.dres3,3), '          ', num2str(varianceval{2,1}.compslack3,3), '            ', num2str(varianceval{2,1}.time3,3), '           ', num2str(varianceval{2,1}.gradnum3,3), '         '...
        num2str(varianceval{2,1}.pres2,3), '        ', num2str(varianceval{2,1}.dres2,3), '       ', num2str(varianceval{2,1}.compslack2,3), '       ', num2str(varianceval{2,1}.time2,3), '           ', num2str(varianceval{2,1}.gradnum2,3), '         ',...
        num2str(varianceval{2,1}.pres4,3), '        ', num2str(varianceval{2,1}.dres4,3), '       ', num2str(varianceval{2,1}.compslack4,3), '             ', num2str(varianceval{2,1}.time4,3), '    ', num2str(varianceval{2,1}.gradnum4,3), '         ']);

fprintf(fileID, '%s\n', ['mean for rho = 10     ', num2str(meanvalues{3,1}.pres1,3), '       ', num2str(meanvalues{3,1}.dres1,3), '       ', num2str(meanvalues{3,1}.compslack1,3), '           ', num2str(meanvalues{3,1}.time1,3), '           ', num2str(meanvalues{3,1}.gradnum1,3), '        '...
        num2str(meanvalues{3,1}.pres3,3), '        ', num2str(meanvalues{3,1}.dres3,3), '          ', num2str(meanvalues{3,1}.compslack3,3), '            ', num2str(meanvalues{3,1}.time3,3), '           ', num2str(meanvalues{3,1}.gradnum3,3), '         '...
        num2str(meanvalues{3,1}.pres2,3), '        ', num2str(meanvalues{3,1}.dres2,3), '       ', num2str(meanvalues{3,1}.compslack2,3), '       ', num2str(meanvalues{3,1}.time2,3), '         ', num2str(meanvalues{3,1}.gradnum2,3), '         ',...
        num2str(meanvalues{3,1}.pres4,3), '        ', num2str(meanvalues{3,1}.dres4,3), '       ', num2str(meanvalues{3,1}.compslack4,3), '              ', num2str(meanvalues{3,1}.time4,3), '         ', num2str(meanvalues{3,1}.gradnum4,3), '         ']);
fprintf(fileID, '%s\n', ['mean for rho = 10     ', num2str(varianceval{3,1}.pres1,3), '       ', num2str(varianceval{3,1}.dres1,3), '       ', num2str(varianceval{3,1}.compslack1,3), '           ', num2str(varianceval{3,1}.time1,3), '      ', num2str(varianceval{3,1}.gradnum1,3), '        '...
        num2str(varianceval{3,1}.pres3,3), '        ', num2str(varianceval{3,1}.dres3,3), '          ', num2str(varianceval{3,1}.compslack3,3), '            ', num2str(varianceval{3,1}.time3,3), '            ', num2str(varianceval{3,1}.gradnum3,3), '         '...
        num2str(varianceval{3,1}.pres2,3), '        ', num2str(varianceval{3,1}.dres2,3), '       ', num2str(varianceval{3,1}.compslack2,3), '       ', num2str(varianceval{3,1}.time2,3), '          ', num2str(varianceval{3,1}.gradnum2,3), '         ',...
        num2str(varianceval{3,1}.pres4,3), '        ', num2str(varianceval{3,1}.dres4,3), '          ', num2str(varianceval{3,1}.compslack4,3), '             ', num2str(varianceval{3,1}.time4,3), '     ', num2str(varianceval{3,1}.gradnum4,3), '         ']);






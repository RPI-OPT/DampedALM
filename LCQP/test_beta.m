% % min 1/2*x'Qx+c'x s.t. Ax=b; 0 <= x_i <= x_u. % added box for nonconvex!
% Data: c; A,b ~ U[0,1]; Q ~ random PSD matrix with mean = I, covariance = I.
% nonconvex Q case!
clear; close all;
rand('seed',20220927); randn('seed',20220927); % same seed for debug!

m         = 10;
n         = 1000;
num_trial = 1;
rhonum    = 3;
tol       = 1e-3;
out6      = cell(rhonum,num_trial); 
out7      = cell(rhonum,num_trial);
out8      = cell(rhonum,num_trial);

x_l = -5 * ones(n,1); % lower box. -5. 0 for aipp to work.
x_u = 5 * ones(n,1); % upper box. 5
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

    %% Gaussian random case

    for num = 1:num_trial

        A     = [randn(m, n-m), eye(m)];
        b     = randn(m,1) + 0.1; % check if b in the box constraint
        c     = randn(n,1);
        [U,~] = qr(randn(n));
        v     = max(0, randn(n,1)*5); % v ~ Gaussian.
        Q     = U*diag(v)*U' + temp*eye(n);
        Q     = (Q + Q') / 2; % symmetrization is important
        rho   = min(eig(Q));
        if rho > 0
            error('rho > 0');
        end
        %% obtain a feasible initial point

        xinit = quadprog(eye(n),randn(n,1), [], [], A, b, x_l, x_u);
        %% Compare to iPALM
        opts1     = [];
        opts1.Q   = Q;
        opts1.r   = c;
        opts1.A   = A;
        opts1.b   = b;
        opts1.ell = -5;
        opts1.u   = 5;
        x0        = xinit;
        y0        = zeros(m,1);
        I         = eye(n);

        if hardness == 3
            beta0 = 5;
        elseif hardness == 2
            beta0 = 0.001;
        elseif hardness == 1
            beta0 = 0.0001;
        end
%         beta0            = -10*temp;
        rho      = -1*min(eig(Q));
        epsilon  = tol;
        maxouter = 1e4;
        mu       = min(eig(Q)) + 2*rho;

        t0 = tic;
        [x,y,out6{hardness,num}] = iPALM_QP(opts1,x0, y0,I, beta0,rho,epsilon,maxouter,mu);
        time2 = toc(t0);

        out6{hardness,num}.time = time2;
        %% Compare to NLNL-IAPIAL
        
        mf     = -1*min(eig(opts1.Q));
        sigma  = 0.5/sqrt(2);
        lambda = 0.5/(-1*temp);
        if hardness == 3
            tol1 = 1e3;
        else
            tol1 = 1e2;
        end
        t0    = tic;
        [x,out7{hardness,num}] = NLIPIAL_QP(opts1, ...
            x0,y0,I,mf, sigma, beta0, epsilon,maxouter,lambda,tol1);
        time3 = toc(t0);
        out7{hardness,num}.time = time3;
        
        tol1 = 1e-5;
        beta0 = -10*temp;
        t0    = tic;
        [x,out8{hardness,num}] = NLIPIAL_QP(opts1, ...
            x0,y0,I,mf, sigma, beta0, epsilon,maxouter,lambda, tol1);
        time4 = toc(t0);
        out8{hardness,num}.time = time4;
    end
end

%%

save result_beta out6 out7 out8

%%

colors = [0 0 0; 255 0 0; 0 0 255; 0 255 0; 255 0 255; 192 192 192; 128 128 105; 244 164 96; 8 46 84; 
    210 105 30; 61 89 171; 0 210 87; 0 255 127; 65 105 225;221 160 221; 227 23 13]/255;
lines   = {'-.' '--' ':' '-'};

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
semilogy(out6{1,1}.betaval, lines{1},'Color',colors(1,:), 'LineWidth',2);
hold on
semilogy(out7{1,1}.betaval(1:3000), lines{3},'Marker','diamond','Color','r', 'LineWidth',2)
hold on
semilogy(out8{1,1}.betaval(1:3000), lines{4},'Marker','*','Color','b', 'LineWidth',2)
hold off
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Outer iteration number", 'FontSize',16)
ylabel('Beta','FontSize',16)
legend('DPALM \beta_0 = 1e-4','NL-IAPIAL \beta_0 = 1e-4', 'NL-IAPIAL \beta_0 = 1', 'box','off','location','best', 'Orientation','vertical');
print(fig,'figures/beta_comparison1.pdf','-dpdf')

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
semilogy(out6{2,1}.betaval, lines{1},'Color',colors(1,:), 'LineWidth',2);
hold on
semilogy(out7{2,1}.betaval(1:3000), lines{4},'Marker','diamond','Color','r', 'LineWidth',2)
hold on
semilogy(out8{2,1}.betaval(1:3000), lines{4},'Marker','*','Color','b', 'LineWidth',2)
hold off
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Outer iteration number", 'FontSize',16)
ylabel('Beta','FontSize',16)
legend('DPALM \beta_0 = 1e-3','NL-IAPIAL \beta_0 = 1e-3', 'NL-IAPIAL \beta_0 = 10', 'box','off','location','best', 'Orientation','vertical');
print(fig,'figures/beta_comparison2.pdf','-dpdf')

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
semilogy(out6{3,1}.betaval, lines{1},'Color',colors(1,:), 'LineWidth',2);
hold on
semilogy(out7{3,1}.betaval(1:3000), lines{4},'Marker','diamond', 'Color', 'r', 'LineWidth',2)
hold on
semilogy(out8{3,1}.betaval(1:3000), lines{4},'Marker','*','Color','b', 'LineWidth',2)
hold off
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Outer iteration number", 'FontSize',16)
ylabel('Beta','FontSize',16)
legend('DPALM \beta_0 = 5','NL-IAPIAL \beta_0 = 5', 'NL-IAPIAL \beta_0 = 100' ,'box','off','location','best', 'Orientation','vertical');
print(fig,'figures/beta_comparison3.pdf','-dpdf')



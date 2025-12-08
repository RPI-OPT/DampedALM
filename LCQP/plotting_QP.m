colors = [0 0 0; 255 0 0; 0 0 255; 0 255 0; 255 0 255; 192 192 192; 128 128 105; 244 164 96; 8 46 84; 
    210 105 30; 61 89 171; 0 210 87; 0 255 127; 65 105 225;221 160 221; 227 23 13]/255;
lines   = {'-.' '--' ':' '-'};

load result_compare.mat

i     = 1;
delta = 25;

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
[a, id] = max(out1{1,i}.gradnum);
semilogy(out1{1,i}.gradnum([1:2*delta:id,id]),out1{1,i}.primalviol([1:2*delta:id,id]),'Marker','x','Color',colors(5,:), 'LineWidth',2);hold on
[a, id] = max(out2{1,i}.gradnum);
semilogy(out2{1,i}.gradnum([1:2*delta:id,id]),out2{1,i}.primalviol([1:2*delta:id,id]),lines{4},'Marker','diamond','Color','r', 'LineWidth',2); hold on
[a, id] = max(out3{1,i}.gradnum);
semilogy(out3{1,i}.gradnum(1:1*delta:id),out3{1,i}.primalviol(1:delta*1:id),lines{2},'Marker','square','Color',colors(3,:), 'LineWidth',2); hold on
[a, id] = max(out4{1,i}.gradnum);
semilogy(out4{1,i}.gradnum(1:delta:id),out4{1,i}.primalviol(1:delta:id),lines{1},'Color',colors(1,:), 'LineWidth',2); hold on
[a, id] = max(out5{1,i}.gradnum);
semilogy(out5{1,i}.gradnum(11:delta*10:id),out5{1,i}.primalviol(11:delta*10:id),lines{3},'Marker','+','Color',colors(7,:), 'LineWidth',2); hold on
hold off
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to PF','FontSize',16)
print(fig,'figures/primal_feasibility_qp1.pdf','-dpdf')

delta = 25;
fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
[a, id] = max(out1{1,i}.gradnum);
semilogy(out1{1,i}.gradnum([1:delta:id,id]),out1{1,i}.dualviol([1:delta:id,id]),'Marker','x','Color',colors(5,:), 'LineWidth',2);hold on
[a, id] = max(out2{1,i}.gradnum);
semilogy(out2{1,i}.gradnum([1:2*delta:id,id]),out2{1,i}.dualviol([1:2*delta:id,id]),lines{4},'Marker','diamond','Color','r', 'LineWidth',2); hold on
[a, id] = max(out3{1,i}.gradnum);
semilogy(out3{1,i}.gradnum(1:1*delta:id),out3{1,i}.dualviol(1:delta*1:id),lines{2},'Marker','square','Color',colors(3,:), 'LineWidth',2); hold on
[a, id] = max(out4{1,i}.gradnum);
semilogy(out4{1,i}.gradnum(1:delta:id),out4{1,i}.dualviol(1:delta:id),lines{1},'Color',colors(1,:), 'LineWidth',2); hold on
[a, id] = max(out5{1,i}.gradnum);
semilogy(out5{1,i}.gradnum(11:delta*10:id),out5{1,i}.dualviol(11:delta*10:id),lines{3},'Marker','+','Color',colors(7,:), 'LineWidth',2); hold on
hold off
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to DF','FontSize',16)
print(fig,'figures/dual_feasibility_qp1.pdf','-dpdf')

delta = 100;
fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
[a, id] = max(out1{2,i}.gradnum);
semilogy(out1{2,i}.gradnum([1:delta:id,id]),out1{2,i}.primalviol([1:delta:id,id]),'Marker','x','Color',colors(5,:), 'LineWidth',2);hold on
[a, id] = max(out2{2,i}.gradnum);
semilogy(out2{2,i}.gradnum([1:2*delta:id,id]),out2{2,i}.primalviol([1:2*delta:id,id]),lines{4},'Marker','diamond','Color','r', 'LineWidth',2); hold on
[a, id] = max(out3{2,i}.gradnum);
semilogy(out3{2,i}.gradnum(1:1*delta:id),out3{2,i}.primalviol(1:delta*1:id),lines{2},'Marker','square','Color',colors(3,:), 'LineWidth',2); hold on
[a, id] = max(out4{2,i}.gradnum);
semilogy(out4{2,i}.gradnum(1:delta:id),out4{2,i}.primalviol(1:delta:id),lines{1},'Color',colors(1,:), 'LineWidth',2); hold on
[a, id] = max(out5{2,i}.gradnum);
semilogy(out5{2,i}.gradnum(11:delta*10:id),out5{2,i}.primalviol(11:delta*10:id),lines{3},'Marker','+','Color',colors(7,:), 'LineWidth',2); hold on
hold off
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to PF','FontSize',16)
print(fig,'figures/primal_feasibility_qp2.pdf','-dpdf')

delta = 50;
fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
[a, id] = max(out1{2,i}.gradnum);
semilogy(out1{2,i}.gradnum([1:delta:id,id]),out1{2,i}.dualviol([1:delta:id,id]),'Marker','x','Color',colors(5,:), 'LineWidth',2);hold on
[a, id] = max(out2{2,i}.gradnum);
semilogy(out2{2,i}.gradnum([1:2*delta:id,id]),out2{2,i}.dualviol([1:2*delta:id,id]),lines{4},'Marker','diamond','Color','r', 'LineWidth',2); hold on
[a, id] = max(out3{2,i}.gradnum);
semilogy(out3{2,i}.gradnum(1:1*delta:id),out3{2,i}.dualviol(1:delta*1:id),lines{2},'Marker','square','Color',colors(3,:), 'LineWidth',2); hold on
[a, id] = max(out4{2,i}.gradnum);
semilogy(out4{2,i}.gradnum(1:delta:id),out4{2,i}.dualviol(1:delta:id),lines{1},'Color',colors(1,:), 'LineWidth',2); hold on
[a, id] = max(out5{2,i}.gradnum);
semilogy(out5{2,i}.gradnum(11:delta*10:id),out5{2,i}.dualviol(11:delta*10:id),lines{3},'Marker','+','Color',colors(7,:), 'LineWidth',2); hold on
hold off
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to DF','FontSize',16)
print(fig,'figures/dual_feasibility_qp2.pdf','-dpdf')

delta = 5;
fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
[a, id] = max(out1{3,i}.gradnum);
semilogy(out1{3,i}.gradnum([1:delta:id,id]),out1{3,i}.primalviol([1:delta:id,id]),'Marker','x','Color',colors(5,:), 'LineWidth',2);hold on
[a, id] = max(out2{3,i}.gradnum);
semilogy(out2{3,i}.gradnum([1:1*delta:id,id]),out2{3,i}.primalviol([1:1*delta:id,id]),lines{4},'Marker','diamond','Color','r', 'LineWidth',2); hold on
[a, id] = max(out3{3,i}.gradnum);
semilogy(out3{3,i}.gradnum(1:1*delta:id),out3{3,i}.primalviol(1:delta*1:id),lines{2},'Marker','square','Color',colors(3,:), 'LineWidth',2); hold on
[a, id] = max(out4{3,i}.gradnum);
semilogy(out4{3,i}.gradnum(1:delta:id),out4{3,i}.primalviol(1:delta:id),lines{1},'Color',colors(1,:), 'LineWidth',2); hold on
[a, id] = max(out5{3,i}.gradnum);
semilogy(out5{3,i}.gradnum(11:delta*10:id),out5{3,i}.primalviol(11:delta*10:id),lines{3},'Marker','+','Color',colors(7,:), 'LineWidth',2); hold on
hold off
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to PF','FontSize',16)
print(fig,'figures/primal_feasibility_qp3.pdf','-dpdf')

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
[a, id] = max(out1{3,i}.gradnum);
semilogy(out1{3,i}.gradnum([1:delta:id,id]),out1{3,i}.dualviol([1:delta:id,id]),'Marker','x','Color',colors(5,:), 'LineWidth',2);hold on
[a, id] = max(out2{3,i}.gradnum);
semilogy(out2{3,i}.gradnum([1:1*delta:id,id]),out2{3,i}.dualviol([1:1*delta:id,id]),lines{4},'Marker','diamond','Color','r', 'LineWidth',2); hold on
[a, id] = max(out3{3,i}.gradnum);
semilogy(out3{3,i}.gradnum(1:1*delta:id),out3{3,i}.dualviol(1:delta*1:id),lines{2},'Marker','square','Color',colors(3,:), 'LineWidth',2); hold on
[a, id] = max(out4{3,i}.gradnum);
semilogy(out4{3,i}.gradnum(1:delta:id),out4{3,i}.dualviol(1:delta:id),lines{1},'Color',colors(1,:), 'LineWidth',2); hold on
[a, id] = max(out5{3,i}.gradnum);
semilogy(out5{3,i}.gradnum(11:delta*10:id),out5{3,i}.dualviol(11:delta*10:id),lines{3},'Marker','+','Color',colors(7,:), 'LineWidth',2); hold on
hold off
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to DF','FontSize',16)
print(fig,'figures/dual_feasibility_qp3.pdf','-dpdf')
%%
%legend
fig = figure('papersize',[10,1],'paperposition',[0,0,9,1.5]);
[a, id] = max(out1{3,i}.gradnum);
semilogy(out1{3,i}.gradnum([1:delta:id,id]),out1{3,i}.dualviol([1:delta:id,id]),'Marker','x','Color',colors(5,:), 'LineWidth',2);hold on
[a, id] = max(out2{3,i}.gradnum);
semilogy(out2{3,i}.gradnum([1:1*delta:id,id]),out2{3,i}.dualviol([1:1*delta:id,id]),lines{4},'Marker','diamond','Color','r', 'LineWidth',2); hold on
[a, id] = max(out3{3,i}.gradnum);
semilogy(out3{3,i}.gradnum(1:1*delta:id),out3{3,i}.dualviol(1:delta*1:id),lines{2},'Marker','square','Color',colors(3,:), 'LineWidth',2); hold on
[a, id] = max(out4{3,i}.gradnum);
semilogy(out4{3,i}.gradnum(1:delta:id),out4{3,i}.dualviol(1:delta:id),lines{1},'Color',colors(1,:), 'LineWidth',2); hold on
[a, id] = max(out5{3,i}.gradnum);
semilogy(out5{3,i}.gradnum(11:delta*10:id),out5{3,i}.dualviol(11:delta*10:id),lines{3},'Marker','+','Color',colors(7,:), 'LineWidth',2); hold on
hold off
legend('DPALM','NL-IAPIAL', 'IPL(A)', 'HiAPeM', 'LiMEAL' ,'box','off','location','southoutside', 'Orientation','horizontal');
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Dual infeasibility','FontSize',16)

print(fig,'figures/legend.pdf','-dpdf')
%%
% semilogy(idxsgdx,sgdobj(idxsgd), lines{1},'Color',colors(1,:),'linewidth',2);hold on
% % semilogy(sgsfullobj, lines{2},'Color',colors(9,:),'linewidth',2);hold on
% semilogy(idxsgdx,adamobj(idxsgd), lines{3},'Color',colors(3,:),'linewidth',2);hold on
% %semilogy(idxsgdx,delta2obj(idxsgd),lines{2},'Marker','h','Color',colors(4,:),'linewidth',2);
% %semilogy(idxsgdx,adagraobj(idxsgd), lines{3},'Marker','o','Color',colors(5,:),'linewidth',2);hold on
% %semilogy(idxsgdx,adagra2obj(idxsgd), lines{2},'Marker','x','Color',colors(6,:),'linewidth',2);hold on
% semilogy(idxsgdx,deltaobj(idxsgd), lines{3},'Marker','+','Color',colors(7,:),'linewidth',2);hold on
% %semilogy(idxsgdx,proxsgdobj(idxsgd), lines{2},'Marker','p','Color',colors(8,:),'linewidth',2);hold on
% semilogy(idxour,loss3(idxour),lines{4},'Marker','diamond','Color','r','linewidth',2);
%  option = ['-' markers(2)];
%  semilogy((1:size(lossadmm2,2))*ttt1/size(lossadmm2,2)+Adadeltatt/20, lossadmm2, lines{1},'Color','b');
%plot((1:70)*aa/70, log10(lossadmm(1:70)));
%% Plot results -- Convergence plot -- Purely stochastic gradient
% hold off
% grid on
% xlabel('Iteration\timesbatchsize/N');xlim([0,min(maxsize,200)]);
% ylabel('TrainErr');ylim([min(loss3),max(sgdobj)*50]);
% solvers = { ...
%     'Vanilla SGD (batch)','Adam',...%'Adammax', 'AdaGrad', 'AdaGradDecay', 
%     'Adadelta', ...%'ProxSGD',
%     'IALAM'
%     };
% legend(solvers,'box','off','location','best');
% set(gca,'linewidth',1.6);
% set(gca,'fontsize',16,'fontweight','bold');
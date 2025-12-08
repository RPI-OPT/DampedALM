%%
colors = [0 0 0; 255 0 0; 0 0 255; 0 255 0; 255 0 255; 192 192 192; 128 128 105; 244 164 96; 8 46 84; 
    210 105 30; 61 89 171; 0 210 87; 0 255 127; 65 105 225;221 160 221; 227 23 13]/255;
lines   = {'-.' '--' ':' '-'};

%plotting the first instance of the result
fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
semilogy(out1{1,1}.gradnum,out1{1,1}.primalviol, lines{1},'Color',colors(1,:), 'LineWidth',2)
hold on
semilogy(out2{1,1}.gradnum,out2{1,1}.primalviol,'Marker','x','Color',colors(5,:), 'LineWidth',2)
hold on 
semilogy(out3{1,1}.gradnum,out3{1,1}.primalviol,lines{4},'Marker','diamond','Color','r', 'LineWidth',2)
hold on 
semilogy(out4{1,1}.gradnum,out4{1,1}.primalviol,lines{2},'Marker','diamond','Color',colors(3,:), 'LineWidth',2)
hold off
%legend('HiAPeM', 'DPALM', 'NL-IAPIAL','LiMEAL' ,'box','off','location','best');
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to PF','FontSize',16)
% saveas(gcf,'primal_feasibility_qcqp1.png')
print(fig,'figures/newprimal_feasibility_qcqp1.pdf','-dpdf')

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
semilogy(out1{1,1}.gradnum,out1{1,1}.dualviol, lines{1},'Color',colors(1,:), 'LineWidth',2)
hold on
semilogy(out2{1,1}.gradnum,out2{1,1}.dualviol,'Marker','x','Color',colors(5,:), 'LineWidth',2)
hold on
semilogy(out3{1,1}.gradnum,out3{1,1}.dualviol,'Marker','diamond','Color', 'r', 'LineWidth',2)
hold on 
semilogy(out4{1,1}.gradnum,out4{1,1}.dualviol,lines{2},'Marker','diamond','Color',colors(3,:), 'LineWidth',2)
hold off
%legend('HiAPeM', 'DPALM', 'NL-IAPIAL','LiMEAL' ,'box','off','location','best');
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to DF','FontSize',16)
% saveas(gcf,'dual_feasibility_qcqp1.png')
print(fig,'figures/newdual_feasibility_qcqp1.pdf','-dpdf')

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
semilogy(out1{1,1}.gradnum,out1{1,1}.compslack, lines{1},'Color',colors(1,:), 'LineWidth',2)
hold on
semilogy(out2{1,1}.gradnum,out2{1,1}.compslack,'Marker','x','Color',colors(5,:), 'LineWidth',2)
hold on
semilogy(out3{1,1}.gradnum,out3{1,1}.compslack,'Marker','diamond','Color', 'r', 'LineWidth',2)
hold on 
semilogy(out4{1,1}.gradnum,out4{1,1}.compslack,lines{2},'Marker','diamond','Color',colors(3,:), 'LineWidth',2)
hold off
%legend('HiAPeM', 'DPALM', 'NL-IAPIAL','LiMEAL' ,'box','off','location','best');
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to CS','FontSize',16)
% saveas(gcf,'CompSlack_qcqp1.png')
print(fig,'figures/newCompSlack_qcqp1.pdf','-dpdf')

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
semilogy(out1{2,1}.gradnum,out1{2,1}.primalviol, lines{1},'Color',colors(1,:), 'LineWidth',2)
hold on
semilogy(out2{2,1}.gradnum,out2{2,1}.primalviol,'Marker','x','Color',colors(5,:), 'LineWidth',2)
hold on 
semilogy(out3{2,1}.gradnum,out3{2,1}.primalviol,lines{4},'Marker','diamond','Color','r', 'LineWidth',2)
hold on 
semilogy(out4{2,1}.gradnum,out4{2,1}.primalviol,lines{2},'Marker','diamond','Color',colors(3,:), 'LineWidth',2)
hold off
%legend('HiAPeM', 'DPALM', 'NL-IAPIAL','LiMEAL' ,'box','off','location','best');
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to PF','FontSize',16)
% saveas(gcf,'primal_feasibility_qcqp1.png')
print(fig,'figures/newprimal_feasibility_qcqp2.pdf','-dpdf')

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
semilogy(out1{2,1}.gradnum,out1{2,1}.dualviol, lines{1},'Color',colors(1,:), 'LineWidth',2)
hold on
semilogy(out2{2,1}.gradnum,out2{2,1}.dualviol,'Marker','x','Color',colors(5,:), 'LineWidth',2)
hold on
semilogy(out3{2,1}.gradnum,out3{2,1}.dualviol,'Marker','diamond','Color', 'r', 'LineWidth',2)
hold on 
semilogy(out4{2,1}.gradnum,out4{2,1}.dualviol,lines{2},'Marker','diamond','Color',colors(3,:), 'LineWidth',2)
hold off
%legend('HiAPeM', 'DPALM', 'NL-IAPIAL','LiMEAL' ,'box','off','location','best');
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to DF','FontSize',16)
% saveas(gcf,'dual_feasibility_qcqp1.png')
print(fig,'figures/newdual_feasibility_qcqp2.pdf','-dpdf')

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
semilogy(out1{2,1}.gradnum,out1{2,1}.compslack, lines{1},'Color',colors(1,:), 'LineWidth',2)
hold on
semilogy(out2{2,1}.gradnum,out2{2,1}.compslack,'Marker','x','Color',colors(5,:), 'LineWidth',2)
hold on
semilogy(out3{2,1}.gradnum,out3{2,1}.compslack,'Marker','diamond','Color', 'r', 'LineWidth',2)
hold on 
semilogy(out4{2,1}.gradnum,out4{2,1}.compslack,lines{2},'Marker','diamond','Color',colors(3,:), 'LineWidth',2)
hold off
%legend('HiAPeM', 'DPALM', 'NL-IAPIAL','LiMEAL' ,'box','off','location','best');
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to CS','FontSize',16)
% saveas(gcf,'CompSlack_qcqp1.png')
print(fig,'figures/newCompSlack_qcqp2.pdf','-dpdf')

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
loglog(out1{3,1}.gradnum,out1{3,1}.primalviol, lines{1},'Color',colors(1,:), 'LineWidth',2)
hold on
loglog(out2{3,1}.gradnum,out2{3,1}.primalviol,'Marker','x','Color',colors(5,:), 'LineWidth',2)
hold on 
loglog(out3{3,1}.gradnum,out3{3,1}.primalviol,lines{4},'Marker','diamond','Color','r', 'LineWidth',2)
hold on 
semilogy(out4{3,1}.gradnum,out4{3,1}.primalviol,lines{2},'Marker','diamond','Color',colors(3,:), 'LineWidth',2)
hold off
%legend('HiAPeM', 'DPALM', 'NL-IAPIAL','LiMEAL' ,'box','off','location','best');
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to PF','FontSize',16)
% saveas(gcf,'primal_feasibility_qcqp1.png')
print(fig,'figures/newprimal_feasibility_qcqp3.pdf','-dpdf')

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
loglog(out1{3,1}.gradnum,out1{3,1}.dualviol, lines{1},'Color',colors(1,:), 'LineWidth',2)
hold on
loglog(out2{3,1}.gradnum,out2{3,1}.dualviol,'Marker','x','Color',colors(5,:), 'LineWidth',2)
hold on
loglog(out3{3,1}.gradnum,out3{3,1}.dualviol,'Marker','diamond','Color', 'r', 'LineWidth',2)
hold on 
semilogy(out4{3,1}.gradnum,out4{3,1}.dualviol,lines{2},'Marker','diamond','Color',colors(3,:), 'LineWidth',2)
hold off
%legend('HiAPeM', 'DPALM', 'NL-IAPIAL','LiMEAL' ,'box','off','location','best');
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to DF','FontSize',16)
% saveas(gcf,'dual_feasibility_qcqp1.png')
print(fig,'figures/newdual_feasibility_qcqp3.pdf','-dpdf')

fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
loglog(out1{3,1}.gradnum,out1{3,1}.compslack, lines{1},'Color',colors(1,:), 'LineWidth',2)
hold on
loglog(out2{3,1}.gradnum,out2{3,1}.compslack,'Marker','x','Color',colors(5,:), 'LineWidth',2)
hold on
loglog(out3{3,1}.gradnum,out3{3,1}.compslack,'Marker','diamond','Color', 'r', 'LineWidth',2)
hold on 
semilogy(out4{1,1}.gradnum,out4{1,1}.primalviol,lines{2},'Marker','diamond','Color',colors(3,:), 'LineWidth',2)
hold off
%legend('HiAPeM', 'DPALM', 'NL-IAPIAL','LiMEAL' ,'box','off','location','best');
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Violation to CS','FontSize',16)
% saveas(gcf,'CompSlack_qcqp1.png')
print(fig,'figures/newCompSlack_qcqp3.pdf','-dpdf')
%%
fig = figure('papersize',[9,1],'paperposition',[0,0,8,1.5]);
loglog(out1{3,1}.gradnum,out1{3,1}.compslack, lines{1},'Color',colors(1,:), 'LineWidth',2)
hold on
loglog(out2{3,1}.gradnum,out2{3,1}.compslack,'Marker','x','Color',colors(5,:), 'LineWidth',2)
hold on
loglog(out3{3,1}.gradnum,out3{3,1}.compslack,'Marker','diamond','Color', 'r', 'LineWidth',2)
hold on 
semilogy(out4{3,1}.gradnum,out4{3,1}.compslack,lines{2},'Marker','diamond','Color',colors(3,:), 'LineWidth',2)
hold off
%legend('HiAPeM', 'DPALM', 'NL-IAPIAL','LiMEAL' ,'box','off','location','best');
set(gca,'linewidth',1.6);
set(gca,'fontsize',16,'fontweight','bold');
xlabel("Number of gradient","FontSize",16)
ylabel('Complementary Slackness','FontSize',16)
legend('HiAPeM', 'DPALM', 'NL-IAPIAL', 'IPL(A)', 'box','off','location','southoutside','Orientation','horizontal');
print(fig,'figures/legend.pdf','-dpdf')


load('treated_untreated_codon_cluster_022724.mat')
savefig_path = '../subfigures/';
a = 1;

%MG132 DD6 TNF-H  1578
%MG132 DD6 Pam-H  1590 

% CpG: [137 180 66]/255
% Pam: [229 129 56]/255
% TNF: [119 180 202]/255
figure(1)
hold off
plot(1:6,metric_cluster_norm(1578,:),'Color',[119 180 202]/255); hold on% TNF
plot(1:6,metric_cluster_norm(1590,:),'Color',[229 129 56]/255); hold on% PAm

figure(1)
plot(1:6,untreated_metric_norm(3,:),'--','Color',[119 180 202]/255); hold on
plot(1:6,untreated_metric_norm(15,:),'--','Color',[229 129 56]/255); hold on

ylim([-2,4])
xticks(1:6)
xticklabels('')
yticklabels('')

saveas(gcf,strcat(savefig_path,'Figure5_suppMG132_TNF_Pam_codon'),'epsc')



%MG132 DD6 TNF-H  1578
%MG132 DD6 CpG-H  1584 

% CpG: [137 180 66]/255
% Pam: [229 129 56]/255
% TNF: [119 180 202]/255
figure(1)
hold off
plot(1:6,metric_cluster_norm(1590,:),'Color',[229 129 56]/255); hold on% TNF
plot(1:6,metric_cluster_norm(1584,:),'Color',[137 180 66]/255); hold on% PAm

figure(1)
plot(1:6,untreated_metric_norm(15,:),'--','Color',[229 129 56]/255); hold on
plot(1:6,untreated_metric_norm(9,:),'--','Color',[137 180 66]/255); hold on

ylim([-2,4])
xticks(1:6)
xticklabels('')
yticklabels('')

saveas(gcf,strcat(savefig_path,'Figure5_suppMG132_CpG_Pam_codon'),'epsc')


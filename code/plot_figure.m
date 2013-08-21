load('../data/prelim_fch.txt')
load('../data/prelim_ha.txt')
load('../data/prelim_stackSVD.txt')

figure
hold on; grid on;
plot(prelim_ha(:,1),prelim_ha(:,2),'-g.', 'LineWidth',1.5,'MarkerSize',12);
plot(prelim_fch(:,1),prelim_fch(:,2), '-b.', 'LineWidth',1.5,'MarkerSize',12);
plot(prelim_stackSVD(:,1),prelim_stackSVD(:,2),'-r.', 'LineWidth',1.5,'MarkerSize',12);
plot(prelim_stackSVD(:,1),0.1428*ones(9,1),'--k', 'LineWidth',1.5,'MarkerSize',12);
legend('HA','FCH','StackSVD','chance=1/7=0.1428','Location', 'NorthWest');
xlabel('number of subjects')
ylabel('accuracy')
%axis([2 10 0.1 0.5])
saveas(gcf,  '/Users/ChimatChen/Dropbox/Research/Hyperalignment/prelim_accuracy', 'epsc')

figure
hold on; grid on;
plot(prelim_ha(:,1),prelim_ha(:,3),'-g.', 'LineWidth',1.5,'MarkerSize',12);
plot(prelim_fch(:,1),prelim_fch(:,3), '-b.', 'LineWidth',1.5,'MarkerSize',12);
plot(prelim_stackSVD(:,1),prelim_stackSVD(:,3),'-r.', 'LineWidth',1.5,'MarkerSize',12);
legend('HA','FCH','StackSVD','Location', 'NorthWest');
xlabel('number of subjects')
ylabel('alignment time(sec)')
saveas(gcf,  '/Users/ChimatChen/Dropbox/Research/Hyperalignment/prelim_time', 'epsc')

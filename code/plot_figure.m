%load('../data/prelim_fch.txt')
%load('../data/prelim_ha.txt')
%load('../data/prelim_stackSVD.txt')
clear
close all
data_path = '../data/output/'

exp_name = { 'HA' 'FCH' 'JointSVD'  };
nTR_set = [100 200 400 800];
data_accuracy = zeros(2,9,numel(nTR_set),numel(exp_name));
data_time = zeros(9,numel(nTR_set),numel(exp_name));
load('../data/output/data_accuracy.mat')
load('../data/output/data_time.mat')
%{
for k =1:numel(exp_name)
  for i=1:9
    for j=1:numel(nTR_set)
      fprintf('%d %d %d\n',i,j,k);
      load([data_path exp_name{k} '_PRELIM_accuracy_' num2str(i+1) '_' num2str(nTR_set(j)) '.mat']);
      load([data_path exp_name{k} '_PRELIM_time_' num2str(i+1) '_' num2str(nTR_set(j)) '.mat']);
      data_accuracy(1,i,j,k) = mean(accuracy);
      data_accuracy(2,i,j,k) = std(accuracy);
      data_time(i,j,k) = time_spent;
    end
  end
end
%}


for j=1:numel(nTR_set)
  figure
  hold on; grid on;
  errorbar(2:10,data_accuracy(1,:,j,1),data_accuracy(2,:,j,1),'-g.', 'LineWidth',1.5,'MarkerSize',12);
  errorbar(2:10,data_accuracy(1,:,j,2),data_accuracy(2,:,j,2), '-b.', 'LineWidth',1.5,'MarkerSize',12);
  errorbar(2:10,data_accuracy(1,:,j,3),data_accuracy(2,:,j,3),'-r.', 'LineWidth',1.5,'MarkerSize',12);
  plot(2:10,0.1428*ones(9,1),'--k', 'LineWidth',1.5,'MarkerSize',12);
  legend([exp_name 'chance=1/7=0.1428'],'Location', 'NorthWest');
  xlabel('number of subjects')
  ylabel('accuracy')
  title([num2str(nTR_set(j)) 'TRs' ' accuracy'])
  axis([2 10 0 0.6])
  saveas(gcf,  ['/Users/ChimatChen/Dropbox/Research/Hyperalignment/prelim_accuracy_' num2str(nTR_set(j)) 'TRs']  , 'epsc')


  figure
  hold on; grid on;
  plot(2:10,data_time(:,j,1),'-g.', 'LineWidth',1.5,'MarkerSize',12);
  plot(2:10,data_time(:,j,2), '-b.', 'LineWidth',1.5,'MarkerSize',12);
  plot(2:10,data_time(:,j,3),'-r.', 'LineWidth',1.5,'MarkerSize',12);
  legend(exp_name,'Location', 'NorthWest');
  xlabel('number of subjects')
  ylabel('alignment time(sec)')
  title([num2str(nTR_set(j)) 'TRs' ' time'])
  saveas(gcf, [ '/Users/ChimatChen/Dropbox/Research/Hyperalignment/prelim_time_' num2str(nTR_set(j)) 'TRs'], 'epsc')
end

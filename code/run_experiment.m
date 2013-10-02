clear

%%set parameters
fprintf('loading options\n')
options.exptype = 'StackSVD_PRELIM' % PRELIM: preliminary or LOO: leave one out
alignment_handle = @StackSVD;

options.nvoxel = 1000
options.nTR_all =[50:50:1000]


if strcmp(options.exptype,'RHA_PRELIM') | strcmp(options.exptype,'RHAchol_PRELIM') 
  options.alpha_all = [0.7]
else
  options.alpha_all = [1]
end

for alpha_idx = 1:numel(options.alpha_all) 
for a=1:numel(options.nTR_all)

options.alpha = options.alpha_all(alpha_idx)
options.nTR = options.nTR_all(a)

options.time = clock;
options.time = [date '-' num2str(options.time(4)) num2str(options.time(5))]
%options.input_path = '/mnt/cd/ramadge/pohsuan/hyperalignment/data/input/'; 
%options.working_path = '/mnt/cd/ramadge/pohsuan/hyperalignment/data/working/' ; 
%options.output_path = '/mnt/cd/ramadge/pohsuan/hyperalignment/data/output/' ; 
options.input_path = ['/mnt/cd/fastscratch/pohsuan/X_matrix_vs' num2str(options.nvoxel) '_' num2str(options.nTR) '/']; 
options.working_path = ['/mnt/cd/fastscratch/pohsuan/X_matrix_vs' num2str(options.nvoxel) '_' num2str(options.nTR) '/working/'] ; 
options.output_path = ['/mnt/cd/fastscratch/pohsuan/X_matrix_vs' num2str(options.nvoxel) '_' num2str(options.nTR) '/output/' ]; 
options.random_seed = 1;
options.training_data_all = {'t_cb_vt_rh','t_dm_vt_rh','t_hj_vt_rh','t_kd_vt_rh',...
                    't_kl_vt_rh','t_mh_vt_rh','t_ph_vt_rh','t_rb_vt_rh',...
                   't_se_vt_rh','t_sm_vt_rh'}

options.testing_data = 'mkdg_vt_tlrc' ;

mkdir(options.input_path)
mkdir(options.working_path)
mkdir(options.output_path)

libsvm_path = '/mnt/cd/ramadge/pohsuan/libsvm-3.17-rondo/matlab/';
origin_path = '/mnt/cd/ramadge/pohsuan/hyperalignment/code/';

%libsvmpath = '../../libsvm-3.17-mac/matlab'
%origin_path = '../../hyperalignment/code/';

% create diary to save log
%filename = [options.time '_diary.txt'];
%eval(['diary ',filename]);
%diary on;

for m=10

options.training_data = options.training_data_all(1:m)

rng(options.random_seed);
data = zeros(options.nTR, options.nvoxel, m);
%load alignment data "raider of the lost track" (TRAINING)
for i=1:size(options.training_data,2),
  fprintf('loading data: %_s\n',options.training_data{i});
  tmp=load([options.input_path options.training_data{i}]);
  %X(:,:,i)=tmp.(options.training_data{i});
  data(:,:,i) = tmp.X;
%  tmp_data = tmp.(options.training_data{i}); 
%  data(:,:,i) = tmp_data(options.nTR,options.nvoxel);
end
  data=data(1:options.nTR,:,:);

%load testing data "monkey-dog" (TRAINING/TESTING)
load([options.input_path options.testing_data]);

clear tmp;

%%hyperalignment (TRAINING)
fprintf('start hyperalignment\n');

if strcmp(options.exptype,'RHA_PRELIM') | strcmp(options.exptype,'RHAchol_PRELIM') 
  RG_filename = [options.exptype '_RG_data_' num2str(size(options.training_data,2)) '_' num2str(options.nTR) '_' num2str(options.alpha*100)];
  time_filename =  [options.exptype '_time_' num2str(size(options.training_data,2)) '_' num2str(options.nTR) '_' num2str(options.alpha*100)];
else
  RG_filename = [options.exptype '_RG_data_' num2str(size(options.training_data,2)) '_' num2str(options.nTR)];
  time_filename =  [options.exptype '_time_' num2str(size(options.training_data,2)) '_' num2str(options.nTR)];
end



if exist([options.output_path RG_filename '.mat'], 'file')
  fprintf('load RG data from file\n');
  load([options.output_path RG_filename '.mat']);
  load([[ options.output_path time_filename] '.mat']);
else
  fprintf('calculate RG\n');
  tic
  if strcmp(options.exptype,'RHA_PRELIM') | strcmp(options.exptype,'RHAchol_PRELIM')
    [R G]=alignment_handle(data,options.alpha);
  else
    [R G]=alignment_handle(data);
  end
  time_spent = toc
  save( [ options.output_path RG_filename ], 'R', 'G');
end

fprintf('Use R to rotate data\n');

XRotate = zeros(size(XRtst,1),size(R(:,:,1),2));

for subj=1:size(options.training_data,2),
  idx = find(subjid==subj);
  for i=1:numel(idx)
    %XRtst(idx(i),:) = XRtst(idx(i),:)*R(:,:,subjid(idx(i)));
    XRotate(idx(i),:) = XRtst(idx(i),:)*R(:,:,subjid(idx(i)));
  end
end

fprintf('DO SVM\n');
cd(libsvm_path)

accuracy =[]; allCount =0;
for subj=1:size(options.training_data,2),
  fprintf('SVM subject:%d\n',subj);
  idx_trn=find(subjid~=subj);
  idx_tst=find(subjid==subj);

  svm_model=svmtrain2(lblid(idx_trn),XRotate(idx_trn,:),...
                      '-s 1 -t 0 -n 0.5 -p 0.001 -b 1 -h 0');
  [p1 p2 p3] = svmpredict(lblid(idx_tst),XRotate(idx_tst,:),svm_model,'-b 1');
  
  accuracy = [ accuracy sum(p1 == lblid(idx_tst))/length(p1)] ;
end

accuracy 
time_spent
if strcmp(options.exptype,'RHA_PRELIM') | strcmp(options.exptype,'RHAchol_PRELIM')
  accuracy_filename =  [options.exptype '_accuracy_' num2str(size(options.training_data,2)) '_' num2str(options.nTR) '_' num2str(options.alpha*100)];
else
  accuracy_filename =  [options.exptype '_accuracy_' num2str(size(options.training_data,2)) '_' num2str(options.nTR)];
end

save( [ options.output_path accuracy_filename], 'accuracy');
save( [ options.output_path time_filename], 'time_spent');

cd(origin_path)
%diary off;
end

end
end

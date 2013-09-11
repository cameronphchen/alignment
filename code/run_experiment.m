clear

%%set parameters
fprintf('loading options\n')
options.exptype = 'HA_PRELIM' % PRELIM: preliminary or LOO: leave one out
alignment_handle = @HA;

options.nvoxel = 400
options.nTR =800
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

for m=2:10

options.training_data = options.training_data_all(1:m)

rng(options.random_seed);

%load alignment data "raider of the lost track" (TRAINING)
for i=1:size(options.training_data,2),
  fprintf('loading data: %_s\n',options.training_data{i});
  tmp=load([options.input_path options.training_data{i}]);
%  X(:,:,i)=tmp.(options.training_data{i});
  data(:,:,i) = tmp.X;
end
  data=data(1:options.nTR,:,:);

%load testing data "monkey-dog" (TRAINING/TESTING)
load([options.input_path options.testing_data]);

clear tmp;

%%hyperalignment (TRAINING)
fprintf('start hyperalignment\n');

RG_filename = [options.exptype '_RG_data_' num2str(size(options.training_data,2)) '_' num2str(options.nTR)];
time_filename =  [options.exptype '_time_' num2str(size(options.training_data,2)) '_' num2str(options.nTR)];

if exist([options.output_path RG_filename '.mat'], 'file')
  fprintf('load RG data from file\n');
  load([options.output_path RG_filename '.mat']);
  load([[ options.output_path time_filename] '.mat']);
else
  fprintf('calculate RG\n');
  tic
  [R G]=alignment_handle(data);
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
                      '-s 1 -t 0 -n 0.5 -p 0.001 -b 1');
  [p1 p2 p3] = svmpredict(lblid(idx_tst),XRotate(idx_tst,:),svm_model,'-b 1');
  
  accuracy = [ accuracy sum(p1 == lblid(idx_tst))/length(p1)] ;
end

accuracy 
time_spent
accuracy_filename =  [options.exptype '_accuracy_' num2str(size(options.training_data,2)) '_' num2str(options.nTR)];

save( [ options.output_path accuracy_filename], 'accuracy');
save( [ options.output_path time_filename], 'time_spent');

cd(origin_path)
%diary off;
end


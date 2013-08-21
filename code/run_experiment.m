clear

%%set parameters
fprintf('loading options\n')
options.exptype = 'CFCH2_PRELIM' % PRELIM: preliminary or LOO: leave one out
alignment_handle = @CFCH2;

options.time = clock
options.time = [date '-' num2str(options.time(4)) num2str(options.time(5))]
options.input_path = '../data/input/'; 
options.working_path = '../data/working/' ; 
options.output_path = '../data/output/' ; 
options.random_seed = 1;
options.training_data = {'t_cb_vt_rh','t_dm_vt_rh'}%,'t_hj_vt_rh','t_kd_vt_rh',...
%                    't_kl_vt_rh','t_mh_vt_rh','t_ph_vt_rh','t_rb_vt_rh',...
%                   't_se_vt_rh','t_sm_vt_rh'};

options.testing_data = 'mkdg_vt_tlrc' ;

options.nTR = 400;

libsvm_path = '/mnt/cd/ramadge/pohsuan/libsvm-3.17-rondo/matlab/';
origin_path = '/mnt/cd/ramadge/pohsuan/hyperalignment/code/';

%libsvmpath = '../../libsvm-3.17-mac/matlab'
%origin_path = '../../hyperalignment/code/';



% create diary to save log
filename = [options.time '_diary.txt'];
eval(['diary ',filename]);
diary on;

rng(options.random_seed);

%load alignment data "raider of the lost track" (TRAINING)
for i=1:size(options.training_data,2),
  fprintf('loading data: %s\n',options.training_data{i});
  tmp=load([options.input_path options.training_data{i}]);
  X(:,:,i)=tmp.(options.training_data{i});
end
  X=X(1:options.nTR,:,:);

%load testing data "monkey-dog" (TRAINING/TESTING)
load([options.input_path options.testing_data]);

clear tmp;

%%hyperalignment (TRAINING)
fprintf('start hyperalignment\n');

RG_filename = [options.exptype '_RG_data_' num2str(size(options.training_data,2))];

if exist([options.output_path RG_filename '.mat'], 'file')
  fprintf('load RG data from file\n');
  load([options.output_path RG_filename '.mat']);
else
  fprintf('calculate RG\n');
  tic
  [R G]=alignment_handle(X);
  toc
  save( [ options.output_path RG_filename ], 'R', 'G');
end

error('stop')

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

correctNum =0; allCount =0;
for subj=1:size(options.training_data,2),
  fprintf('SVM subject:%d\n',subj);
  idx_trn=find(subjid~=subj);
  idx_tst=find(subjid==subj);

  svm_model=svmtrain2(lblid(idx_trn),XRotate(idx_trn,:),...
                      '-s 1 -t 0 -n 0.5 -p 0.001 -b 1');
  [p1 p2 p3] = svmpredict(lblid(idx_tst),XRotate(idx_tst,:),svm_model,'-b 1');
  
  correctNum = correctNum + sum(p1 == lblid(idx_tst));
  allCount = allCount + length(p1);
end

bsc_acc = correctNum/allCount

cd(origin_path)
diary off;
%for i=1:size(options.training_data,2),
%  idx = find(subjid!=i);
%  for j=1:numel(idx)
%    training_data(j,:)=XRtst(idx(j),:)*R(:,:,subjid(idx(j)));
%  end
%  svm_model = svm_train(@@@@);
%  end


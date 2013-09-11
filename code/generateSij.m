clear

nTR    = [100 200 400 800];
nvoxel = [100 200 400 800];

for a = 1:numel(nvoxel)
  for b = 1:numel(nTR)
    fprintf('nTR %d nvoxel %d',nTR(a), nvoxel(b));
    options.input_path = ['/mnt/cd/fastscratch/pohsuan/X_matrix_vs' num2str(nvoxel(a)) '_' num2str(nTR(b)) '/'];
    working_path = ['/mnt/cd/fastscratch/pohsuan/S_matrix_vs' num2str(nvoxel(a)) '_' num2str(nTR(b)) '/'];

    if exist(working_path,'dir')
      continue
    end

    mkdir(working_path)

    options.training_data = {'t_cb_vt_rh','t_dm_vt_rh','t_hj_vt_rh','t_kd_vt_rh',...
                    't_kl_vt_rh','t_mh_vt_rh','t_ph_vt_rh','t_rb_vt_rh',...
                   't_se_vt_rh','t_sm_vt_rh'};

    nsubj  = size(options.training_data,2);

    for i=1:nsubj,
      fprintf('load subject %d\n',i);
      load([options.input_path options.training_data{i}]);
      %input=load([options.input_path options.training_data{i} num2str(nvoxel)]);
      %X=input.(options.training_data{i});
      X=X(1:nTR(b),:);
      for j=1:nTR(b),
%        fprintf('generate S%d%d \n', i,j);
        S = X(j,:)'*X(j,:);
        save( [working_path 'S' num2str(i) num2str(j) ], 'S'); 
      end
    end
  end
end


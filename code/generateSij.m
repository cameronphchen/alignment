clear

working_path = '/mnt/cd/fastscratch/pohsuan/S_matrix/'
options.input_path = '../data/input/'; 


options.training_data = {'t_cb_vt_rh','t_dm_vt_rh','t_hj_vt_rh','t_kd_vt_rh',...
                    't_kl_vt_rh','t_mh_vt_rh','t_ph_vt_rh','t_rb_vt_rh',...
                   't_se_vt_rh','t_sm_vt_rh'};


  nTR    = 400;
  nvoxel = 2997;
  nsubj  = size(options.training_data,2);

  for i=7:nsubj,
    fprintf('load subject %d\n',i);
    input=load([options.input_path options.training_data{i}]);
    X=input.(options.training_data{i});
    X=X(1:nTR,:);
    for j=1:nTR,
      fprintf('generate S%d%d \n', i,j);
      S = X(j,:)'*X(j,:);
      save( [working_path 'S' num2str(i) num2str(j) ], 'S'); 
    end
  end


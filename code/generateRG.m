clear 

options.working_path = '/mnt/cd/fastscratch/pohsuan/'
options.input_path = '../data/input/'; 

options.training_data = {'t_cb_vt_rh','t_dm_vt_rh','t_hj_vt_rh','t_kd_vt_rh',...
                    't_kl_vt_rh','t_mh_vt_rh','t_ph_vt_rh','t_rb_vt_rh',...
                   't_se_vt_rh','t_sm_vt_rh'};
options.random_seed = 1;
rng(options.random_seed);

nTR    = 400;
nvoxel = 2997;
nsubj  = size(options.training_data,2);
nv2 = nvoxel*nvoxel;

R = zeros(nv2,nv2);
for i=1:nsubj,
  fprintf('generate random orthgonal matrix R%d\n',i);
  R=orth(randn(nv2,nv2));
  save([options.working_path 'RG/' 'R' num2str(i)], 'R');
end

G = zeros(nTR,nv2);

fprintf('calculate G\n');
for i=1:nsubj,
  fprintf('including S%d',i);
  S = zeros(nTR,nv2);
  for j=1:nTR,
    fprintf('loading S%d%d',i,j);
    load([options.working_path 'S' num2str(i) '/' 'S' num2str(i) num2str(j)]);
    S(j,:) = tmp;
    clear tmp;
  end
  fprintf('loading R%i\n');
  load([options.working_path 'RG/' 'R' num2str(i)]);
  fprintf('calculating G +  S*R\n'); 
  G = G +  S*R;
end

save([otions.working_path 'RG/' 'G' num2str(numbj)],'G');

G=G/nsubj;
G_tmp = zeros(nTR,nv2);



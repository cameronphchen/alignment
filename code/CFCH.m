%function [Q,G] =  CFCH()

%  nTR    = size(X,1);
%  nvoxel = size(X,2);
%  nsubj  = size(X,3);
  tic
  nTR    = 100
  nvoxel = 100
  nsubj  = 10
  N=1
  %use generateSij to precompute S 
   S_datapath = ['/mnt/cd/fastscratch/pohsuan/S_matrix_vs' num2str(nvoxel)]
  G_datapath = ['/mnt/cd/fastscratch/pohsuan/G_matrix_vs' num2str(nvoxel) num2str(N) '_' num2str(nsubj) '_' num2str(nTR) '/']
  Q_datapath = ['/mnt/cd/fastscratch/pohsuan/Q_matrix_vs' num2str(nvoxel) num2str(N) '_' num2str(nsubj) '_' num2str(nTR) '/']
  setting_path = ['/mnt/cd/fastscratch/pohsuan/setting_vs'num2str(nvoxel) num2str(N) '_' num2str(nsubj) '_' num2str(nTR) '/']
  mkdir(G_datapath)
  mkdir(Q_datapath)
  mkdir(setting_path)
  
  clear X
  iter = 0;
  if exist([setting_path 'iter.mat'], 'file')
    load([setting_path 'iter.mat']);
  end

  options.random_seed = 1;
  rng(options.random_seed);

  if iter~=0
    fprintf('load Q data from file\n');
    tmp=load([Q_datapath 'iter' num2str(iter) 'Q' '.mat']);
    Q =tmp.Q;
  else
    if exist([Q_datapath 'iter0Q.mat'], 'file')
        fprintf('load Q data from file\n');
        tmp=load([Q_datapath 'iter' num2str(iter) 'Q' '.mat']);
        Q =tmp.Q;
    else
        fprintf(['generate random Q\n']);
        Q = zeros(nvoxel,nvoxel,nsubj);
        for i=1:nsubj,
          fprintf('generating Q%d\n',i)
          Q(:,:,i) = orth(randn(nvoxel,nvoxel));
        end
        save( [ Q_datapath 'iter' num2str(iter) 'Q' ],'Q');
    end
  end
  G_diff = 0;

  
  if iter~=0
    fprintf('G exists, Use G data from file\n');
  else
    for t=1:nTR
      if ~exist([ G_datapath 'iter0G' num2str(t) '.mat' ], 'file')
          fprintf(['calculate G' num2str(t) '\n']);
          G = zeros(nvoxel,nvoxel);
          for i=1:nsubj,
            fprintf(['loading data' 'S' num2str(i) num2str(t) '\n']);
            load([S_datapath 'S' num2str(i) num2str(t)]);
            G = G + Q(:,:,i)'*S*Q(:,:,i);
          end
            G=G/nsubj;
          save( [ G_datapath 'iter' num2str(iter) 'G' num2str(t) ],'G');
      else
          fprintf('load G from file\n');
      end
    end
  end
  
  obj_val = 0;
  for t=1:nTR
     %fprintf('calculate obj_val, t= %d \n',t);
     load([G_datapath 'iter' num2str(iter) 'G' num2str(t) ]);
     G_diff = G_diff + norm(G,'fro');
     for i=1:nsubj               
        load([S_datapath  'S' num2str(i) num2str(t)]);
        obj_val = obj_val + trace(Q(:,:,i)'*S*Q(:,:,i)*G);
     end
  end


  fprintf('G_diff %f\n', G_diff);
  while G_diff>1e-2
    fprintf('=======================================\n');
    fprintf('%dth iteration : ',iter);
    fprintf('obj_val : %d \n', obj_val);
    fprintf('G_diff %f\n', G_diff);

    for t=1:nTR
      load([G_datapath 'iter' num2str(iter) 'G' num2str(t) ]);
      G_tmp = G;
      save( [ G_datapath 'iter' num2str(iter) 'G_tmp' num2str(t) ],'G_tmp');
    end

    for i=1:nsubj,
      fprintf('Calculate D%d \n',i);
      [U,D,V]=svd(Q(:,:,i),'econ');
%      fprintf('Calculate SG\n');
      SG = zeros(nvoxel,nvoxel);
      for t=1:nTR
        %fprintf('t=%d \n',t);
        load([S_datapath 'S' num2str(i) num2str(t)]);
        load([G_datapath 'iter' num2str(iter) 'G' num2str(t) ]);
        SG = SG + (U'*S*U) .* (V'*G*V);
        %assert(norm(SG-SG','fro')<1e-3,'SG not symmetric')
      end
      %SOLVE SDP PROBLEM WITH SG AND GET D, size nvoxel*nvoxel
%      fprintf('Calculate SG Done\n');
      [U_sg D_sg V_sg] = svd(SG,'econ');
      SD=U_sg(:,1);
      SD= sign(SD(1))*SD;
      D=diag(sign(SD));

      W = eye(nvoxel,nvoxel);
      for n=1:N %N is the amount of pivot
%        fprintf('%dth iteration / %d for rotation pivot\n',n,N);
        GS = zeros(nvoxel,nvoxel);
%        fprintf('Calculate GS\n');
        for t=1:nTR
          %fprintf('t=%d \n',t);
          load([S_datapath 'S' num2str(i) num2str(t)]);
          load([G_datapath 'iter' num2str(iter) 'G' num2str(t) ]);
          GS = GS + (D*V'*G*V*D')*(W'*U'*S*U*W);
        end
%        fprintf('Calculate GS Done\n');

        %[p q] = argmax |GS(p,q)-GS(q,p)|
        GS_tmp = GS - GS';

        [C I] = max(GS_tmp);
        [C q] = max(max(GS_tmp));
        p = I(q);
        %solve partial f / partial theta =0 to get 4 roots
        eq = zeros(nvoxel,1);
        eq(q,1) = 1;
        ep = zeros(nvoxel,1);
        ep(p,1) = 1;
        Z = eq*ep' - ep*eq';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
         fun_partial_2 = @(x) 0;
         fun_2 = @(x) 0;
        
         for t=1:nTR
%           fprintf('%d th for generate f\n',t);
           load([S_datapath  'S' num2str(i) num2str(t)]);
           load([G_datapath 'iter' num2str(iter) 'G' num2str(t) ]);
           S_t = W'*U'*S*U*W;
           G_t = D*V'*G*V*D';
           fun_tmp1=@(x) trace((cos(x)*Z+sin(x)*Z*Z)'*S_t ...
                              *(eye(nvoxel)+sin(x)*Z+(1-cos(x))*Z*Z)*G_t ...
                              +  (eye(nvoxel)+sin(x)*Z+(1-cos(x))*Z*Z)'*S_t...
                              *(cos(x)*Z+sin(x)*Z*Z)*G_t);
           fun_partial_2 = @(x) fun_partial_2(x) + fun_tmp1(x);
           fun_tmp2=@(x) trace((eye(nvoxel)+sin(x)*Z+(1-cos(x))*Z*Z)'*S_t...
                              *(eye(nvoxel)+sin(x)*Z+(1-cos(x))*Z*Z)*G_t);
           fun_2 = @(x) fun_2(x) + fun_tmp2(x);
         end
         fprintf('finish calculate the function, start to find zeros\n');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%        fprintf('generate function f\n');
        SG = zeros(nvoxel,nvoxel);
        ZSG = zeros(nvoxel,nvoxel);
        ZZSG = zeros(nvoxel,nvoxel);
        ZSZG= zeros(nvoxel,nvoxel);
        ZZSZG= zeros(nvoxel,nvoxel);
        ZSZZG= zeros(nvoxel,nvoxel);
        ZZSZZG= zeros(nvoxel,nvoxel);
        SZG= zeros(nvoxel,nvoxel);
        SZZG= zeros(nvoxel,nvoxel);
        ZZSZZG= zeros(nvoxel,nvoxel);
     
        ZZ = Z*Z;
        for t=1:nTR
          %fprintf('t=%d for preparing f(theta)\n',t);
          load([S_datapath 'S' num2str(i) num2str(t)]);
          load([G_datapath 'iter' num2str(iter) 'G' num2str(t) ]);
          %fprintf('1');
          St = (W'*U'*S*U*W);
          Gt = (D*V'*G*V*D');
          %fprintf('2');
          SG = SG + St*Gt;
          ZSG = ZSG + Z'*St*Gt;
          ZZSG = ZZSG + ZZ'*St*Gt;
          %fprintf('3');
          ZSZG = ZSZG + Z'*St*Z*Gt;
          ZZSZG = ZZSZG + ZZ'*St*Z*Gt;
          %fprintf('4');
          ZSZZG = ZSZZG + Z'*St*ZZ*Gt;
          ZZSZZG = ZZSZZG + ZZ'*St*ZZ*Gt;
          %fprintf('5',t);
          SZG = SZG + St*Z*Gt;
          SZZG = SZZG + St*ZZ*Gt;
        end
          SG = trace(SG);
          ZSG = trace(ZSG);
          ZZSG = trace(ZZSG);
          ZSZG = trace(ZSZG);
          ZZSZG = trace(ZZSZG);
          ZSZZG = trace(ZSZZG);
          ZZSZZG = trace(ZZSZZG);
          SZG = trace(SZG);
          SZZG = trace(SZZG);

         fun_partial = @(x) cos(x)*ZSG + sin(x)*ZZSG + sin(x)*cos(x)*ZSZG ...
                   + sin(x)^2*ZZSZG + cos(x)*(1-cos(x))*ZSZZG + sin(x)*(1-cos(x))*ZZSZZG...
                   + cos(x)*SZG + sin(x)*cos(x)*ZSZG + (1-cos(x))*cos(x)*ZZSZG...
                   + sin(x)*SZZG + sin(x)^2*ZSZZG + (1-cos(x))*sin(x)*ZZSZZG;

          fun = @(x) SG + sin(x)*ZSG + (1-cos(x))*ZZSG ...
                   + sin(x)*SZG + sin(x)^2*ZSZG + (1-cos(x))*sin(x)*ZZSZG...
                   + (1-cos(x))*SZZG + sin(x)*(1-cos(x))*ZSZZG + (1-cos(x))*(1-cos(x))*ZZSZZG;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         theta = findAllZeros(fun_partial,4);%%
%         fprintf('finish finding the zeros\n');
         max_val = fun(theta(1));%%
         max_idx = 1;
         for y=2:size(theta,2)
            if( fun(theta(y))> max_val)%%
              max_idx = y;
              max_val = fun(theta(y));%%
            end
         end
         theta = theta(max_idx);
%         fprintf('calculate W \n');
         W = W * (eye(nvoxel,nvoxel)+sin(theta)*Z + (1-cos(theta))*Z*Z);
      end
      Q(:,:,i) = U*W*D*V';
    end

    %Update G
    G_diff=0;
    fprintf('update G\n');
    for t=1:nTR,
    G = zeros(nvoxel,nvoxel);
      for i=1:nsubj,
        load([S_datapath 'S' num2str(i) num2str(t)]);
        G = G + Q(:,:,i)'*S*Q(:,:,i);
      end
      G = G/nsubj;
      save( [ G_datapath 'iter' num2str(iter+1) 'G' num2str(t) ],'G');
      load( [ G_datapath 'iter' num2str(iter) 'G_tmp' num2str(t) ],'G_tmp');
      G_diff = G_diff + norm(G-G_tmp,'fro');
    end
    save( [ Q_datapath 'iter' num2str(iter+1) 'Q' ],'Q');

    iter = iter+1;
    save( [ setting_path 'iter'],'iter');
    obj_val = 0;
    for t=1:nTR
      load([G_datapath 'iter' num2str(iter) 'G' num2str(t) ]);
      for i=1:nsubj
        load([S_datapath  'S' num2str(i) num2str(t)]);
        obj_val = obj_val + trace(Q(:,:,i)'*S*Q(:,:,i)*G);
      end
    end


  end
  
  fprintf('finish CFCH alignment');
  toc

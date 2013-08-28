%function [Q,G] =  CFCH()

%  nTR    = size(X,1);
%  nvoxel = size(X,2);
%  nsubj  = size(X,3);
  nTR    = 3;
  nvoxel = 2997;
  nsubj  = 3;
  N=3 
  %use generateSij to precompute S 
  S_datapath = '/mnt/cd/fastscratch/pohsuan/S_matrix/'
  G_datapath = '/mnt/cd/fastscratch/pohsuan/G_matrix/'
  Q_datapath = '/mnt/cd/fastscratch/pohsuan/Q_matrix/'

  clear X

  options.random_seed = 1;
  rng(options.random_seed);

  if exist([Q_datapath 'Q' '.mat'], 'file')
    fprintf('load Q data from file\n');
    tmp=load([Q_datapath 'Q' '.mat']);
    Q =tmp.Q;
  else
    fprintf(['generate random Q\n']);
    Q = zeros(nvoxel,nvoxel,nsubj);
    for i=1:nsubj,
      fprintf('generating Q%d\n',i)
      Q(:,:,i) = orth(randn(nvoxel,nvoxel));
    end
      save( [ Q_datapath 'Q'],'Q');
  end
  G_diff = 0;
  for t=1:nTR,
    if exist([G_datapath 'G' num2str(t) '.mat'], 'file')
      fprintf('G%d exists, Use G data from file\n',t);
      load([G_datapath 'G' num2str(t) ]);
    else
      fprintf(['calculate G' num2str(t) '\n']);
      G = zeros(nvoxel,nvoxel);
      for i=1:nsubj,
        fprintf(['loading data' 'S' num2str(i) num2str(t) '\n']);
        load([S_datapath 'S' num2str(i) num2str(t)]);
        G = G + Q(:,:,i)'*S*Q(:,:,i);
      end
      save( [ G_datapath 'G' num2str(t) ],'G');
    end
    G_diff = G_diff + norm(G,'fro');
  end
  

  %Q_tmp=zeros(nvoxel,nvoxel,nsubj);
  fprintf('G_diff %f\n', G_diff);
  count = 0;
  while G_diff>1e-2
    count = count+1;
    fprintf('=======================================');
    fprintf('%dth iteration : ',count);
    fprintf('G_diff %f\n', G_diff);
    %Q_tmp = Q;
    for t=1:nTR
      G_tmp = G;
      save( [ G_datapath 'G_tmp' num2str(t) ],'G_tmp');
    end


    for i=1:nsubj,
      fprintf('Calculate D%d \n',i);
      [U,D,V]=svd(Q(:,:,i),'econ');
      SG = zeros(nvoxel,nvoxel);
      fprintf('Calculate SG\n');
      for t=1:nTR
        fprintf('t=%d \n',t);
        load([S_datapath 'S' num2str(i) num2str(t)]);
        load([G_datapath 'G' num2str(t) ]);
        SG = SG + U'*S*U .* V'*G*V;
      end
      fprintf('Finshed Calculate SG\n'); 
        %SOLVE SDP PROBLEM WITH SG AND GET D, size nvoxel*nvoxel
      fprintf('solve sdp\n');
%      cvx_begin sdp
%        variable D_tilde(nvoxel,nvoxel) symmetric
%        maximize (trace(SG'*D_tilde))
%        for i=1:nvoxel
%          D_tilde(i,i) == 1
%        end
%        D_tilde>=0
%      cvx_end
      
      SD=SG(:,1);
      D=diag(sign(SD));
      
      W = eye(nvoxel,nvoxel);
      for n=1:N %N is the amount of pivot
        fprintf('%dth iteration / %d for rotation pivot\n',n,N);
        GS = zeros(nvoxel,nvoxel);
        fprintf('Calculate GS\n');
        for t=1:nTR
          fprintf('t=%d \n',t);
          load([S_datapath 'S' num2str(i) num2str(t)]);
          load([G_datapath 'G' num2str(t) ]);
          GS = GS + (D*V'*G*V*D')*(W'*U'*S*U*W);
        end
        %[p q] = argmax |GS(p,q)-GS(q,p)|
        [C I] = max(GS);
        [C q] = max(max(GS));
        p = I(q);
        %solve partial f / partial theta =0 to get 4 roots
        eq = zeros(nvoxel,1);
        eq(q,1) = 1;
        ep = zeros(nvoxel,1);
        ep(p,1) = 1;
        Z = eq*ep' - ep*eq';
        fprintf('generate function f\n');
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
          fprintf('t=%d for preparing f(theta)\n',t);
          load([S_datapath 'S' num2str(i) num2str(t)]);
          load([G_datapath 'G' num2str(t) ]);
          fprintf('1',t);
          St = (W'*U'*S*U*W);
          Gt = (D*V'*G*V*D');
          fprintf('2',t);
          ZSG = ZSG + Z*St*Gt;
          ZZSG = ZZSG + ZZ*St*Gt;
          fprintf('3',t);
          ZSZG = ZSZG + Z*St*Z*Gt;
          ZZSZG = ZZSZG + ZZ*St*Z*Gt;
          fprintf('4',t);
          ZSZZG = ZSZZG + Z*St*ZZ*Gt;
          ZZSZZG = ZZSZZG + ZZ*St*ZZ*Gt;
          fprintf('5',t);
          SZG = SZG + St*Z*Gt;
          SZZG = SZZG + St*ZZ*Gt;
          ZZSZZG = ZZSZZG + ZZ*St*ZZ*Gt;
        end
          ZSG = trace(ZSG);
          ZZSG = trace(ZZSG);
          fprintf('3',t);
          ZSZG = trace(ZSZG);
          ZZSZG = trace(ZZSZG);
          fprintf('4',t);
          ZSZZG = trace(ZSZZG);
          ZZSZZG = trace(ZZSZZG);
          fprintf('5',t);
          SZG = trace(SZG);
          SZZG = trace(SZZG);
          ZZSZZG = trace(ZZSZZG);




         func = @(x) cos(x)*ZSG + sin(x)*ZZSG + sin(x)*cos(x)*ZSZG ...
                   + sin(x)^2*ZZSZG + cos(x)*(1-cos(x))*ZSZZG + sin(x)*(1-cos(x))*ZZSZZG...
                   + cos(x)*SZG + sin(x)*cos(x)*ZSZG + (1-cos(x))*cos(x)*ZZSZG...
                   + sin(x)*SZZG + sin(x)^2*ZSZZG + (1-cos(x))*sin(x)*ZZSZZG;
        

         theta = findAllZeros(func,4);
         display(theta); 
         %theta_star = argmax_theta f()
         max_val = func(theta(1));
         max_idx = 1;
         for y=2:size(theta,2)
            if( func(theta(y))> max_val)
              max_idx = y;
              max_val = func(theta(y));
            end
         end
         theta = theta(max_idx);

         W = W * exp(theta * Z);
      end
      Q(:,:,i) = U'*W*D*V;
    end

    %Update G
    G_diff=0;
    fprintf('update G\n');
    for t=1:nTR,
    G = zeros(nvoxel,nvoxel);
      for i=1:nsubj,
        fprintf(['loading data' 'S' num2str(i) num2str(t) '\n']);
        load([S_datapath 'S' num2str(i) num2str(t)]);
        G = G + Q(:,:,i)'*S*Q(:,:,i);
      end
      save( [ G_datapath 'G' num2str(t) ],'G');
      load( [ G_datapath 'G_tmp' num2str(t) ],'G_tmp');
      G_diff = G_diff + norm(G-G_tmp);
    end

  end
  save( [ Q_datapath 'Q' ],'Q');
  fprintf('finish CFCH alignment');


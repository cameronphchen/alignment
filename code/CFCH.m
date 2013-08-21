function [Q,G] =  CFCH2(X)

%  nTR    = size(X,1);
%  nvoxel = size(X,2);
%  nsubj  = size(X,3);
  nTR    = 50;
  nvoxel = 2997;
  nsubj  = 10;
  N=10 
  %use generateSij to precompute S 
  S_datapath = '/mnt/cd/fastscratch/pohsuan/S_matrix/'
  G_datapath = '/mnt/cd/fastscratch/pohsuan/G_matrix/'

  clear X

  options.random_seed = 1;
  rng(options.random_seed);

  Q = zeros(nvoxel,nvoxel,nsubj);
  for i=1:nsubj,
    Q(:,:,i) = orth(randn(nvoxel,nvoxel));
  end
  for t=1:nTR,
    G = zeros(nvoxel,nvoxel);
    for i=1:nsubj,
      fprintf(['loading data' 'S' num2str(i) num2str(t)]);
      load([S_datapath 'S' num2str(i) num2str(t)]);
      G = G + Q(:,:,i)'*S*Q(:,:,i);
    end
      save( [ G_datapath 'G' num2str(t) ],'G');
  end


  Q_tmp=zeros(nvoxel,nvoxel,nsubj);

  count = 0;
  while norm(Q-Q_tmp,'fro')>1e-3 
    count = count+1;
    fprintf('%dth iteration : ',count);
    fprintf('Q-Q_tmp %f\n', norm(Q-Q_tmp,'fro'));
    for i=1:nsubj,
      [U,D,V]=svd(Q(:,:,i),'econ');
      for t=1:nTR
        load([S_datapath 'S' num2str(i) num2str(t)]);
        load([[ G_datapath 'G' num2str(t) ]);
        SG = SG + U'*S*U .* V'*G*V; 
        %SOLVE SDP PROBLEM WITH SG AND GET D, size nvoxel*nvoxel
        cvx_begin sdp
          variable D_tilde(nvoxel,nvoxel) symmetric
          maximize (trace(M'*D_tilde))
          for i=1:nvoxel
            D_tilde(i,i) == 1
          end
          D_tilde>=0
        cvx_end
        SD=D_tilde(:,1);
        D=diag(sign(SD));
      end
      W = eye(nvoxel,nvoxel);
      for n=1:N %N is the amount of pivot
        GS = zeros(nvoxel,nvoxel);
        for t=1:nTR
          load([S_datapath 'S' num2str(i) num2str(t)]);
          load([G_datapath 'G' num2str(t) ]);
          GS = GS + (D*V'*G*V*D')*(W'*U'*S*U*W);
        end
         %[p q] = argmax |GS(p,q)-GS(q,p)|
         [C I] = max(GS);
         [C q] = max(max(GS));
         p = I(q);
         %solve partial f / partial theta =0 to get 4 roots
         func = @(x) 
         theta = findAllZeros(func,4);
 
         %theta_star = argmax_theta f()

         eq = zeros(nvoxel,1);
         eq(q,1) = 1;
         ep = zeros(nvoxel,1);
         ep(p,1) = 1;
         W = W * exp(theta * (eq*ep' - ep*eq'));
      end
      Q(:,:,i) = U'*W*D*V;
    end

    %Update G
    for t=1:nTR,
    G = zeros(nvoxel,nvoxel);
      for i=1:nsubj,
        fprintf(['loading data' 'S' num2str(i) num2str(t)]);
        load([S_datapath 'S' num2str(i) num2str(t)]);
        G = G + Q(:,:,i)'*S*Q(:,:,i);
      end
      save( [ G_datapath 'G' num2str(t) ],'G');
    end

  end
  fprintf('finish CFCH alignment');



clear

  nTR    = 3;
  nvoxel = 5;
  nsubj  = 3;
  N=1
  %use generateSij to precompute S 
  S_datapath = '/mnt/cd/fastscratch/pohsuan/tmp_S_matrix/'
  G_datapath = '/mnt/cd/fastscratch/pohsuan/tmp_G_matrix/'
  Q_datapath = '/mnt/cd/fastscratch/pohsuan/tmp_Q_matrix/'

  clear X

  options.random_seed = 1;
  rng(options.random_seed);
  
  X = zeros(nTR,nvoxel,nsubj);
  S = zeros(nvoxel,nvoxel,nsubj,nTR); 
  Q = zeros(nvoxel,nvoxel,nsubj);
  G = zeros(nvoxel,nvoxel,nTR);

  for i=1:nsubj,
    fprintf('generating Q%d\n',i)
%    Q(:,:,i) = orth(randn(nvoxel,nvoxel));
    Q(:,:,i) = eye(nvoxel,nvoxel);

    X(:,:,i) = randn(nTR,nvoxel);
    for t=1:nTR
      S(:,:,i,t) = X(t,:,i)'*X(t,:,i);
    end
  end
  obj_val = zeros(nsubj,1);
  G_diff = 0;
  for t=1:nTR,
    fprintf(['calculate G' num2str(t) '\n']);
    for i=1:nsubj,
      G(:,:,t) = G(:,:,t) + Q(:,:,i)'*S(:,:,i,t)*Q(:,:,i);
    end
    G_diff = G_diff + norm(G(:,:,t),'fro');
  end
  
  for i=1:nsubj
     for t=1:nTR
        obj_val(i) = obj_val(i) + trace(Q(:,:,i)'*S(:,:,i,t)*Q(:,:,i)*G(:,:,t));
     end
  end
  

  fprintf('G_diff %f\n', G_diff);
  count = 0;
  while G_diff>1e-2
    count = count+1;
    fprintf('=======================================\n');
    fprintf('%dth obj_val=%d %d, %d, %d\n',count,sum(obj_val),obj_val(1),obj_val(2),obj_val(3));
    fprintf('G_diff %f\n', G_diff);

    G_tmp = G;

    for i=1:nsubj,
      fprintf('----subject %d----\n',i);
      %fprintf('Calculate D%d \n',i);
      [U,D,V]=svd(Q(:,:,i),'econ');
      SG = zeros(nvoxel,nvoxel);
      %fprintf('Calculate SG\n');
      for t=1:nTR
        SG = SG + (U'*S(:,:,i,t)*U) .* (V'*G(:,:,t)*V);
        assert(norm(SG-SG','fro')<1e-3,'SG not symmetric')
      end
      %fprintf('Finshed Calculate SG\n'); 
      %fprintf('solve sdp\n');
      [U_sg D_sg V_sg] = svd(SG,'econ');
      SD=U_sg(:,1);
      SD= sign(SD(1))*SD;
      D=diag(sign(SD));
      
      W = eye(nvoxel,nvoxel);
      for n=1:N %N is the amount of pivot
        %fprintf('%d/%d for rotation pivot\n',n,N);
        GS = zeros(nvoxel,nvoxel);
        %fprintf('Calculate GS\n');
        for t=1:nTR
          GS = GS + (D*V'*G(:,:,t)*V*D')*(W'*U'*S(:,:,i,t)*U*W);
        end
        GS_tmp = GS - GS';
        %[p q] = argmax |GS(p,q)-GS(q,p)|
        %[sortedValues,sortIndex] = sort(GS_tmp(:),'descend');
        %idx_count = 1;
        %p = mod((sortIndex(idx_count)-1),nvoxel)+1;
        %q = floor((sortIndex(idx_count)-1)/nvoxel)+1;
        %while(p==q) 
        %    idx_count = idx_count +1;
        %    p = mod((sortIndex(idx_count)-1),nvoxel)+1;
        %    q = floor((sortIndex(idx_count)-1)/nvoxel)+1;
        %end
        [C I] = max(GS_tmp);
        [C q] = max(max(GS_tmp));
        p = I(q);
        %solve partial f / partial theta =0 to get 4 roots
        eq = zeros(nvoxel,1);
        eq(q,1) = 1;
        ep = zeros(nvoxel,1);
        ep(p,1) = 1;
        Z = eq*ep' - ep*eq';
        
         fun_partial = @(x) 0;
         for t=1:nTR
           S_t = W'*U'*S(:,:,i,t)*U*W;
           G_t = D*V'*G(:,:,t)*V*D';
           fun_tmp=@(x) trace((cos(x)*Z+sin(x)*Z*Z)'*S_t*(eye(nvoxel)+sin(x)*Z+(1-cos(x))*Z*Z)*G_t ...
                       +(eye(nvoxel)+sin(x)*Z+(1-cos(x))*Z*Z)'*S_t*(cos(x)*Z+sin(x)*Z*Z)*G_t);
           fun_partial = @(x) fun_partial(x) + fun_tmp(x);
         end
        
         fun = @(x) 0;
         for t=1:nTR
           S_t = W'*U'*S(:,:,i,t)*U*W;
           G_t = D*V'*G(:,:,t)*V*D';
           fun_tmp=@(x) trace((eye(nvoxel)+sin(x)*Z+(1-cos(x))*Z*Z)'*S_t*(eye(nvoxel)+sin(x)*Z+(1-cos(x))*Z*Z)*G_t);
           fun = @(x) fun(x) + fun_tmp(x);
         end
         
        
         %figure 
         %figure(1)
         %fplot(fun_partial,[0 2*pi],1000)
         %figure(2)
         %fplot(fun,[0 2*pi],1000)

         theta = findAllZeros(fun_partial,4)
         %display(theta); 
         %theta_star = argmax_theta f()
         max_val = fun(theta(1));
         max_idx = 1;
         for y=2:size(theta,2)
            if( fun(theta(y))> max_val)
              max_idx = y;
              max_val = fun(theta(y));
            end
         end
         theta = theta(max_idx);
         
         W = W * (eye(nvoxel,nvoxel)+sin(theta)*Z + (1-cos(theta))*Z*Z);
      end
      
      %
      obj_val_s1=0;
      for t=1:nTR
         obj_val_s1 = obj_val_s1 + trace(Q(:,:,i)'*S(:,:,i,t)*Q(:,:,i)*G(:,:,t));
      end
      
      Q(:,:,i) = U*W*D*V';
      
      %
      obj_val_s2=0;
      for t=1:nTR
        obj_val_s2 = obj_val_s2 + trace(Q(:,:,i)'*S(:,:,i,t)*Q(:,:,i)*G(:,:,t));
      end
      fprintf('sub %d obj_val difference: %d \n',i,obj_val_s2-obj_val_s1);

    end

    %Update G
    G_diff=0;
    obj_val = zeros(nsubj,1);
    G = zeros(nvoxel,nvoxel,nTR);
    %fprintf('update G\n');
    for t=1:nTR,
      for i=1:nsubj,
        G(:,:,t) = G(:,:,t) + Q(:,:,i)'*S(:,:,i,t)*Q(:,:,i);
      end
      G_diff = G_diff + norm(G(:,:,t)-G_tmp(:,:,t),'fro');
    end
    
    for t=1:nTR
      for i=1:nTR
        obj_val(i) = obj_val(i) + trace(Q(:,:,i)'*S(:,:,i,t)*Q(:,:,i)*G(:,:,t));
      end
    end

  end
  save( [ Q_datapath 'Q' ],'Q');
  fprintf('finish CFCH alignment');


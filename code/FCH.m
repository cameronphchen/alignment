function [R,G] =  FCH(X)
  
  fprintf('DO FCH\n');
  nTR    = size(X,1);
  nvoxel = size(X,2);
  nsubj  = size(X,3);

  fprintf('calculate S\n');
  for i=1:nsubj,
    S(:,:,i) = X(:,:,i)'*X(:,:,i);
  end

  R = zeros(nvoxel,nvoxel,nsubj);
  % random orthogonal matrix initialization
  for i=1:nsubj,
    R(:,:,i) = orth(randn(nvoxel,nvoxel));
  end
  
  fprintf('calculate G\n');
  G = zeros(nvoxel,nvoxel);
  for i=1:nsubj,
    G = G + R(:,:,i)'*S(:,:,i)*R(:,:,i);
  end
  G = G/nsubj;
  G_tmp=zeros(nvoxel,nvoxel);
  
  fprintf('start iterative solving\n');
  count = 0;
  while norm(G-G_tmp,'fro')>1e-3 
    count = count+1;
    fprintf('%dth iteration : ',count);
    fprintf('G-G_tmp %f\n', norm(G-G_tmp,'fro'));
%    display(norm(G-G_tmp,'fro'));
    G_tmp = G;
%    display(G_tmp)
    for i=1:nsubj,
      [UG,SG,VG]=svd(G,'econ');
      [US,SS,VS]=svd(S(:,:,i),'econ');
      
      R(:,:,i) = VS*VG';
    end
    G = zeros(nvoxel,nvoxel);
    for j=1:nsubj,
      G = G + R(:,:,i)'*S(:,:,i)*R(:,:,i);
    end
    G = G/nsubj;

%    display(G);
  end

  fprintf('finish hyperalignment : ');
  fprintf('G-G_tmp %f\n', norm(G-G_tmp,'fro'));
return


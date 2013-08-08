function [R,G] =  HA(s)

  nTR    = size(s,1);
  nvoxel = size(s,2);
  nsubj  = size(s,3);

  R = zeros(nvoxel,nvoxel,nsubj);
  % random orthogonal matrix initialization
  for i=1:nsubj,
    R(:,:,i) = orth(randn(nvoxel,nvoxel));
  end

  G = zeros(nTR,nvoxel);
  for i=1:nsubj,
    G = G + s(:,:,i)*R(:,:,i);
  end
  G = G/nsubj;
  G_tmp=zeros(nTR,nvoxel);

  count = 0;
  while norm(G-G_tmp,'fro')>1e-3 
    count = count+1;
    fprintf('%dth iteration : ',count);
    fprintf('G-G_tmp %f\n', norm(G-G_tmp,'fro'));
%    display(norm(G-G_tmp,'fro'));
    G_tmp = G;
%    display(G_tmp)
    for i=1:nsubj,
      [U,S,V]=svd(s(:,:,i)'*G,'econ');
      R(:,:,i) = U*V';
    end
    G = zeros(nTR,nvoxel);
    for j=1:nsubj,
      G = G + s(:,:,j)*R(:,:,j);
    end
    G = G/nsubj;

%    display(G);
  end
  fprintf('finish hyperalignment');
return


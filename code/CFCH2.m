function [R,G] =  CFCH2(X)

  nTR    = size(X,1);
  nvoxel = size(X,2);
  nsubj  = size(X,3);

  nv2 = nvoxel*nvoxel;

  S = zeros(nTR,nv2,nsubj);
  for i=1:nsubj,
    for j=1:nTR,
      fprintf('generate S%d%d \n', i,j);
      tmp = X(j,:,i)'*X(j,:,i);
      S(j,:,i) = tmp(:)';
    end
  end

  error

  R = zeros(nv2,nv2,nsubj);
  % random orthogonal matrix initialization
  for i=1:nsubj,
    R(:,:,i) = orth(randn(nv2,nv2));
  end

  G = zeros(nTR,nv2);
  for i=1:nsubj,
    G = G + S(:,:,i)*R(:,:,i);
  end
  G = G/nsubj;
  G_tmp=zeros(nTR,nv2);

  count = 0;
  while norm(G-G_tmp,'fro')>1e-3 
    count = count+1;
    fprintf('%dth iteration : ',count);
    fprintf('G-G_tmp %f\n', norm(G-G_tmp,'fro'));
%    display(norm(G-G_tmp,'fro'));
    G_tmp = G;
%    display(G_tmp)
    for i=1:nsubj,
      [U,M,V]=svd(S(:,:,i)'*G,'econ');
      R(:,:,i) = U*V';
    end
    G = zeros(nTR,nv2);
    for j=1:nsubj,
      G = G + S(:,:,j)*R(:,:,j);
    end
    G = G/nsubj;

%    display(G);
  end
  fprintf('finish hyperalignment');
return


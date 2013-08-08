function [R,G] =  StackSVD(X)

  nTR    = size(X,1);
  nvoxel = size(X,2);
  nsubj  = size(X,3);

  StackX = zeros(nTR, nvoxel*nsubj);
  fprintf('stackX\n');
  for i=1:nsubj,
    StackX(:,(1+nvoxel*(i-1)):nvoxel*i)=X(:,:,i);
  end
  
  [U,S,V] = svd(StackX,'econ');
  R = zeros(nTR,nvoxel,nsubj);

  for i=1:nsubj,
    display(size(R(:,:,i)))
    display(size(V((1+nvoxel*(i-1)):nvoxel*i,:)'))
    R(:,:,i) = V((1+nvoxel*(i-1)):nvoxel*i,:)';
  end

  G=U*S/nsubj;

  fprintf('Stack SVD done\n');
return


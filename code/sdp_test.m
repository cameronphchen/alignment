M= randn(5);
n=5

cvx_begin sdp
    variable D(n,n) symmetric
    maximize (trace(M'*D))
    for i=1:n
        D(i,i) == 1
        D(i,i) == 1
        D(i,i) == 1
    end
    D>=0
cvx_end
S2 = [1 2; 2 4;]
S1 = [9 12; 12 16;]
e = [1; 1;]


cvx_begin
  variable Q
  minimize trace(S2*Q'*S1*Q)
  subject to
    Q*e==e
    e'*Q==e'
cvx_end


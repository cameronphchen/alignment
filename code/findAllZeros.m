function x = findAllZeros( func, numroot )

lb = 0; % Set a lower bound for the function.
ub = 2*pi; % Set an upper bound for the function.
x = []; % Initializes x.

for i=lb:1:ub
%    fprintf('theta = %f\n',i);
    % Look for the zeros in the function's current time window.
    if ~isequal(sign(func(i-.01)), sign(func(i+1)))
      x = [x fzero(func, i,optimset('display','off'))];
    end
end

% Make sure that there are no duplicates.
x = unique(x);
DUPE = diff(x) < 1e-2;
x(DUPE) = [];
%assert(size(x,2) <= numroot);

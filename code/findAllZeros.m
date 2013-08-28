function x = findAllZeros( func, numroot )

lb = 0; % Set a lower bound for the function.
ub = 2*pi; % Set an upper bound for the function.
x = []; % Initializes x.

for i=lb:0.1:ub
    % Look for the zeros in the function's current time window.
    x = [x fzero(func, i)];
end

% Make sure that there are no duplicates.
x = unique(x);
DUPE = diff([x NaN]) < 1.5e-1;
x(DUPE) = [];
%assert(size(x,2) <= numroot);

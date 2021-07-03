% ==========================================================================
% Check if the input argument is invalid, which means it is either nan or 
%   inf
% ==========================================================================

function output = IsInvalid(x)
    output = isnan(x) | isinf(x);
end


function M = helper_mat(x, R, C)
% Takes a length R*C vector x, returns an NxM matrix M.  Undoes 'vec()'.

if ( length(x) ~= R*C )
    error('Matrix size inconsistent with provided vector');
end

M = zeros(R,C);

for c = 1:C
    for r = 1:R
        M(r,c) = x( r + R*(c-1) );
    end
end

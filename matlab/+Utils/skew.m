function S = skew(w)
%SKEW skew symmetric matrix
    assert(length(w) == 3, 'skew symmetric matrix for 3d vector only')
    S = [0, -w(3), w(2); ...
         w(3), 0.0, -w(1); ...
         -w(2), w(1), 0];
end


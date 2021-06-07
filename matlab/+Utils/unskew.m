function w = unskew(S)
%SKEW skew symmetric matrix
    assert(S == -S', 'input not a skew symmetric matrix')
    w = [S(3, 2); S(1, 3); S(2, 1)];
end



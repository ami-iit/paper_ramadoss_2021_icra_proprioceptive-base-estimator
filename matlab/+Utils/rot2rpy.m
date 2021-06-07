function rpy = rot2rpy(R)
arguments
    R (3, 3) double;
end
Rdyn = iDynTree.Rotation();
Rdyn.fromMatlab(R);
rpy = Rdyn.asRPY().toMatlab();
end


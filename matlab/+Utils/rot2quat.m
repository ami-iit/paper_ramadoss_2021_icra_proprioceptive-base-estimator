function q = rot2quat(R)
arguments
    R (3, 3) double;
end
Rdyn = iDynTree.Rotation();
Rdyn.fromMatlab(R);
q = Rdyn.asQuaternion().toMatlab();
end


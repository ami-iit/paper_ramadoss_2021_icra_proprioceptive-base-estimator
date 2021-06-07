function R = quat2rot(q)
    qdyn = iDynTree.Vector4();
    qdyn.fromMatlab(q);
    Rdyn = iDynTree.Rotation();
    Rdyn.fromQuaternion(qdyn);
    R = Rdyn.toMatlab();
end


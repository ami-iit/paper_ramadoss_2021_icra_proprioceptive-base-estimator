function R = rpy2rot(roll, pitch, yaw)
R = iDynTree.Rotation.RPY(roll, pitch, yaw).toMatlab();
end

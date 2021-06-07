function pose = idynPose(R,p)
pose = iDynTree.Transform;
Rdyn = iDynTree.Rotation();
Rdyn.fromMatlab(R);
pdyn = iDynTree.Position();
pdyn.fromMatlab(p)
pose.setRotation(Rdyn);
pose.setPosition(pdyn);
end


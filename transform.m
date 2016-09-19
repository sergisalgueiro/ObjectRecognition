function npose = transform(pose,oldPose)

T=[cos(pose(3)) -sin(pose(3)) 0 pose(1);
   sin(pose(3)) cos(pose(3)) 0 pose(2);
   0 0 1 0;
   0 0 0 1];
oPose=[oldPose 1]';
npose=T*oPose;
end
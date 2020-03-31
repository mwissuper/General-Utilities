function ThetaInDegrees = getAngVec(u,v,p)
% Get angle between two 3D vectors on the same plane 
% Project vectors to plane by remove element defined by p. Data in (x,y,z)
% (e.g. if p == 2, remove y component and project onto frontal plane with
% person facing y dir). Does not work well for small angles.
u(p) = 0; v(p) = 0;
ThetaInDegrees = atan2d(norm(cross(u,v)),dot(u,v));
% quiver([0,0],[0,0],u([1 3]),v([1 3]));
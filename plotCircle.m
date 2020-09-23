function h = plotCircle(x,y,d)

% plot circle centered at [x,y] with diameter d
h = rectangle('position',[x-d/2 y-d/2 d d],'curvature',1);
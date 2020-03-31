function [] = plotij(m,n,i,j)
% function [] = plotij(m,n,i,j)
% Makes the (i,j)th subplot of an (m,n) subplot figure active.  I have no idea
% why Matlab orders the subplots this way by default.

p = sub2ind([n,m],j,i);
subplot(m,n,p)
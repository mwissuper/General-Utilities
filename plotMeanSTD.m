function plotMeanSTD(ax,x,y,std,style,color)
%plotMeanSTD: Plots mean and standard deviation of a curve.

%   Luke Drnach
%   November 14, 2018

% First plot the standard deviation as a fill plot in the background
upper = y + std;
lower = y - std;
fill(ax,[x,x(end:-1:1)],[upper, lower(end:-1:1)],color,'EdgeColor','none','FaceAlpha',0.4);
hold on;
% Then plot the line
plot(x,y,style,'Color',color);

end


%% for lab meeting adobe illustrator tutorial
close all; clear
load('ForIllustratorTutorial.mat')

%% Getting a basic plot into illustrator
close all
figure(1)
plot(time,mean(Cz),'Color',[0 0 0]);
title('Cz EEG data for one participant')
ylabel('Cz-EEG (uV)')
xlabel('Time (s)')
xl=[-0.2 1];
xlim(xl)


%% something with more things we can format
figure(2);clf
yl1=[-60 30]; % y axis limits for EEG data
yl2=[-0.6 0.6]; % y axis limits for acceleration data
for_legend = [];
for x = 1:3 % for each perturbation magnitude
   inds = find(flags(:,3)==x); % index the perturbations of that magnitude
   plotij(2,1,1,1); plot(time,mean(Cz(inds,:)));hold on
   plotij(2,1,2,1); plot(time,mean(accel(inds,:)));hold on
   for_legend = ([for_legend,{['Level ', num2str(x)]}]);
end
plotij(2,1,1,1);ylim(yl1);xlim(xl);plot([0.1 0.1],yl1,'k');plot([0.2 0.2],yl1,'k')
plotij(2,1,2,1);ylim(yl2);xlim(xl);plot([0.1 0.1],yl2,'k');plot([0 0],yl2,'k')

legend(for_legend)

%% to send it to illustrator
print_all_figures_to_eps

% sometimes, if you have an unusual amount of data in your figure (e.g. all
% the single trial data from all your subjects in the same plot), matlab
% will decide to make it out of uneditable tiles, the work around for this
% is the following statement, but will only print the most recently active
% figure and you have to manually name it. 

%print -painters -depsc figureName.eps

%%
error % handy trick to force run by section if you hit run at the top

%% second example scatter plot
close all; clear
load('forDryad.mat')
clearvars -except Subject_data Subject_data_Names
%%
figure(1);
plot(Subject_data(:,3),Subject_data(:,4),'k.')

print_all_figures_to_eps

%% Full model 

close all; clear all; clc;

set(0, 'defaultaxesfontsize',16)

% initialize the model
Setup_StageI

%% Optimization
arSimu(true,true,true);

n_fit = 1000; % number of starts

if n_fit > 1
    arFitLHS(n_fit); % multistart optimization
else
    arFit;
end

% Print optimization results
arPrint;     % display parameter values
arSimu(true,true,true);

% %Profile likelihood
% arPLEInit
% ple


arPlot

figure(11)
plot(ar.model.condition(1).tFine,ar.model.condition(1).xFineSimu(:,1), '-b','LineWidth',2);
title('S')
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Susc','eps')
saveas(fig,'Susc','fig')

figure(12)
plot(ar.model.condition(1).tFine,ar.model.condition(1).xFineSimu(:,2), '-b','LineWidth',2);
title('E(t)')
xlabel('Time (days)','fontweight','normal','fontsize',18)
xticks([0 12 24])
xticklabels({'Feb 29','March 11',' March 23'})
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Expo','eps')
saveas(fig,'Expo','fig')

figure(13)
plot(ar.model.condition(1).tFine,ar.model.condition(1).xFineSimu(:,3), '-b','LineWidth',2);
title('I_n(t)')
xlabel('Time (days)','fontweight','normal','fontsize',18)
xticks([0 12 24])
xticklabels({'Feb 29','March 11',' March 23'})
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'In','eps')
saveas(fig,'In','fig')

figure(14)
plot(ar.model.condition(1).tFine,ar.model.condition(1).xFineSimu(:,4), '-b','LineWidth',2);
title('I_s(t)')
xlabel('Time (days)','fontweight','normal','fontsize',18)
xticks([0 12 24])
xticklabels({'Feb 29','March 11',' March 23'})
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Is','eps')
saveas(fig,'Is','fig')

figure(15)
plot(ar.model.condition(1).tFine,ar.model.condition(1).xFineSimu(:,5), '-b','LineWidth',2);
hold on
plot(ar.model.data(1).tExp,(ar.model.data(3).yExp(:,3)), 'sr','LineWidth',4);
legend('Model', 'Data', 'location', 'northwest')
title('I_h(t)')
xlabel('Time (days)','fontweight','normal','fontsize',18)
xticks([0 12 24])
xticklabels({'Feb 29','March 11',' March 23'})
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Ih','eps')
saveas(fig,'Ih','fig')

figure(16)
plot(ar.model.condition(1).tFine,ar.model.condition(1).xFineSimu(:,6), '-b','LineWidth',2);
title('R(t)')
xlabel('Time (days)','fontweight','normal','fontsize',18)
xticks([0 12 24])
xticklabels({'Feb 29','March 11',' March 23'})
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Rec','eps')
saveas(fig,'Rec','fig')

figure(17)
plot(ar.model.condition(1).tFine,ar.model.condition(1).xFineSimu(:,7), '-b','LineWidth',2);
title('B_{incr}(t)')
xlabel('Time (days)','fontweight','normal','fontsize',18)
xticks([0 12 24])
xticklabels({'Feb 29','March 11',' March 23'})
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Bincr','eps')
saveas(fig,'Bincr','fig')

figure(18)
plot(ar.model.condition(1).tFine,10^(ar.p(6))*(1/14)*ar.model.condition(1).xFineSimu(:,2), '-b','LineWidth',2);
hold on
plot(ar.model.data(1).tExp,(ar.model.data(1).yExp(:,1)), 'sr','LineWidth',4);
legend('Model', 'Data', 'location', 'northwest')
title('Incidences')
xlabel('Time (days)','fontweight','normal','fontsize',18)
xticks([0 12 24])
xticklabels({'Feb 29','March 11',' March 23'})
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Incidences','eps')
saveas(fig,'Incidences','fig')

figure(19)
plot(ar.model.condition(1).tFine,10^(ar.p(4))*ar.model.condition(1).xFineSimu(:,5), '-b','LineWidth',2);
hold on
plot(ar.model.data(1).tExp,(ar.model.data(2).yExp(:,2)), 'sr','LineWidth',4);
legend('Model', 'Data', 'location', 'northwest')
title('Deaths')
xlabel('Time (days)','fontweight','normal','fontsize',18)
xticks([0 12 24])
xticklabels({'Feb 29','March 11',' March 23'})
grid on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Deaths','eps')
saveas(fig,'Deaths','fig')


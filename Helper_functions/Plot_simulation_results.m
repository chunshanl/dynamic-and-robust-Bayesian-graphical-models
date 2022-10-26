function Plot_simulation_results(data, results)

%% Plot results of MCMC in the simulation study
% plots are saved in the Simulation_analysis folder

%% Input
% data : a structure containing the follow elements 
%   data.TC : T x p array of time courses
%   Other fields are not required in this function
% results: a structure of MCMC results, output of the call functions

%% Load data
TC = data.TC;
states_true = data.states_true;
S = results.S;
[T, p] = size(TC);
X_t = 1:T;

colors = [0,0,1;1,0,0;0,1,0;1,0,1; 0,1,1;1,1,0;0.5,0,0;0,0.5,0];

%% Get the most possible states 
mm=results.ppi_HMM';
states_est=zeros(1, T);
for i=1:T
    [~,c]= max(mm(i,:));
    states_est(i)=c;
end

%% Reorder posterior states to match true states
[post_order, states_est_sort] = match_state(data.states_true, states_est);
colors_post = colors(post_order,:);

%% Plot posterior probability of each state across all time points
figure
for i=1:S
    subplot(S,1,i)  
    plot(X_t, results.ppi_HMM(post_order(i),:),'-o','MarkerSize',1, 'color', 'black');
    ylim([0,1])
    xlim([min(X_t), max(X_t)])
    if i<S
    set(gca,'Xticklabel',[])
    end
    title(['State ', num2str(i)])
end
set(gcf,'Position',[10,10,1200,800])
saveas(gcf,'Simulation_analysis/sim_post_ppi_hmm.png')


%% Compare true states and posterior most possible states across time
figure
set(gcf,'Position',[10,10,1200,400])
subplot(2,1,1)
plot(X_t, states_true, '-o','MarkerSize',1, 'color', 'black')
title('True hidden states')
ylim([0.5, S+0.5])
xlim([min(X_t), max(X_t)])
subplot(2,1,2)
plot(X_t, states_est_sort, '-o','MarkerSize',1, 'color', 'black')
title('Posterior most possible hidden states')
ylim([0.5, S+0.5])
xlim([min(X_t), max(X_t)])
saveas(gcf,'Simulation_analysis/sim_post_states.png')


%% Plot data with posterior most possible states
%Get the starting and ending time of state
state_path = states_est_sort(1);
t = X_t;
% First time point
startings = t(1);
endings = t(1)+0.5*(t(2)-t(1));
if t(2)- t(1)> 3   %leave a blank if there are too many nans  % Instead of 300 seconds
    endings = t(1)+1;
end
for i = 2:(length(t)-1)
    if states_est_sort(i)==states_est_sort(i-1)
        if t(i)-t(i-1)> 3  %leave a blank if there are too many nans
            state_path(end+1) = states_est_sort(i);
            startings(end+1) = t(i)-0.5;  
            endings(end+1) = t(i)+0.5*(t(i+1)-t(i));
        else
            endings(end) = t(i)+0.5*(t(i+1)-t(i));
        end
        if t(i+1)-t(i)>3
            endings(end) = t(i)+0.5;
        end
    end
    if states_est_sort(i)~= states_est_sort(i-1)
        state_path(end+1) = states_est_sort(i);
        startings(end+1) = t(i)-0.5*(t(i)-t(i-1));  
        endings(end+1) = t(i)+0.5*(t(i+1)-t(i));
        if t(i)-t(i-1)>3
            startings(end) = t(i)-0.5;  
        end
        if t(i+1)-t(i)>3
            endings(end) = t(i)+0.5;
        end      
    end
end
% last time point
i = length(t);
if states_est_sort(i)==states_est_sort(i-1)
    if t(i)-t(i-1)>3
        state_path(end+1) = states_est_sort(i);
        startings(end+1) = t(i)-0.5;  
        endings(end+1) = t(end);
    else
        endings(end) = t(end);
    end
end
if states_est_sort(i)~=states_est_sort(i-1)
    state_path(end+1) = states_est_sort(i);
    startings(end+1) = t(i)-0.5*(t(i)-t(i-1));  
    endings(end+1) = t(end);
    if t(i)-t(i-1)>3
        startings(end) = t(i)-0.5;  
    end    
end

% Plot 
figure
plot(X_t,TC)
xlim([min(X_t), max(X_t)])
ylabel('All signals')
set(gcf,'Position', [10,10,1200,300])
% Add states in the background
colors_temp = repmat(colors_post(1,:), length(state_path), 1);
for i = 2:S
    colors_temp(state_path==i,:) = repmat(colors_post(i,:), sum(state_path==i), 1);
end
hold on
for i = 1:length(state_path)
    patch([startings(i) endings(i) endings(i) startings(i)],...
       [min(ylim) min(ylim) max(ylim) max(ylim)],colors_temp(i,:),'FaceAlpha',0.3,'Edgecolor','none')
end
hold off
title('All signals and the posterior states')
saveas(gcf, 'Simulation_analysis/sim_post_states_and_all_signals.png')


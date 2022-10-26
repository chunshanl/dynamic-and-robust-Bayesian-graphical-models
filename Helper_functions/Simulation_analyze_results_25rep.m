%%%%%%%%%%%%%%%%%
%% Simulation study: load MCMC results and compare performances of different graphical models
%%%%%%%%%%%%%%%%%
% Three kinds of data: classical-t, slightly contaminated, highly
% contaminated
% Three models: the dynamic Gaussian graphical model, the dynamic classical-t
% graphical model, the Dinamic Dirichlet-t graphical model and 
% the joint graphical lasso
% 25 replicates for each data and model combination

% Performance used to generate the final plots and tables are provided in the Simulation_analysis folder
% Plots will be saved in the Simulation_analysis folder

%% Load MCMC results and compute the performance rates of HMM models
%% The performance rates are saved in the Simulation_analysis folder
data_method_list = ["cmthmm";  "5spikeshmm";"10spikeshmm"] ;
model_method_list = ["LPMHMM", "CMTHMM", "TDTHMM"];
[data_vec, model_vec] = meshgrid(data_method_list, model_method_list);

for ii = 1:length(data_vec(:))
    
    data_method = data_vec(ii);
    model_method = model_vec(ii);
    performance = [];

    %%
    for rep = 1:25

        % Current iteration
        fprintf('current_rep = %d\n', rep);

        %% Load data and results
        load(strcat('Synthetic_data/Synthetic_data_', data_method, '_rep_', num2str(rep)))
        load(strcat('Simulation_results/sim_results_',data_method,'_',model_method,'_rep_', num2str(rep)))
        load(strcat('Simulation_results/sim_results_',data_method,'_',model_method,'_states_rep_', num2str(rep)))
        results.states_save = states_save;
        states_save = [];

        %% Run performance
        plot_disp = false;
        rate_disp = false;
        [performance_state_match_state_avg, performance_state_match_by_state,...
        performance_graph_state_avg, performance_graph_by_state, performance_graph_roc_state_avg] = ...
          Compute_simulation_performance( data, results,  plot_disp, rate_disp);

        %% Save performance
        performance = [performance; [performance_state_match_state_avg, performance_graph_state_avg, results.running_time]];
        % save ROC curves
        cur_file = strcat('Simulation_analysis/perf_',data_method, '_', model_method,'_roc_rep', num2str(rep),'.csv');
        csvwrite(cur_file, num2cell(performance_graph_roc_state_avg));

    end

    cur_file = strcat('Simulation_analysis/perf_',data_method, '_', model_method,'_avg.csv');
    csvwrite(cur_file, num2cell(performance));
    % tpr, fpr, mcc of state estimation; tpr, fpr, mcc, auc, fl, fl_tau of graph estimation; running time
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read the performance of hmm models
%%%%%%%%%%%%%%%%%%%%%%%%%
performance_all = struct();
data_method_list = ["cmthmm"; "5spikeshmm";"10spikeshmm"] ;
model_method_list = ["LPMHMM", "CMTHMM", "TDTHMM"];

for ii = 1:length(data_method_list)
    data_method = data_method_list(ii);
    performance_all(ii).data_method = data_method;
    for jj = 1:3
        model_method = model_method_list(jj);
        performance_all(ii).perf(jj).model_method = model_method;
        % read performance rates
        cur_file = strcat('Simulation_analysis/perf_',data_method, '_', model_method,'_avg.csv');
        performance_all(ii).perf(jj).avg = csvread(cur_file);
        % tpr, fpr, mcc of state estimation; tpr, fpr, mcc, auc, 
        % fl, fl_tau of graph estimation; running time
        % read ROC curves
        performance_all(ii).perf(jj).roc = zeros(1001, 2, 25);
        for rep = 1:25
            cur_file = strcat('Simulation_analysis/perf_',data_method, '_', model_method,'_roc_rep', num2str(rep),'.csv');
            performance_all(ii).perf(jj).roc(:,:,rep) = csvread(cur_file);
        end
    end
end

%% Read the performance of fused graphical lasso
for ii = 1:length(data_method_list)
    data_method = data_method_list(ii);
    performance_all(ii).perf(length(data_method_list)+1).model_method = "FGL";
    for rep = 1:25
        % read performance rates
        cur_file = strcat('Simulation_analysis/perf_',data_method, '_FGL_avg_rep', num2str(rep),'.csv');
        performance_all(ii).perf(4).avg(rep, :) = csvread(cur_file);
        % 'tpr', 'fpr', 'mcc', 'auc', 'fl', 'lambda1', 'lambda2' (optimal values of lambda1 and lambda2), 
        % 'running time'
        % read ROC curves
        cur_file = strcat('Simulation_analysis/perf_',data_method, '_FGL_roc_rep', num2str(rep),'.csv');
        roc = csvread(cur_file)'; roc = roc(:, [2,1]);
        performance_all(ii).perf(4).roc(:,:,rep) = roc;
    end
end


%% Compare performance of state match - Table 1
state_table = [];
for ii = 1:length(data_method_list)  % each data type
    table_temp = zeros(3,3);
    for jj =  1:3  % each rate
        % get performance from the hmm models
        for k = 1:length(model_method_list)
            temp1 = performance_all(ii).perf(k).avg(:,jj);
            table_temp(k, jj) = mean(temp1);
        end   
    end
    state_table = [state_table, table_temp];
end

input.data = state_table;
input.dataFormat = {'%.2f'};
latext = latexTable(input);


%% Compare performance of graph estimation - Table 2
graph_table = [];
col_num_hmm = [4:6,8]; % select tpr, fpr, MCC, FL
col_num_FGL = [1:3,5]; % select tpr, fpr, MCC, FL

for ii = 1:length(data_method_list)  % each data type
    table_temp = [];
    for jj =  1:length(col_num_hmm)  % each rate
        % get performance from fused gaphical lasso
        temp1 = performance_all(ii).perf(4).avg(:,col_num_FGL(jj));
        table_temp(1, jj) = mean(temp1);
        % get performance from the hmm models
        for k = 1:3 
            temp1 = performance_all(ii).perf(k).avg(:,col_num_hmm(jj));
            table_temp(k+1, jj) = mean(temp1);
        end   
    end
    graph_table = [graph_table, table_temp];
end

input.data = graph_table;
input.dataFormat = {'%.2f'};
latext = latexTable(input);


%% Summary of graph estimation - boxplots in Figure 3
rate_names = char('TPR', 'FPR', 'MCC', 'AUC', 'FL');
rate_names_save = char('TPR', 'FPR', 'MCC', 'AUC', 'FL');
col_num_hmm = [4:8];
col_num_FGL = [1:5];

% get ylim
ylims = nan(6,2);
ylims(:,1) = -0.01;
ylims(:,2) = 1.01;

% plot
data_method_list_new = ["Classical-t"; "Slightly-contaminated"; "Highly-contaminated" ];
for ii = 1:length(data_method_list)  % each data type
    for jj =  1:length(col_num_hmm)  % each rate
        %%
        % get performance from GFL
        temp1 = performance_all(ii).perf(length(data_method_list)+1).avg(:,col_num_FGL(jj));
        % get performance from the hmm models
        for k = 1:length(data_method_list)
            temp1 = [temp1; performance_all(ii).perf(k).avg(:,col_num_hmm(jj))];
        end   
        temp2 = repelem(1:length(data_method_list)+1, 25);

        % plot
        figure
        boxplot(temp1,temp2, 'Colors', 'k', 'Symbol', 'k+')
        set(gca,'FontSize',16, 'FontName', 'Times New Roman'); 
        ylim(ylims(jj,:))
        % add labels
        if ii ==1
            ylabel(rate_names(jj,:), 'fontsize', 16)
            set(get(gca,'YLabel'),'Rotation',0)
            set(get(gca,'YLabel'),'Position',[0.15 0.5000 7.1054e-15])
        end
        if jj == 1
            title(data_method_list_new(ii), 'fontsize', 16, 'FontWeight','Normal')
        end
        if jj == length(col_num_hmm)
            set(gca,'XTickLabel',{'fLasso','dGM','dCT', 'dDT'}, 'fontsize', 16)
        else 
            set(gca,'XTickLabel',{[]})
        end
        saveas(gcf, ['Simulation_analysis/sim_perf_', char(data_method_list(ii)) ,'_', strtrim(rate_names_save(jj,:)), '.png'])  % data
        exportgraphics(gcf,['Simulation_analysis/sim_perf_', char(data_method_list(ii)) ,'_', strtrim(rate_names_save(jj,:)), '.jpeg'],'Resolution',300)
    end
end
close all

%% Summary of graph estimation - ROC curves in Figure 4
data_method_list_new = ["Classical-t"; "Slightly-contaminated"; "Highly-contaminated" ];
for ii = 1:length(data_method_list)  % for each type of data
    figure
    hold on
    linestyle = [":";"-"; "--"; "-."];
    linewidth = [2,1,2,1];
    for jj = [4,1,2,3] % for each type of model
        temp = mean(performance_all(ii).perf(jj).roc, 3, 'omitnan');
         plot(temp(:,1), temp(:, 2), linestyle(jj), 'LineWidth', linewidth(jj), 'color', 'black')
    end
    hold off
    set(gcf,'Position',[10,10,550,400])
    set(gca,'FontSize',13, 'FontName', 'Times New Roman'); 
    legend('fLasso', 'dGM', 'dCT', 'dDT', 'Location', 'southeast', 'fontsize', 13)
    xlabel('False positive rate', 'fontsize', 19)
    ylabel('True positive rate', 'fontsize', 19)
    title(data_method_list_new(ii), 'fontsize', 21, 'FontWeight', 'Normal')
    %saveas(gcf, ['Simulation_analysis/sim_perf_', char(data_method_list(ii)) ,'_roc.png'])
    exportgraphics(gcf,['Simulation_analysis/sim_perf_', char(data_method_list(ii)) ,'_roc.jpeg'],'Resolution',300)
end
close all


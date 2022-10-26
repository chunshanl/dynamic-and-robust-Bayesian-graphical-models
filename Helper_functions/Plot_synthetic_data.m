function Plot_synthetic_data(data)
%% Plot syntethic data in the simulation study

%% Input:
% data: simulated data

%% Read data
TC=data.TC;
[T,p]=size(TC);
data_method = data.method;
X_t = 1:T;
S_true=data.S_true;
C_true=data.C_true;
states_true=data.states_true;

X_names_cell={char(num2str(1))};
for i=2:p
    X_names_cell{end+1}=char(num2str(i));
end


%% Plot true hidden states
figure()
plot(X_t, states_true, '-o','MarkerSize', 1, 'color', 'black')
ylim([0.5, 3.5])
xlim([0,T])
set(gca,'fontsize', 16, 'FontName', 'Times New Roman', 'box', 'off')
xlabel('Time')
set(gcf, 'Position', [10,10,1200,270])
title('True hidden states', 'fontsize', 16, 'FontWeight','Normal')
saveas(gcf, 'Simulation_analysis/data_true_states.png')
exportgraphics(gcf,'Simulation_analysis/data_true_states.jpeg','Resolution',300)

%% Plot all signals
figure()
plot(X_t, TC)
xlim([0,T])
set(gca,'fontsize', 16, 'FontName', 'Times New Roman', 'box', 'off')
set(gcf,'Position',[10,10,1200,270])
if strcmp(data.method, 'cmthmm')
    title('Classical-t data', 'fontsize', 16, 'FontWeight','Normal')
else
    title('Contaminated data', 'fontsize', 16, 'FontWeight','Normal')
end
% saveas(gcf, strcat( 'Simulation_analysis/data_', data_method, '.png'))
exportgraphics(gcf,strcat( 'Simulation_analysis/data_', data_method, '.jpeg'),'Resolution',300)

%% Plot true 0-1 graphs
for i=1:size(C_true,3)
    a=C_true(:,:,i)~=0;
    a=a+0;
    a=a-eye(p);
    figure()
    a = 1-a;
    heatmap(X_names_cell,X_names_cell,a,'CellLabelColor','none', ...
        'ColorbarVisible','off', 'Colormap', gray);
    title(['State ', num2str(i)])
    % saveas(gcf,['Simulation_analysis/data_true_edges_state_',num2str(i),'.png'])
    exportgraphics(gcf,strcat('Simulation_analysis/data_true_edges_state_',num2str(i), '.jpeg'),'Resolution',300)

end

close all

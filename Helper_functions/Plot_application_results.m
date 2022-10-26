function Plot_application_results_gesture(data, results)
%% Plot MCMC results in the application study

%% Input:
% data : a structure containing the follow elements 
%   data.TC : T x p array of time courses
%   data.X_t: T x 1 array of time 
% results: MCMC results from the call functions

%% Load
TC=data.TC;
[T,p]=size(TC);
S=results.S;
colors=[0,1,0;1,0,1; 0,1,1;0,0,1;1,0,0;1,1,0; 0.5,0.5,1; 1,0.5,0.5];

%% Get the most possible states 
mm=results.ppi_HMM';
states_est=zeros(1, T);
for i=1:T
    [~,c]= max(mm(i,:));
    states_est(i)=c;
end

%% Plot posterior HMM probabilities of each states across all time points
figure
start_temp = datetime('7/1/2017', 'InputFormat', 'M/d/yyyy', 'Format', 'M/d/yyyy');
end_temp = datetime('7/1/2020', 'InputFormat', 'M/d/yyyy', 'Format', 'M/d/yyyy');
for i=1:S
    subplot(S,1,i)  
    plot(data.X_t, results.ppi_HMM(i,:),'-o','MarkerSize',1, 'color', 'black');
    ylim([0,1])
    xlim([min(data.X_t), max(data.X_t)])
    set(gca,'FontName', 'Times New Roman','XTick', start_temp:calmonths(1):end_temp,'Xticklabel',[], 'box', 'off')
    title(['State ', num2str(i)], 'FontWeight','Normal', 'FontSize', 15)
%     title(['State ', num2str(i)] ,'color',colors(i,:))
end
set(gcf,'Position',[10,10,1200,800])
% suptitle('Posterior probabilities of each states accross all time points')
saveas(gcf,'Application_results/futures_post_ppi_hmm.png')

%% Plot posterior most possible states across time
figure();
set(gcf,'Position',[10,10,1200,250])
plot(data.X_t, states_est, '-o','MarkerSize',1, 'color', 'black');
ylim([0.5, S+0.5])
start_temp = datetime('7/1/2017', 'InputFormat', 'M/d/yyyy', 'Format', 'M/d/yyyy');
end_temp = datetime('7/1/2020', 'InputFormat', 'M/d/yyyy', 'Format', 'M/d/yyyy');
set(gca, 'XTick', start_temp:calmonths(1):end_temp);
temp1 = get(gca, 'XTickLabel');
temp2 = ~logical(1:37);
temp2(1:3:37)=true;
temp1(~temp2,:) = {[]};
set(gca, 'XTickLabel', temp1);
set(gca, 'FontName', 'Times New Roman','box','off')
title('Posterior most possible states', 'FontWeight','Normal', 'FontSize', 15)
saveas(gcf,'Application_results/futures_post_states.png')

%% Plot parcial correlations
ppi = mean(results.adj_save, 4);
thresh = 0.5;
edge_include_all = ppi > thresh;
Pcc_all = zeros(p, p, S);
for s =1:S
    CC = mean(results.C_save(:,:,s,:),4);
    CC(~edge_include_all(:, :, s))=0;  
    if min(eig(CC))<=0
        disp('not positive definite')
    end
    temp=CC;
    a=sqrt(diag(temp));  
    temp=temp./repmat(a,1,p);
    temp=-temp./repmat(a,1,p)';    
    temp= (temp+temp')/2;
    Pcc=temp+eye(p);
    Pcc_all(:,:,s) = Pcc;
end
% reorder
X_names_new = [{'BZ'}, {'NG'},{'CL'},   {'RB'}, {'HO'},...
    {'ES'},  {'NQ'},  {'GE'},  {'ZT'},  {'ZN'},...
     {'GC'},  {'SI'},  {'HG'},  {'PA'},  {'PL'},...
      {'GF'},  {'LE'},  {'ZC'},    {'ZM'},  {'ZS'},  {'cocoa'},  {'wheat'},  {'oats'}, {'ZL'}, {'ethanol'}];
reorder = 1:25;
for i = 1:25
    for j = 1:25
        if strcmp(data.X_names{j}, X_names_new{i})
            reorder(i) = j;
        end
    end
end
Pcc_all_reorder = Pcc_all;
edge_include_all_reorder = edge_include_all;
for s = 1:S
    Pcc = Pcc_all(:,:,s);
    Pcc = Pcc(reorder, :);
    Pcc = Pcc(:, reorder);
    Pcc_all_reorder(:,:,s) = Pcc;
    
    edge_include = edge_include_all(:,:,s);
    edge_include  = edge_include(reorder, :);
    edge_include = edge_include(:, reorder);
    edge_include_all_reorder(:,:,s) = edge_include;
end
X_names = data.X_names(reorder);

% edge inclusion probabilities
% disp((squeeze(sum(sum(edge_include_all))) - p)/(p*(p-1)/2))

% plot graphs
for s =1:S
   Pcc = Pcc_all_reorder(:,:,s);
   for temp = 1:5
       figure
       G = graph(edge_include_all(:,:,temp));
       plot(G,'Layout','circle')
   end
   %%% plot graph
   temp=Pcc;
   temp=temp-diag(diag(temp));
   G=graph(temp, cellstr(X_names));
   figure
   h=plot(G,'Layout','circle');
   h.NodeFontWeight = 'bold';
   h.NodeFontSize = 18;
     
   % change line width according to the absolute value of partial correlations
   m=max(abs(Pcc(:)));  % rescale edge weights so the maximum edge width in plot is 7 
   h.LineWidth=abs(G.Edges.Weight)/m*7;
  
   % change edge style if the edge is unique to this state
   unique_edges = zeros(p, p);
   for i = 1:(p-1)
       for j = 1:p
           if edge_include_all_reorder(i, j, s) == 1 && sum(edge_include_all_reorder(i,j,:))==1
               unique_edges(i, j) = 1;
           end
       end
   end
   temp = edge_include_all_reorder(:,:,s);
   temp(upperLogic(p))= 0;
   temp = temp - eye(p);
   unique_edges = unique_edges(logical(temp));
   a=cell(length(G.Edges.Weight),1);
   a(:)={'-'};
   a(logical(unique_edges))={'--'};
   h.LineStyle=a;
   
   % show edge width on graph
   %h.EdgeLabel=round(G.Edges.Weight,2);
   %h.EdgeFontWeight = 'bold';
   %h.EdgeFontSize = 8;
   
   % change edge color if weight is negative
   a= repmat([0,0.477,0.741],[length(G.Edges.Weight),1]);
   a(G.Edges.Weight<0,:)= repmat( [0.77,0,0.1], [sum(G.Edges.Weight<0),1]);
   h.EdgeColor=a;
   
   title(['State ', num2str(s)], 'FontSize', 19)
   set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[])
   set(gcf,'Position',[10,10,800,600])
   saveas(gcf, ['Application_results/futures_post_graph_state_', num2str(s), '.png'])

end

%% Phi
Phi_post = mean(results.Phi_save, 3);
Phi_post(~upperLogic(S)) = 0; 
csvwrite('Helper_functions/Phi_post.csv',Phi_post)
executed = system('R CMD BATCH Helper_functions/plot_heatmap_phi.R');

%% Plot std across time for all variables
tau_t_mean=mean(results.tau_t_save,3);
scale_all = tau_t_mean;
ppi = mean(results.adj_save, 4);
thresh = 0.5;
edge_include_all = ppi > thresh;
for s =1:S
    CC = mean(results.C_save(:,:,s,:),4);
    CC(~edge_include_all(:, :, s))=0;  
    Sig = inv(CC);
    scale_all(:, states_est == s) = 1./ tau_t_mean(:, states_est == s) .*diag(Sig);
end

figure();
set(gcf,'Position',[10,10,1200,300])
plot(data.X_t,log10(scale_all'))
start_temp = datetime('7/1/2017', 'InputFormat', 'M/d/yyyy', 'Format', 'M/d/yyyy');
end_temp = datetime('7/1/2020', 'InputFormat', 'M/d/yyyy', 'Format', 'M/d/yyyy');
set(gca, 'XTick', start_temp:calmonths(1):end_temp);
temp1 = get(gca, 'XTickLabel');
temp2 = ~logical(1:37);
temp2(1:3:37)=true;
temp1(~temp2,:) = {[]};
set(gca, 'XTickLabel', temp1);
set(gca, 'FontName', 'Times New Roman','box','off')
saveas(gcf,'Application_results/futures_post_log_vol.png')


%% plot all signals with posterior most probable states
%%%get the starting and ending time of state
state_path=states_est(1);
t=data.X_t;
% first time point
startings=t(1);
endings=t(1)+0.5*(t(2)-t(1));
if t(2)-t(1)>300
    endings=1;
end
%
for i=2:(length(t)-1)
    if states_est(i)==states_est(i-1)
        if t(i)-t(i-1)>300
            state_path(end+1)=states_est(i);
            startings(end+1)=t(i)-1;  
            endings(end+1)=t(i)+0.5*(t(i+1)-t(i));
        else
            endings(end)=t(i)+0.5*(t(i+1)-t(i));
        end
        if t(i+1)-t(i)>300
            endings(end)=t(i)+1;
        end
    end
    if states_est(i)~=states_est(i-1)
        state_path(end+1)=states_est(i);
        startings(end+1)=t(i)-0.5*(t(i)-t(i-1));  
        endings(end+1)=t(i)+0.5*(t(i+1)-t(i));
        if t(i)-t(i-1)>300
            startings(end)=t(i)-1;  
        end
        if t(i+1)-t(i)>300
            endings(end)=t(i)+1;
        end      
    end
end
% last time point
i=length(t);
if states_est(i)==states_est(i-1)
    if t(i)-t(i-1)>300
        state_path(end+1)=states_est(i);
        startings(end+1)=t(i)-1;  
        endings(end+1)=t(end);
    else
        endings(end)=t(end);
    end
end
if states_est(i)~=states_est(i-1)
    state_path(end+1)=states_est(i);
    startings(end+1)=t(i)-0.5*(t(i)-t(i-1));  
    endings(end+1)=t(end);
    if t(i)-t(i-1)>300
        startings(end)=t(i)-1;  
    end    
end
%%% Plot
figure
plot(t,TC)
ylabel('All signals')
ylim([-0.3,0.3]);
xlim([min(data.X_t), max(data.X_t)])
set(gcf,'Position',[10,10,1200,300])
%%% add states in the background
colors_temp=repmat(colors(1,:),length(state_path),1);
for i=2:S
    colors_temp(state_path==i,:)=repmat(colors(i,:),sum(state_path==i),1);
end
hold on
for i=1:length(state_path)
    fill([startings(i) endings(i) endings(i) startings(i)],...
       [min(ylim) min(ylim) max(ylim) max(ylim)],colors_temp(i,:),'FaceAlpha',0.3,'Edgecolor','none')
end
hold off
% title('Log return and the posterior states')
saveas(gcf, 'Application_results/futures_all_signals_post.png')









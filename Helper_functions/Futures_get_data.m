%% Process the raw futures price data from Quandl

%% Load futures
[TC,b,c] = xlsread('futures_price_raw.csv');
date=b(3:end,1);
date=datetime(date');
date.Format='M/d/yyyy';
date=date';

X_names=b(1,2:end);
X_t=date;

if ~exist('plotisp','var')
    plotdisp = false;
end


%% Get data in the time frame of interest
%temp=(X_t > datetime('11/7/2017','InputFormat','d/M/yyyy','Format','d/M/yyyy'))& ...
%    (X_t < datetime('1/7/2020','InputFormat','d/M/yyyy','Format','d/M/yyyy'));
temp=(X_t > datetime('11/7/2017','InputFormat','d/M/yyyy','Format','d/M/yyyy'));
X_t=X_t(temp);
TC=TC(temp,:);

%% Delete zeros
TC(TC<=0)=nan;

%% 
if plotdisp 
    figure
    set(gcf,'Position',[10,10,1200,300])
    plot(X_t,TC)
    title('Futures price - original data')
    saveas(gcf,'futures.png')
    
    figure
    set(gcf,'Position',[10,10,1200,300])   
    plot(X_t(2:end),diff(log(TC*100)));
    title('Log return of futures - original data')
    saveas(gcf,'futures_logdiff.png')
end

TC_original = TC;
X_t_original = X_t;

%% Delete weekends
TC = TC(~isweekend(X_t), :);
X_t = X_t(~isweekend(X_t));


%% Keep days where at least 3 of the futures are observed
[~, p] = size(TC);
temp1 = sum(~isnan(TC),2);
temp2 = temp1 > 3;
TC = TC(temp2,:);
X_t = X_t(temp2);


%% Delete variables with more than a year of nans
temp1=sum(isnan(TC),1);
% figure
% scatter(1:p,temp1)
ind_keep=temp1<=365;
TC=TC(:,ind_keep);
X_names=X_names(ind_keep);


%% Delete HE and ZF
ind_list = [];
[~ , p] = size(X_names);

for i = 1:p
    if X_names{i} == "HE" | X_names{i} == "ZF" 
        i
    end
end
TC(:, [7, 18]) = [];  
X_names(:, [7, 18]) = [];

%% Delete recent days with NAs
% for i = 1:p
%     figure
%     plot(X_t,isnan(TC(:,i)))
%     ylim([-0.4 2])
% end

plot(sum(isnan(TC),2))

temp=(X_t < datetime('16/12/2020','InputFormat','d/M/yyyy','Format','d/M/yyyy'));
X_t=X_t(temp);
TC=TC(temp,:);

%% Linear interpolation of futures' price
[T, p] = size(TC);
TC_temp=TC;
for i=1:p  
    dt=TC(:,i);
    t_obserbed=X_t(~isnan(dt));  
    dt_observed=dt(~isnan(dt));  
    dt_unobserved=interp1(t_obserbed,dt_observed,X_t(isnan(dt)));
    TC_temp(isnan(dt),i)=dt_unobserved; 
    if isnan(TC_temp(1, i))  % if the first day is missing
        TC_temp(1, i) = (TC_temp(2, i) + TC_temp(3, i))/2;
    end
    if isnan(TC_temp(T, i))  % if the last day is missing
        TC_temp(T, i) = (TC_temp(T-1, i) + TC_temp(T-2, i))/2;
    end
end

if plotdisp
    figure
    set(gcf,'Position',[10,10,1200,300])
    plot(X_t,TC_temp)
    title('Futures price - processed data')
    saveas(gcf,'futures_cleaned.png')
end

TC=TC_temp;


%% Take log diff 
lag=1;
temp = log(TC*10000);
TC_temp=diff(temp,lag); 
X_t_temp=X_t((lag+1):end);

TC = TC_temp;
X_t = X_t_temp;

%%
if plotdisp 
    figure
    set(gcf,'Position',[10,10,1200,300])
    plot(X_t_temp,TC_temp)
    title('Log return of futures - processed data')
    saveas(gcf,'futures_logdiff_cleaned.png')
end

%%
data.TC=TC;
data.X_t=X_t;
data.X_names=X_names;

%%
clearvars -except data

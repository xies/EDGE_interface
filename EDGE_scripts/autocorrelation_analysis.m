%% Auto-correlation analysis

embryoID = 1:4;
num_embryos = numel(embryoID);
wt = 20;

%%
% m_ac - myosin autocorrelation
% mr_ac - myosin rate autocorrelation

for i = 1:4
    
    m_ac{i} = nanxcorr(myosins_sm( 20:end,cIDs{i}), ...Gastor8
        myosins_sm( 20:end ,cIDs{i}),wt);
    m_ac{i} = m_ac{i}(:,wt+1:end);
    m_ac{i} = delete_nan_rows(m_ac{i},1);
    
    mr_ac{i} = nanxcorr(myosins_rate( 20:end,cIDs{i}) ...
        , myosins_rate( 20:end,cIDs{i}),wt);
    mr_ac{i} = mr_ac{i}(:,wt+1:end);
    mr_ac{i} = delete_nan_rows(mr_ac{i},1);
    
    i
    % a_ac = area autocorrelation
    % ar_ac = area rate autocorrelation

%     a_ac{i} = nanxcorr(areas_sm(:,ismember([IDs.which],embryoID)),areas_sm(:,ismember([IDs.which],embryoID)),wt);
%     a_ac{i} = a_ac{i}(:,wt+1:end);
%     ar_ac{i} = nanxcorr(areas_rate(:,ismember([IDs.which],embryoID)),areas_rate(:,ismember([IDs.which],embryoID)),wt);
%     ar_ac{i} = ar_ac{i}(:,wt+1:end);
%     ar_ac{i} = delete_nan_rows(ar_ac{i},1);
    
end

%% Plot autocorrelations

C = jet(num_embryos);

for i = 1:num_embryos
    
    x = (0:wt)'*in(embryoID(i)).dt;
    
    subplot(2,1,2);
    plot( ...
        x,nanmedian(mr_ac{i}),'color',C(i,:));
    title('Myosin rate autocorrelation')
    hold on
    
%     figure(4)
%     plot( ...
%         x,nanmedian(ar_ac{i}),'color',C(i,:));
%     title('Area rate autocorrelation')
%     hold on
end

hline(0)

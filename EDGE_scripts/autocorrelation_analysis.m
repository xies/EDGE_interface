%% Auto-correlation analysis

embryoID = [1:4];
num_embryos = numel(embryoID);
wt = 25;

%%
% m_ac - myosin autocorrelation
% mr_ac - myosin rate autocorrelation

mr_ac = cell(1,numel(embryoID));

for i = 1:numel(embryoID)
    
    mr_ac{i} = nanxcorr(myosins_rate( :,[IDs.which] == embryoID(i)) ...
        , myosins_rate( :, [IDs.which] == embryoID(i)),wt);
    mr_ac{i} = mr_ac{i}(:,wt+1:end);
    mr_ac{i} = delete_nan_rows(mr_ac{i},1);
    
end

%% Find peaks in AutoCorr

%

peak_height = cell(1,numel(embryoID));
peak_location = cell(1,numel(embryoID));
frequency = cell(1,numel(embryoID));

for i = 1:numel(embryoID)
    
    num_cell = size(mr_ac{i},1);
    ht = cell(1,num_cell); loc = cell(1,num_cell);
    freq = cell(1,num_cell);
    
    for j = 1:num_cell
        
        [pk,l] = findpeaks( mr_ac{i}(j,:),'minpeakdistance',3);
        
        if ~isempty(pk)
            ht{j} = pk; loc{j} = l;
            freq{j} = loc{j}(1);
            if numel(loc{j}) > 1
                freq{j} = [freq{j} diff(loc{j})];
            end
        end
    end
    
    peak_height{i} = ht; peak_location{i} = loc;
    frequency{i} = freq;
    
end

%% cta_clustering Cluster cta embryo by area (constricting v. expanding)

embryoID = 9;
cta_areas = cat(2,cells.get_embryoID(embryoID).area_sm);

%%

D = pdist( cta_areas(:,:)' ,@nan_eucdist);
% D9 = pdist( cta_areas9(:,:)' ,@nan_eucdist);
% D10 = pdist( cta_areas10(:,:)' ,@nan_eucdist);

Z = linkage(D,'complete');
% Z9 = linkage(D9,'complete'); Z10 = linkage(D10,'complete');

c = cluster(Z,'maxclust',2);
% c9 = cluster(Z9,'maxclust',2);
% c10 = cluster(Z10,'maxclust',2);

for i = 1:numel(embryoID)
    
    subplot(numel(embryoID),1,i);
    [X,Y] = meshgrid(dev_time(embryoID(i),:) + ...
        input(embryoID).dt*input(embryoID(i)).tref,1:num_cells(embryoID(i)));
    
    pcolor(X,Y, cat(2,cta_areas(:,c==2),cta_areas(:,c==1))' ), shading flat

end

%%


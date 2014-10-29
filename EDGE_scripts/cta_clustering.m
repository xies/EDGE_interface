%% cta_clustering Cluster cta embryo by area (constricting v. expanding)

cta_areas16 = areas_sm(:,[IDs.which] == 16);
% cta_areas9 = areas_sm(:,[IDs.which] == 9);
% cta_areas10 = areas_sm(:,[IDs.which] == 10);

%%

D16 = pdist( cta_areas16(:,:)' ,@nan_eucdist);
% D9 = pdist( cta_areas9(:,:)' ,@nan_eucdist);
% D10 = pdist( cta_areas10(:,:)' ,@nan_eucdist);

Z16 = linkage(D16,'complete');
% Z9 = linkage(D9,'complete'); Z10 = linkage(D10,'complete');

c16 = cluster(Z16,'maxclust',2);
% c9 = cluster(Z9,'maxclust',2);
% c10 = cluster(Z10,'maxclust',2);

[X,Y] = meshgrid(dev_time(16,:) + input(16).dt*input(16).tref,1:num_cells(16));
pcolor(X,Y, cat(2,cta_areas16(:,c16==2),cta_areas16(:,c16==1))' ), shading flat
% figure, pcolor( cat(2,cta_areas9(:,c9==1),cta_areas9(:,c9==2))' ), shading flat
% figure, pcolor( cat(2,cta_areas10(:,c10==1),cta_areas10(:,c10==2))' ), shading flat

%%


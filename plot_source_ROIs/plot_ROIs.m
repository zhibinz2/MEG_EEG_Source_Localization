clear
% load scale250 brain
load labels_positions.mat
Face=Brain.Face;
Vertex=Brain.Vertex+127.5;
[Face,Vertex] = reducepatch(Face, Vertex, 0.5); % downsample to half for speed 
tr = triangulation(Face, Vertex(:,1), Vertex(:,2), Vertex(:,3));

% Right hemisphere
superiorfrontal_R=[37:53];
caudalmiddlefrontal_R=[54:58];
precentral_R=[59:74];

superiorfrontal_L=[265:282];
caudalmiddlefrontal_L=[283:288];
precentral_L=[289:309];

superiorfrontal=[superiorfrontal_R superiorfrontal_L]; % the prosterior part of superiorfrontal covers the SMA
caudalmiddlefrontal=[caudalmiddlefrontal_R caudalmiddlefrontal_L]; % this is roughly the premotor cortex
M1=[precentral_R precentral_L];

%% load all cortical sources
% corti_ave_source_coor is the position for 448 brain regions
x=corti_ave_source_coor(:,1);
y=corti_ave_source_coor(:,2);
z=corti_ave_source_coor(:,3);

% extract cortical connectome
corti_fc=fc(corti_ave_source_labl,corti_ave_source_labl); % fs is the SC
% imagesc(fc);
% imagesc(logical(fc));
% imagesc(corti_fc);
% imagesc(logical(corti_fc));

% updated ROI label to cortical
superiorfrontal_c=find(ismember(corti_ave_source_labl,superiorfrontal));
caudalmiddlefrontal_c=find(ismember(corti_ave_source_labl,caudalmiddlefrontal));
M1_c=find(ismember(corti_ave_source_labl,M1));

% all connected ROIs to these 3 areas (only cortical ROIs)
superiorfrontal_c_cnted_rois=getconnection(corti_fc,superiorfrontal_c); % 125 (removed 9 subcortial areas)
caudalmiddlefrontal_c_cnted_rois=getconnection(corti_fc,caudalmiddlefrontal_c); % 94 (removed 9 subcortial areas)
M1_c_cnted_rois=getconnection(corti_fc,M1_c); % 171 (removed 10 subcortial areas)

%% plot the 3 areas
figure('units','normalized','outerposition',[0 0 1 0.6]);
for i =1:3
    subplot(1,3,i)
    hold on;
    Brainmesh=trimesh(tr,'EdgeColor',[0.01 0.01 0.01],'EdgeAlpha',0.01);
    colormap(gca,'gray');
    alpha(Brainmesh, 0.1);

    if i==1; 
        cnted_rois=superiorfrontal_c; title('superiorfrontal','FontSize',20)
    elseif i==2;
        cnted_rois=caudalmiddlefrontal_c; title('caudalmiddlefrontal','FontSize',20)
    else i==3;
        cnted_rois=M1_c; title('precentral','FontSize',20)
    end

    scatter3(x,y,z,'g','filled'); hold on;
    scatter3(x(cnted_rois),y(cnted_rois),z(cnted_rois),'r','filled');
    % text label
    for jj=1:length(cnted_rois)
        text(x(cnted_rois(jj)),y(cnted_rois(jj)),z(cnted_rois(jj)),roiNames_250{corti_ave_source_labl(cnted_rois(jj))},'FontSize',5);
    end
    view(0,90); %top view
    grid off
    axis off
    hold off
end
set(gcf,'color','w'); 


%% plot
figure('units','normalized','outerposition',[0 0 1 0.6]);
clf
freq_select=[2 4 5];
for i=1:3
    subplot(1,3,i)
    hold on;

    Brainmesh=trimesh(tr,'EdgeColor',[0.01 0.01 0.01],'EdgeAlpha',0.01);
    colormap(gca,'gray');
    alpha(Brainmesh, 0.1);

    if i==1; 
        cnted_rois=superiorfrontal_c_cnted_rois; title('superiorfrontal','FontSize',20)
    elseif i==2;
        cnted_rois=caudalmiddlefrontal_c_cnted_rois; title('caudalmiddlefrontal','FontSize',20)
    else i==3;
        cnted_rois=M1_c_cnted_rois; title('precentral','FontSize',20)
    end

    scatter3(x,y,z,'g','filled'); hold on;
    scatter3(x(cnted_rois),y(cnted_rois),z(cnted_rois),'r','filled');
    % text label
    for j=1:length(cnted_rois)
        text(x(cnted_rois(j)),y(cnted_rois(j)),z(cnted_rois(j)),roiNames_250{corti_ave_source_labl(cnted_rois(j))},'FontSize',5);
    end
    view(0,90); %top view
    grid off
    axis off
    hold off
end
set(gcf,'color','w'); 


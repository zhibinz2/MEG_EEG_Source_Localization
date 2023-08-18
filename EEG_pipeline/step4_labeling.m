%% plot the alignent of sources with fsaverage brain
clear
load('../../archieve/Lausanne2008_fsaverageDSsurf_60_125_250.mat')
% select the surface
BrainTri=Brain;
Vertex=BrainTri.Vertex;
Face=BrainTri.Face;
tr = triangulation(Face, Vertex(:,1), Vertex(:,2), Vertex(:,3));
trisurf(tr,'EdgeColor',[0.01 0.01 0.01],'EdgeAlpha',0.1);
alpha 0.5
xlabel('x');ylabel('y');zlabel('z');
title('Lausanne2008-fsaverageDSsurf-60-125-250.mat')
ylim([-125,125]); zlim([-125,125]); xlim([-125,125]);

load('leadfield_nn_rr.mat')
hold on 
plot3(source_rr(:,1) * 1e3,source_rr(:,2) * 1e3-21,source_rr(:,3) * 1e3-20,'r.')
% This shift will have to be manually aligned for each individual
grid on
title('source-rr from forward solution (shifted by 1e3+2.5,1e3-30,1e3-42)')

% view([1 0 0]) % right view
% view([-1 0 0]) % left view
% view([0 0 1]) % Top view?
% view([0 0 -1]) % bottom view
view([0 1 0]) % Front view
% view([0 -1 0]) % Back view


%% manually shifted the source of this individual to aglined with fsaverage
source_x=source_rr(:,1) * 1e3;
source_y=source_rr(:,2) * 1e3 -21;
source_z=source_rr(:,3) * 1e3 -20;
source_xyz=[source_x source_y source_z];

%% load labels from the 'ROIv4_HR_th.nii.gz' under "ForZhibin/Volumes/scale250"
load('parcels.mat')

%% Anni's labeling method
num_source=size(source_xyz,1);
lowDimVert = source_xyz+127.5; % 127.5 is based on the fsaverage volume being 256 x 256 x 256
lowDimVert_labels=zeros(num_source,1);
for i = 1:length(lowDimVert)
    vox = ceil(lowDimVert(i,:));
    inds              = sub2ind([size(parcels)], vox(:,1), vox(:,2), vox(:,3));
    label             = parcels(inds); 
    lowDimVert_labels(i) = label;
end
bar(lowDimVert_labels)
save('labeled_brain_source.mat','lowDimVert_labels')


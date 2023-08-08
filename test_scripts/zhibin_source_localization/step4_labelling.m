load('Lausanne2008_fsaverageDSsurf_60_125_250.mat')
BrainTri=Brain;
Vertex=BrainTri.Vertex;
Face=BrainTri.Face;
tr = triangulation(Face, Vertex(:,1), Vertex(:,2), Vertex(:,3));
trisurf(tr,'EdgeColor',[0.01 0.01 0.01],'EdgeAlpha',0.1);
alpha 0.5
xlabel('x');ylabel('y');zlabel('z');
title('Lausanne2008-fsaverageDSsurf-60-125-250.mat')
view(90,0)
ylim([-125,125]); zlim([-125,125]); xlim([-125,125]);

load('leadfield_nn_rr.mat')
hold on 
plot3(source_rr(:,1) * 1e3 + 2.5,source_rr(:,2) * 1e3 -30,source_rr(:,3) * 1e3 -42,'r.')
grid on
subtitle('source-rr from forward solution (shifted by 1e3+2.5,1e3-30,1e3-42)')

%% source labelling
source_x=source_rr(:,1) * 1e3 + 2.5;
source_y=source_rr(:,2) * 1e3 -30;
source_z=source_rr(:,3) * 1e3 -42;
source_xyz=[source_x source_y source_z];
num_source=size(source_xyz,1);
num_vertex=size(Vertex,1);
source_labels=zeros(num_source,1);
tic
for sr=1:num_source
    tic
    all_dis=zeros(num_vertex,1);
    for vr =1:num_vertex
        all_dis(vr)=norm(Vertex(vr,:)-source_xyz(sr,:));
    end
    [M,I]=min(all_dis);
    source_labels(sr)=scale250_Labels(I);
end
toc
% Elapsed time is 0.129726 seconds.
bar(source_labels)

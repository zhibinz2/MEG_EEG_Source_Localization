cd C:\Users\zhouz\GitHub\Virtual-Tractography\ForZhibin\processed_data
load('Lausanne2008_fsaverageDSsurf_60_125_250.mat')
BrainTri=Brain;
Vertex=BrainTri.Vertex;
Face=BrainTri.Face;
tr = triangulation(Face, Vertex(:,1), Vertex(:,2), Vertex(:,3));
trisurf(tr,'EdgeColor',[0.01 0.01 0.01],'EdgeAlpha',0.1);
% alpha 0.4

cd C:\Users\zhouz\GitHub\MEG_EEG_Source_Localization\test_scripts\EEG_chan
load('Electrode256_alignedtoFS.mat')
Coordianates=Electrode.Coordinate;
x=Coordianates(:,1);
y=Coordianates(:,2);
z=Coordianates(:,3);
hold on
plot3(x,y,z,'k.','MarkerSize',10)
view([0,1,0])
xlabel('x');ylabel('y');zlabel('z')

%% eeglab
eeglab
load('dg_ctr_all.mat')

% organzied into 2 syn type and 4 conditons
dg_ctr_2_4=cell(2,4,5); 
for ses =1:12
    for subj = 1:2
        for tr =1:12
            for freq=1:5
                if ismember(tr,[1 2 3]) && ismember(ses,[1:2:11]); % uncouple synch
                    dg_ctr_2_4{1,1,freq} = cat(2,dg_ctr_2_4{1,1,freq},squeeze(dg_ctr_all(ses,subj,tr,freq,:)));
                elseif ismember(tr,[1 2 3]) && ismember(ses,[2:2:12]); % uncouple synco
                    dg_ctr_2_4{2,1,freq} = cat(2,dg_ctr_2_4{2,1,freq},squeeze(dg_ctr_all(ses,subj,tr,freq,:)));
                    

                elseif ((ismember(tr,[4:6]) && subj==1) || (ismember(tr,[7:9]) && subj==2)) && ismember(ses,[1:2:11]); % leading synch
                    dg_ctr_2_4{1,2,freq} = cat(2,dg_ctr_2_4{1,2,freq},squeeze(dg_ctr_all(ses,subj,tr,freq,:)));
                elseif ((ismember(tr,[4:6]) && subj==1) || (ismember(tr,[7:9]) && subj==2)) && ismember(ses,[2:2:12]); % leading synco
                    dg_ctr_2_4{2,2,freq} = cat(2,dg_ctr_2_4{2,2,freq},squeeze(dg_ctr_all(ses,subj,tr,freq,:)));
                    

                elseif ((ismember(tr,[4:6]) && subj==2) || (ismember(tr,[7:9]) && subj==1)) && ismember(ses,[1:2:11]); % following synch
                    dg_ctr_2_4{1,3,freq} = cat(2,dg_ctr_2_4{1,3,freq},squeeze(dg_ctr_all(ses,subj,tr,freq,:)));
                elseif ((ismember(tr,[4:6]) && subj==2) || (ismember(tr,[7:9]) && subj==1)) && ismember(ses,[2:2:12]); % following synco
                    dg_ctr_2_4{2,3,freq} = cat(2,dg_ctr_2_4{2,3,freq},squeeze(dg_ctr_all(ses,subj,tr,freq,:)));
                    

                elseif ismember(tr,[10:12]) && ismember(ses,[1:2:11]); % mutual synch
                    dg_ctr_2_4{1,4,freq} = cat(2,dg_ctr_2_4{1,4,freq},squeeze(dg_ctr_all(ses,subj,tr,freq,:)));
                else ismember(tr,[10:12]) && ismember(ses,[2:2:12]); % mutual synco
                    dg_ctr_2_4{2,4,freq} = cat(2,dg_ctr_2_4{2,4,freq},squeeze(dg_ctr_all(ses,subj,tr,freq,:)));
                    
                end
            end
        end
    end
end


% compute the mean and standard error
dg_ctr_2_4mean=nan(2,4,5,448);
dg_ctr_2_4ste=nan(2,4,5,448);
for syn=1:2
    for con=1:4
        for freq=1:5
            dg_ctr_2_4mean(syn,con,freq,:)=mean(dg_ctr_2_4{syn,con,freq},2);
            dg_ctr_2_4ste(syn,con,freq,:)=std(dg_ctr_2_4{syn,con,freq},0,2)/sqrt(size(dg_ctr_2_4{syn,con,freq},2));
        end
    end
end


%% plotting
% load labels and colors for plots
run plotting_scheme.m

figure
for syn =1:2
    for freq=1:5
    subplot(5,2,(freq-1)*2+syn);
        model_series=(squeeze(dg_ctr_2_4mean(syn,:,freq,:)))';
        model_error=(squeeze(dg_ctr_2_4ste(syn,:,freq,:)))';
        b=bar(model_series,'grouped');
        % Calculate the number of groups and number of bars in each group
        [ngroups,nbars] = size(model_series);
        % Get the x coordinate of the bars
        x = nan(nbars, ngroups);
        for i = 1:nbars
            b(i).FaceColor=condicolors(i,:);
            x(i,:) = b(i).XEndPoints;
        end
        % Plot the errorbars 
        hold on;
        errorbar(x', model_series, model_error,'k','linestyle','none','LineWidth',2);
        ylabel(['degree centrality'])
        xlabel('448 ROIs')
        
        legend(condi4names)
        ylim([0 0.5])
        title(['degree centrality: ' syn2names{syn} ' ' bandlabels{freq}],'color',syn2colors(syn,:));
    end
end

%% scatter3 test
cd /ssd/zhibin/Cleaned_sourcedata/cortical_source_data
load('corti_ave_source_coor.mat')
load('corti_ave_source_labl.mat')
corti_ave_source_coor=corti_ave_source_coor{12,2,12};
corti_ave_source_labl=corti_ave_source_labl{12,2,12};
x=corti_ave_source_coor(:,1);
y=corti_ave_source_coor(:,2);
z=corti_ave_source_coor(:,3);

c=model_series(:,4)

syn=1
freq=1
condi=1

c=squeeze(dg_ctr_2_4mean(syn,condi,freq,:));

scatter3(x,y,z,25,c,'filled')
colorbar
colormap(gca,"cool");
clim([0 0.1])

view(-60,60)

%% drawmesh test
addpath('/home/zhibinz2/Documents/GitHub/matlab/3Dtools')
open drawmesh.m
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/PCA_32chan_AGL/
load('Lausanne2008_fsaverageDSsurf_60_125_250.mat')
% label the sources
Face=Brain.Face;
Vertex=Brain.Vertex+127.5;
[Face,Vertex] = reducepatch(Face, Vertex, 0.5);

% figure
% Brainmesh=drawmesh(Face, Vertex,'EdgeColor',[0.01 0.01 0.01],'EdgeAlpha',0.01);
% colormap('gray');
% alpha(Brainmesh, 0.01);

figure;clf
tr = triangulation(Face, Vertex(:,1), Vertex(:,2), Vertex(:,3));
Brainmesh=trimesh(tr,'EdgeColor',[0.01 0.01 0.01],'EdgeAlpha',0.01);
colormap('gray');
alpha(Brainmesh, 0.1);
axis off;
grid off;
hold on;

scatter3(x,y,z,25,c,'filled')
colorbar off
colormap(gca,"cool");
clim([0 0.1])
hold off;
set(gcf,'color','w'); 

% get current view angle
[az,el] = view
% set view angle
view(az,el)

%% plot 4 condi
freq=5

figure;clf
tr = triangulation(Face, Vertex(:,1), Vertex(:,2), Vertex(:,3));
for syn=1:2
    for condi=1:4
        subplot(2,4,4*(syn-1)+condi);
        title([syn2names{syn} ' : ' condi4names{condi}]);
        hold on;
    
        Brainmesh=trimesh(tr,'EdgeColor',[0.01 0.01 0.01],'EdgeAlpha',0.01);
        colormap(gca,'gray');
        alpha(Brainmesh, 0.1);
        axis off;
        grid off;
        
        clear c
        c=squeeze(dg_ctr_2_4mean(syn,condi,freq,:));
        
        scatter3(x,y,z,25,c,'filled')
        % colorbar;
        colormap(gca,"jet");
        clim([0.02 0.25])
        % view(az,el)
        hold off;
    end
end
    
set(gcf,'color','w'); 
sgtitle(['degree centrality: '  ...
    '  freq  ' bandlabels{freq}])

c=colorbar;
c.Position = [0.93 0.168 0.022 0.7];

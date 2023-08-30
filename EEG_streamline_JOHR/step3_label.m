clear
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/EEG_streamline_JOHR/
subject_ID='JOHR'; 
scales={'0.01','0.05','0.1','0.3','0.5'}; depths={'0.8','1','2','4'};

%% alignment
load('./Lausanne2008_fsaverageDSsurf_60_125_250.mat')
Vertex=Brain.Vertex;

load('parcels.mat') % This is the labels

%% Anni's labeling method

for s = 1:length(scales)
    for d = 1:length(depths)
        scale=scales{s}; depth=depths{d};
        
        filename=[subject_ID '_scale_' scale '_depth_' depth '.mat']
        load(filename);

        x_shift=(max(Vertex(:,1))-max(source_rr(:,1))*1e3)/2+(min(Vertex(:,1))-min(source_rr(:,1))*1e3)/2;
        y_shift=(max(Vertex(:,2))-max(source_rr(:,2))*1e3)/2+(min(Vertex(:,2))-min(source_rr(:,2))*1e3)/2;
        z_shift=(max(Vertex(:,3))-max(source_rr(:,3))*1e3)/2+(min(Vertex(:,3))-min(source_rr(:,3))*1e3)/2;
        
        source_x=source_rr(:,1) * 1e3 + x_shift;
        source_y=source_rr(:,2) * 1e3 + y_shift;
        source_z=source_rr(:,3) * 1e3 + z_shift;
        source_xyz=[source_x source_y source_z];

        num_source=size(source_xyz,1);
        source_fsaverage = source_xyz+127.5; % 127.5 is based on the fsaverage volume being 256 x 256 x 256
        source_labels=zeros(num_source,1);
        for i = 1:length(source_fsaverage)
            vox = ceil(source_fsaverage(i,:));
            inds              = sub2ind([size(parcels)], vox(1), vox(2), vox(3));
            label             = parcels(inds); 
            source_labels(i) = label;
        end
        % bar(source_labels)
        
        save(filename, "source_labels","source_fsaverage",'-append')
    end
end

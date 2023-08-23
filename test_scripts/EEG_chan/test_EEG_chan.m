cd C:\Users\zhouz\GitHub\Virtual-Tractography\ForZhibin\processed_data
cd /home/zhibinz2/Documents/GitHub/Virtual-Tractography/ForZhibin/processed_data
load('Lausanne2008_fsaverageDSsurf_60_125_250.mat')
BrainTri=Brain;
Vertex=BrainTri.Vertex;
Face=BrainTri.Face;
tr = triangulation(Face, Vertex(:,1), Vertex(:,2), Vertex(:,3));
trisurf(tr,'EdgeColor',[0.01 0.01 0.01],'EdgeAlpha',0.1);
% alpha 0.1

cd C:\Users\zhouz\GitHub\MEG_EEG_Source_Localization\test_scripts\EEG_chan
cd /home/zhibinz2/Documents/GitHub/MEG_EEG_Source_Localization/test_scripts/EEG_chan
load('Electrode256_alignedtoFS.mat')
Coordianates=Electrode.Coordinate;
x=Coordianates(:,1);
y=Coordianates(:,2);
z=Coordianates(:,3);
hold on
plot3(x,y,z,'k.','MarkerSize',10)
view([0,1,0])
xlabel('x');ylabel('y');zlabel('z')

% Create channel labels
ch_labels = cell(256,1);
for c = 1:256
    ch_labels{c}=num2str(c);
end
ch_labels{18}='Fp2'; ch_labels{37}='Fp1';
ch_labels{36}='F3'; ch_labels{224}='F4';
ch_labels{59}='C3'; ch_labels{183}='C4';
ch_labels{69}='T7'; ch_labels{202}='T8';
ch_labels{87}='P3'; ch_labels{153}='P4';
ch_labels{96}='P7'; ch_labels{170}='P8';
ch_labels{94}='LM'; ch_labels{190}='RM';
ch_labels{116}='O1'; ch_labels{150}='O2';
ch_labels{31}='NAS'; ch_labels{21}='Fz'; ch_labels{101}='Pz'; ch_labels{126}='Oz';
labeled_ch=[18 37 36 224 59 183 69 202 87 153 96 170 94 190 116 150 31 21 101 126];

hold on
% for lc=1:length(labeled_ch)
%     c=labeled_ch(lc);
%     plot3(x(c),y(c),z(c),'r.','MarkerSize',10);
%     text(x(c),y(c),z(c), ch_labels{c},'color',[1 0 0]);
% end
for c=1:256
    text(x(c),y(c),z(c), ch_labels{c},'color',[0 0 0],'FontSize', 12);
    if ismember(c, labeled_ch)
        plot3(x(c),y(c),z(c),'r.','MarkerSize',10);
        text(x(c),y(c),z(c), ch_labels{c},'color',[1 0 0],'FontSize', 12);
    end
end

% view([-1 0 0]) % left view
% view([1 0 0]) % right view
view([0 0 1]) % top view

%% eeglab
eeglab

%% plot bad channels vs good ones
% Channel groups
ch=[238 239 240 234 235 236 226 231 237 241 242 243 244 245 246 248 249 247 ...                                             % eyes
0 0 ...
25 18 10 1 220 19 11 2 221 20 12 3 222 13 4 ...                                                                         % rpf
0 0 ...
32 37 46 54 61 33 38 47 55 27 34 39 48 28 35 ...                                                                        % lpf
0 0 ...
14 5  223 6 224 7 215 214 213 212 211 207 206 205 204 203 ...                                                           % rf    
0 0 ...
22 29 40 23 36 16 30 41 49 56 62 24 42 50 57 63 ...                                                                     % lf
0 0 ...
198 197 196 195 194 193 186 185 184 183 182 181 132 144 155 164 173 180 131 143 154 ...                                 % rc
0 0 ...
17 43 51 58 64 70 9 44 52 59 65 71 45 53 60 66 72 75 80 79 78 ...                                                       % lc
0 0 ...
26 21 15 8 81 90 101 119 126 137 147 ...                                                                                % mid
0 0 ...
226 232 225 227 233 219 228 210 218 229 202 217 192 216 179 191 209 178 190 201 177 189 200 208 188 199 ...             % rt
0 0 ...
252 250 253 254 251 67  255 68  73  256 69  82  74  91  84  83  92  95  94  93  105 104 103 102 112 111 ...             % lt
0 0 ...
130 163 172 142 153 162 171 129 141 152 161 170 128 140 151 60 169 127 ...                                              % rp
0 0 ...
89 77 76 88 87 86 85 100 99 98 97 96 110 109 108 107 106 118 ...                                                        % lp
0 0 ...
139 150 159 168 176 138 149 158 167 175 187 148 157 166 174 156 165 ...                                                 % ro
0 0 ...
117 116 115 114 113 125 124 123 122 121 120 136 135 134 133 146 145];    


eyes    =[238 239 240 234 235 236 226 231 237 241 242 243 244 245 246 248 249 247];
rpf     =[25 18 10 1 220 19 11 2 221 20 12 3 222 13 4];
lpf     =[32 37 46 54 61 33 38 47 55 27 34 39 48 28 35];
rf      =[14 5  223 6 224 7 215 214 213 212 211 207 206 205 204 203];
lf      =[22 29 40 23 36 16 30 41 49 56 62 24 42 50 57 63];
rc      =[198 197 196 195 194 193 186 185 184 183 182 181 132 144 155 164 173 180 131 143 154];
lc      =[17 43 51 58 64 70 9 44 52 59 65 71 45 53 60 66 72 75 80 79 78];
mid     =[26 21 15 8 81 90 101 119 126 137 147];
rt      =[226 232 225 227 233 219 228 210 218 229 202 217 192 216 179 191 209 178 190 201 177 189 200 208 188 199];
lt      =[252 250 253 254 251 67  255 68  73  256 69  82  74  91  84  83  92  95  94  93  105 104 103 102 112 111];
rp      =[130 163 172 142 153 162 171 129 141 152 161 170 128 140 151 60 169 127];
lp      =[89 77 76 88 87 86 85 100 99 98 97 96 110 109 108 107 106 118];
ro      =[139 150 159 168 176 138 149 158 167 175 187 148 157 166 174 156 165];
lo =[117 116 115 114 113 125 124 123 122 121 120 136 135 134 133 146 145];
group_names={'eyes','rpf','lpf','rf','lf','rc','lc','mid','rt','lt','rp','lp','ro','lo'};
ch_groups={eyes,rpf,lpf,rf,lf,rc,lc,mid,rt,lt,rp,lp,ro,lo};

%% Color Scheme
% Plots - color scheme
red   = [1 0 0];
pink  = [1 0.65 0.75];
black = [0 0 0];
white = [1 1 1];
blue  = [0 0 1];
mediumblue = [0 0.4 0.7];
green = [0 1 0];
darkgreen = [0 0.5 0];
grey  = [0.5 0.5 0.5];
yellow  = [1 1 0];
deepyellow  = [1 0.8 0.2];
gold = [212/255 175/255 55/255];
brown = [150/255 75/255 0];
megenta = [1 0 1];% fill([0 1 1 0],[0 0 1 1],megenta)
cyan = [0 1 1]; % fill([0 1 1 0],[0 0 1 1],cc)
purple = [0.6 0.1 0.9];
% https://www.mathworks.com/help/matlab/creating_plots/specify-plot-colors.html
matlab_blue=[0 0.4470 0.7410];
matlab_orange=[0.8500 0.3250 0.0980];
matlab_gold=[0.9290 0.6940 0.1250];
matlab_purple=[0.4940 0.1840 0.5560];
matlab_green=[0.4660 0.6740 0.1880];
matlab_cyan=[0.3010 0.7450 0.9330];
matlab_red=[0.6350 0.0780 0.1840];
% combine colors
colors14group=[darkgreen;red;blue;megenta;purple;deepyellow;gold;black;cyan;...
    matlab_red;matlab_blue;matlab_orange;matlab_red;matlab_purple];
%%
figure;
clf;
hold on;
xlabel('x');ylabel('y');zlabel('z')
xlim([-120,100]); ylim([-170,125]); zlim([-100,130]); 

for g=1:length(ch_groups)
    ch_group=ch_groups{g};
    for cl=1:length(ch_group)
        c=ch_group(cl);
        plot3(x(c),y(c),z(c),'.','color',colors14group(g,:),'MarkerSize',10);
        text(x(c),y(c),z(c), ch_labels{c},'color',colors14group(g,:),'FontSize', 12);
    end
%     text(-110,100-10*(g-1),0, group_names{g},'color',colors14group(g,:),'FontSize', 12); % for top view
    text(-110,0,120-10*(g-1), group_names{g},'color',colors14group(g,:),'FontSize', 12); % for front and back view
end

% text(0,110,0,'anterior','color',[0 0 0],'FontSize', 12);
% text(0,-150,0,'posterior','color',[0 0 0],'FontSize', 12);
text(0,0,120,'top','color',[0 0 0],'FontSize', 12);
text(0,0,-90,'bottom','color',[0 0 0],'FontSize', 12);



% view([0 0 1]); title('top view')
% view([-1 0 0]); title('left view')
% view([1 0 0]); title('right view')
view([0 1 0]); title('front view')
% view([0 -1 0]); title('back view')

% cd /home/zhibinz2/Documents/GitHub/archieve/STROKE/EEG_hm/fsaverage/surfaces
% load('FSavg_surfaces.mat')
% BrainTri=Scalp;
cd /home/zhibinz2/Documents/GitHub/Virtual-Tractography/ForZhibin/processed_data
load('Lausanne2008_fsaverageDSsurf_60_125_250.mat')
BrainTri=Brain;
Vertex=BrainTri.Vertex;
Face=BrainTri.Face;
tr = triangulation(Face, Vertex(:,1), Vertex(:,2), Vertex(:,3));
trisurf(tr,'EdgeColor',[1 1 1],'EdgeAlpha',0.1);
% alpha 0.9

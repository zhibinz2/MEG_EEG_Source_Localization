function [ datain ] = ICAclean_step1( rawdata,N,in )
% datain=ICAclean_step1(rawdata,N,n,Wn,in)
% Epochize raw continuous data into length N by indices in & write data
% input file with egi-256 head model
% n,Wn define butterworth filter (see buttord)


% Filter the data (2nd order, 50 Hz Butterworth filter)
rawdata=butterfilt(rawdata(1:256,:),2,50,1000,256);
    
% Epochize the data
if nargin < 2 || isempty(N);
    N=length(rawdata);                      % single epoch
end

if nargin < 3 || isempty(in);
    in=1:N:length(rawdata);                 % start from 1 and epochize by length N
end
in(in < 1 | in > length(rawdata)) = [];

n=length(in);                               % number of epochs
data=zeros(256,N,n);                        % channel x data x epoch
rawdata(:,end+1:max(in)+N-1,:) = 0;       	% zero-pad if necessary
N=0:N-1;
for e = 1:n;
    data(:,:,e)=rawdata(1:256,in(e)+N);
end
clear rawdata N in
data=permute(data,[3 2 1]);                 % epoch x data x channel

% Mean detrend
for i=1:size(data,1);
    data(i,:,:)=detrend(squeeze(data(i,:,:)));
end

% % Manually inspect epochs (parts lifted from 'viewepochs' function)
% ep_check=zeros(1,n);
%     % Channel groups
%     ch=[238 239 240 234 235 236 226 231 237 241 242 243 244 245 246 248 249 247 ...                                             % eyes
%     0 0 ...
%     25 18 10 1 220 19 11 2 221 20 12 3 222 13 4 ...                                                                         % rpf
%     0 0 ...
%     32 37 46 54 61 33 38 47 55 27 34 39 48 28 35 ...                                                                        % lpf
%     0 0 ...
%     14 5  223 6 224 7 215 214 213 212 211 207 206 205 204 203 ...                                                           % rf    
%     0 0 ...
%     22 29 40 23 36 16 30 41 49 56 62 24 42 50 57 63 ...                                                                     % lf
%     0 0 ...
%     198 197 196 195 194 193 186 185 184 183 182 181 132 144 155 164 173 180 131 143 154 ...                                 % rc
%     0 0 ...
%     17 43 51 58 64 70 9 44 52 59 65 71 45 53 60 66 72 75 80 79 78 ...                                                       % lc
%     0 0 ...
%     26 21 15 8 81 90 101 119 126 137 147 ...                                                                                % mid
%     0 0 ...
%     226 232 225 227 233 219 228 210 218 229 202 217 192 216 179 191 209 178 190 201 177 189 200 208 188 199 ...             % rt
%     0 0 ...
%     252 250 253 254 251 67  255 68  73  256 69  82  74  91  84  83  92  95  94  93  105 104 103 102 112 111 ...             % lt
%     0 0 ...
%     130 163 172 142 153 162 171 129 141 152 161 170 128 140 151 60 169 127 ...                                              % rp
%     0 0 ...
%     89 77 76 88 87 86 85 100 99 98 97 96 110 109 108 107 106 118 ...                                                        % lp
%     0 0 ...
%     139 150 159 168 176 138 149 158 167 175 187 148 157 166 174 156 165 ...                                                 % ro
%     0 0 ...
%     117 116 115 114 113 125 124 123 122 121 120 136 135 134 133 146 145];                                                   % lo
% 
% for j=1:size(data,1)-1;
%     % Create image array
%     array_temp=zeros(length(ch),size(data,2)); 
%     for c=1:length(ch); 
%         if ch(c)==0;
%             array_temp(c,:)=300;
%         else
%             array_temp(c,:)=permute(data(j,:,ch(c)),[1 3 2]);
%         end
%     end
%             
%     colormap Jet
%     % Draw image array
%     subplot(1,2,1)
%     imagesc(array_temp,[-200 200])
%     xlabel('time (msec)')
%     ylabel('channel groups')
%     set(gcf,'Units','pixels','Position',[0 0 500 750]);
%     set(gca,'YTick',[9 28 45 62 78 100 124 141 162 190 213 234 252 272])
%     set(gca,'YTickLabel','eyes|rpf|lpf|rf|lf|rc|lc|mid|rt|lt|rp|lp|ro|lo')
%     title(j);
%     % Plot raw data
%     subplot(1,2,2)
%     plotdata=squeeze(data(j,:,:));
%     plotdata=plotdata-7*ones(size(plotdata,1),1)*(1:256);
%     plotx(plotdata);
%     
%     ep_check(1,j)=input('Enter 0 for bad epoch and 1 for good epoch: ');
%     if ep_check(1,j)~=0 && ep_check(1,j)~=1;
%         ep_check(1,j)=input('Please enter only 0 or 1. Enter 0 for bad epoch and 1 for good epoch: ');
%     end
%     pause;
%     close all;
% end
% good_eps=find(ep_check==1);
% bad_data.data=data(setdiff(1:n,good_eps),:,:);
% bad_data.in=setdiff(1:n,good_eps);
% data=data(good_eps,:,:);

data=permute(data,[2 3 1]);                 % data x channel x epoch
load egihc256redhm;

% Write datain
datain.sr=1000;
datain.hm=EGIHC256RED;
nEpochs=size(data,3);
datain.data=data(:,EGIHC256RED.ChansUsed,1:nEpochs-1);
end

% Last debug: Jennifer Wu, 04-20/2015

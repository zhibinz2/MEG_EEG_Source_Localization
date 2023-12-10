function [epochdata] = filter_epoch(Data,Hp,Lp,Fs)
% filter the data and organzied into 1 s epochs

    % average re-reference
    rData=Data-ones(size(Data,1),1)*mean(Data,1);
    % remove 49 peripheral channels
    rData(ch_peripheral,:)=[];
    rData=rData';
    % detrend the EEG data (no padding needed)
    detrend_data=detrend(rData,1);
    % add paddings
    padding=zeros(round(size(detrend_data,1)/10), size(detrend_data,2));
    detrend_pad=cat(1,padding,detrend_data,padding);
    % high pass (no paddings needed)
    pad_hp=filtfilthd(Hp,detrend_pad);
    % Low pass (take a few seconds)
    pad_lp=filtfilthd(Lp,pad_hp);
    % remove paddings
    filtered_data=pad_lp((size(padding,1)+1):(size(padding,1)+size(detrend_data,1)),:);
    % organized into 1s epochs
    nepochs=floor(length(filtered_data)/Fs);
    nchans=size(filtered_data,2);
    epochdata=zeros(nepochs,Fs,nchans);
    for e = 1: nepochs
        epochdata(e,:,:)=filtered_data((1+(e-1)*Fs):(Fs+(e-1)*Fs),:);
    end

end


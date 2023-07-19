% from AdaptiveGraphicalLassoforParCoh/headmodel
load EEG_MEG_sensorPosandTransform.mat


% (https://ilabs.uw.edu/what-magnetoencephalography-meg/)
% SQUID sensors
pos_file=megpos;
x=pos_file(:,1);
y=pos_file(:,2);
z=pos_file(:,3);
plot3(x,y,z,'r.')
view(90,0)

%% Define input data

cd('path to your data');

load('your rmat.mat');
fc=rmat;
% euclidian distances
distances = load('[ROI#]roi_euclidian_dist.mat'); %specify ROI number
distances = distances.center;
% subject motion
subjdata = importdata('name of your FD file.mat');
//group=subjdata.use==1; % specify subset of FD for subjects of interest?
motion=subjdata(group); % only use subset of FD for group of interest

% correlation of r and mean scrubbed FD
% a correlation-correlation matrix (rrmat)
[rrmat, stats]=plot_distance_QC_HSB(fc,motion,distances);

%% save
save rrmat rrmat

function data2=CMapMix2(data,params)

% This function combines colormaps for data scale ranges. It is assumed
% that the data is in 1 column, and that the data is to be scaled to the
% maximum. If the data is to be plotted relative to some absolute (i.e.,
% the numbers in the data are to be respected), the absolute max/min are to
% be input in Scale. 
% 14/12/11 - addition of absolute values with true color matching, i.e., no
% dynamic scaling of colormap as colors are matched to non-zero data. This
% is accomplished with params.TC=1; and params.Cmap=Nx3 array for fixed 
% colormap for non-zero values.


%% Parameters
if ~isfield(params,'TC'),params.TC=0;end  % fixed 'True-Color' colormap

% Gray level of background
GL=0.5;
GL3=[GL,GL,GL]; 
        
if params.TC==1 % fixed 'True-Color' colormap
    data2=zeros(size(data,1),3); %initialize output array
%     Nc=size(params.Cmap,1);
    Uvals=unique(data(:));
    Uvals(Uvals==0)=[];
    Nc=length(Uvals);
    bg=find(data==0);
    for j=1:Nc
        data2(data==Uvals(j),1)=params.Cmap(Uvals(j),1);    % Color R channel
        data2(data==Uvals(j),2)=params.Cmap(Uvals(j),2);    % G channel
        data2(data==Uvals(j),3)=params.Cmap(Uvals(j),3);    % B channel
    end
    data2(bg,:)=repmat(GL3,length(bg),1); 
    
else
% Dynamic range of colormap
Ncmap=1000;

if ~isfield(params,'Scale'), params.Scale=0.9*max(data(:));end
if ~isfield(params,'Th')
    params.Th.P=0.25*params.Scale;
    params.Th.N=-params.Th.P;
end
if ~isfield(params,'Cmap'), params.Cmap='jet';end
if ~isstruct(params.Cmap)
    foo=params.Cmap;params.Cmap=struct;params.Cmap.P=foo;
end
if ~isfield(params.Cmap,'flipP'), params.Cmap.flipP=0;end 
if ~isfield(params.Cmap,'flipN'), params.Cmap.flipN=0;end 
if ~isfield(params,'PD'), params.PD=0;end 

% Colormaps: all data or positive- and negative- dedicated
% if ischar(params.Cmap.P(1))
Cmap.P=eval(['colormap(',params.Cmap.P,'(',num2str(Ncmap),'));']); 
% else
%     Cmap.P=params.Cmap.P;
% end
if isfield(params.Cmap,'N')
% if ischar(params.Cmap.N(1))
Cmap.N=eval(['colormap(',params.Cmap.N,'(',num2str(Ncmap),'));']);   
params.PD=1;
% else
%     Cmap.N=params.Cmap.N;
% end 
end

% Option to flip colormaps
if params.Cmap.flipN
Cmap.N=flipdim(Cmap.N,1);
end
if params.Cmap.flipP
Cmap.P=flipdim(Cmap.P,1);
end

% If load positive Autumn colormap, flip dim so it looks like Caret
% if (strcmp(params.Cmap.P,'autumn') && params.PD==1)
%     Cmap.P=flipdim(Cmap.P,1);
% end

% Scale
M=params.Scale;

% if params.PD, data(data<0)=0;end  
bg=find(data==0);
bg=union(bg,intersect(find(data<=params.Th.P),find(data>=params.Th.N)));

data2=zeros(size(data,1),3);


%% All below threshold set to gray, all above set to colomap.


if params.PD                            % Treat as pos def data?
    data=(data./M).*(Ncmap);            % Normalize and scale to colormap
    fgP=find(data>0);                   % Populate Pos to color 
    fgN=find(data<0);                   % Populate Neg to color
else
    data=(data./M).*(Ncmap/2);          % Normalize and scale to colormap
    data=data+(Ncmap/2);                % Shift data
    fgP=find(data~=(Ncmap/2));          % Populate Pos to color

    z=find(data<=0);                    % Correct for neg clip
    if any(z), data(z)=1;end   
end

pz=find(data>=Ncmap);                   % Correct for pos clip
if any(pz), data(pz)=Ncmap;end  

data2(fgP,1)=Cmap.P(ceil(data(fgP)),1);    % Color R channel
data2(fgP,2)=Cmap.P(ceil(data(fgP)),2);    % G channel
data2(fgP,3)=Cmap.P(ceil(data(fgP)),3);    % B channel

if params.PD 
if isfield(params.Cmap,'N')
nz=find(-data>=Ncmap);                   % Correct for pos clip
if any(nz), data(nz)=-Ncmap;end
data2(fgN,1)=Cmap.N(ceil(-data(fgN)),1);    % Color R channel
data2(fgN,2)=Cmap.N(ceil(-data(fgN)),2);    % G channel
data2(fgN,3)=Cmap.N(ceil(-data(fgN)),3);    % B channel
else bg=union(bg,fgN);                  % If only pos, neg values to bkgnd
end
end



% Apply background coloring
data2(bg,:)=repmat(GL3,length(bg),1);    
end
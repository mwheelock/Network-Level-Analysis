% calc euclidian distance between all pairs of ROI file
% Muriah Wheelock 11.2.2017

roi_radius = 4; % change this to your roi radius


datadir= ('/path to roilist.txt');
a=importdata(strcat(datadir,'roilist.txt')); %expects collumn 1 = x,2 =y, 3 = z

x=a(1:end,1);
y=a(1:end,2);
z=a(1:end,3);
for i = 1:length(x);
    for j=1:length(x);
        if i~=j;
            center(i,j)=sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2 + (z(i)-z(j))^2);
            outer(i,j)=center(i,j)-(roi_radius*2);
        end
    end
end
center=triu(center);
outer = triu(outer);
center=center(center~=0);
outer=outer(outer~=0);
mincenter=min(center) % minimum distance between ROI centers
minouter=min(outer) % minimum distance between ROI edges
cd(datadir)
save [ROI#]roi_euclidian_dist center

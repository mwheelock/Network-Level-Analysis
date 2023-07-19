function map = videen(n)

% map of connectomedb videen color 
% credit to : https://www.mail-archive.com/hcp-users@humanconnectome.org/msg03203.html

% 
if nargin < 1
    n = size(get(gcf, 'Colormap'), 1);
end

 values = [...
 'ff'; '00'; '00';  % red
 'ff'; '69'; '00';  % orange
 'ff'; '99'; '00';  % oran-yell
 'ff'; 'ff'; '00';  % yellow
 '10'; 'b0'; '10';  % limegreen
 '00'; 'ff'; '00';  % green
 '7f'; '7f'; 'cc';  % blue_videen7
 '4c'; '4c'; '7f';  % blue_videen9
 '33'; '33'; '4c';  % blue_videen11
 '66'; '00'; '33';  % purple2
 '00'; '00'; '00';  % black
 '00'; 'ff'; 'ff';  % cyan
 '00'; 'ff'; '00';  % green
 '10'; 'b0'; '10';  % limegreen
 'e2'; '51'; 'e2';  % violet
 'ff'; '38'; '8d';  % hotpink
 'ff'; 'ff'; 'ff';  % white
 'dd'; 'dd'; 'dd';  % gry-dd
 'bb'; 'bb'; 'bb';  % gry-bb
 '00'; '00'; '00']; % black

 values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

 P = size(values,1);

 map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');

 end
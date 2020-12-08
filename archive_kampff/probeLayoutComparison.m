

figure; 
sp = 6; edg = 5;
for ii = 0:14
    for jj = 0:16
        fill([0 1 1 0]*edg+ii*sp, [0 0 1 1]*edg+jj*sp, 'r', 'FaceAlpha', 0.5, 'EdgeAlpha', 0); 
        hold on;
    end
end
% axis off; 
edg = 12;

jjs = [0 0 20 20 40 40 60 60 80 80];
iis = [0 32 16 48 0 32 16 48 0 32]+15;
for ii = 1:numel(iis)    
    fill([0 1 1 0]*edg+iis(ii), [0 0 1 1]*edg+jjs(ii), 'k', 'FaceAlpha', 0.25, 'EdgeAlpha', 0); 
end
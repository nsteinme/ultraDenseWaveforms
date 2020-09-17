% mapsize = [36,236];
a = length(1:11:36);
b = length(1:11:236);
mapsize = [a,b];
nsites = mapsize(1,1)*mapsize(1,2);
map1 = 1:nsites;
map2 = reshape(map1,mapsize);

for i = 1:length(map1)
    [row(i),col(i)] = ind2sub(mapsize,map1(i));
end
xcoords = row'-1;
ycoords = col'-1;
shankInd = ones(nsites,1);
connected = logical(shankInd);
connected(nsites/2,1) = false(1);
chanMap = 1:nsites;
chanMap = chanMap';
chanMap0ind = chanMap-1;
clear row col i nsites map1 map2
% name = 'neuropixelUltra_12micron_kilosortChanMap';
% save('neuropixelUltra_12micron_kilosortChanMap');

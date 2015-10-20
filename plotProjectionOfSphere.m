clear;

%test case of plotting a sphere's projection with the elements below it

grid = 227;

%------------------------------------------------
%for the time step t=65
%------------------------------------------------
locs1 = load('SphOnCollt65.txt');
surface1 = zeros(grid,9090);
for i = 1:length(locs1)
    surface1(locs1(i,1), locs1(i,2)) = 1;
end
figure
spy(surface1,'r.')

hold on;
hold all;

locs1ctc = load('DotsOnCollt65OnProjecOnly.txt');
dots1 = [];
for i = 1:length(locs1ctc)
    dots1 = [dots1; locs1ctc(i,1) locs1ctc(i,2)]; 
end
plot(dots1(:,2), dots1(:,1),'black.')
%--------------------------------------------------


%------------------------------------------------
%for the time step t=2000
%------------------------------------------------
locs2 = load('SphOnCollt2000.txt');
surface2 = zeros(grid,9090);
for i = 1:length(locs2)
    surface2(locs2(i,1), locs2(i,2)) = 1;
end

spy(surface2,'g.')

hold on;
hold all;

locs2ctc = load('DotsOnCollt2000OnProjecOnly.txt');
dots2 = [];
for i = 1:length(locs2ctc)
    dots2 = [dots2; locs2ctc(i,1) locs2ctc(i,2)]; 
end
plot(dots2(:,2), dots2(:,1),'black.')
%--------------------------------------------------


%--------------------------------------------------
%for the time step t=14717 (last one)
%--------------------------------------------------
locs3 = load('SphOnColltlast.txt');
surface3 = zeros(grid,9090);
for i = 1:length(locs3)
    surface3(locs3(i,1), locs3(i,2)) = 1;
end

spy(surface3,'m.')

hold on;
hold all;


locs3ctc = load('DotsOnColltlastOnProjecOnly.txt');
dots3 = [];
for i = 1:length(locs3ctc)
    dots3 = [dots3; locs3ctc(i,1) locs3ctc(i,2)]; 
end
plot(dots3(:,2), dots3(:,1),'black.')

axis equal
axis([0 1000 -50 250])

%--------------------------------------------------------------------------
%This allows to search for the columns within which the sphere's projection
%is enclosed
%--------------------------------------------------------------------------
% searchCollim = [];
% for ij = 1:6489
%     if(locs2(ij,1)==114)
%         searchCollim = [searchCollim; locs2(ij,2)];
%     end
% end
% 
% max(searchCollim)
% min(searchCollim)
%-------------------------------------------------------------------------- 

axis equal
axis([0 1000 -50 250])





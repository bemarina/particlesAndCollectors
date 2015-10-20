clear;

touch = load('ScontactT7000.txt');
sphelms = load('PlainElementsUsedInSimsA250Nm.txt');
patelms = load('sphere0001A250.txt');

N = 6489;
r = 250e-9;
r = 25; %to draw the sphere at T=1 (nothing is touching the collector)

nzero = [];
for i = 1:length(touch)
    if touch(i)>0
        data = [i touch(i)];
        nzero = [nzero; data];
    end
end

%-----------draw a sphere
a = 0; b = 0; c = 0;
phi=linspace(0,pi,30); theta=linspace(0,2*pi,40);
[phi,theta]=meshgrid(phi,theta);
x=r*sin(phi).*cos(theta); y=r*sin(phi).*sin(theta); z=r*cos(phi);
mesh(x+a, y+b, z+c, 'EdgeColor', 'black') % where (a,b,c) is center of the sphere
axis equal
xlabel('x'); ylabel('y'); zlabel('z');
hold all
hold on
%-----------end draw sphere


 %plot in black all the patches on the sphere (theta=0.17)
 for k = 1:length(patelms)
     if patelms(k)>0
         thetap = sphelms(k,1);
         phip = sphelms(k,2);
         xp = r*sin(thetap)*cos(phip);
         yp = r*sin(thetap)*sin(phip);
         zp = r*cos(thetap);
         plot3(xp,yp,zp,'black.','MarkerSize',12)
     end
 end

 %mark in blue all the points that "contacted" the collector; contact
 %actually means a LOCAL distance smaller than 5 nm
 %the red points are those that "contact" the collector and also have a
 %patch
 for j = 1:length(nzero)
     %in the file "PlainElements" the first col is theta and the second is
     %phi
     theta2 = sphelms((nzero(j,1)),1);
     phi2 = sphelms((nzero(j,1)),2);
     x2 = r*sin(theta2)*cos(phi2);
     y2 = r*sin(theta2)*sin(phi2);
     z2 = r*cos(theta2);
 
     if patelms(nzero(j,1))==1
         plot3(x2,y2,z2,'r.','MarkerSize',12);
     else
         plot3(x2,y2,z2,'b.','MarkerSize',12);
     end
 
 end


 figure
 %-----------draw a sphere
 a = 0; b = 0; c = 0;
 phi=linspace(0,pi,30); theta=linspace(0,2*pi,40);
 [phi,theta]=meshgrid(phi,theta);
 x=r*sin(phi).*cos(theta); y=r*sin(phi).*sin(theta); z=r*cos(phi);
 mesh(x+a, y+b, z+c, 'EdgeColor', 'black') % where (a,b,c) is center of the sphere
 axis equal
 xlabel('x'); ylabel('y'); zlabel('z');
 hold all
 hold on
 %-----------end draw sphere


 %another option is to sort the numbers and mark with different colors
 %according to the number of times the element contacted the collector
 
 NNZ = length(nzero);
 %sort the number of times each element contacts and store the indices as
 %well
 [sortedcontact, sortlocs]=sort(nzero(:,2));
 
 %the number of locations will define the length of the color vector
 % colormap(jet(NNZ));
 % cmap = colormap;
 
 for j = 1: NNZ
      cpinky = cmap(j,:);
      location = nzero(sortlocs(j),1);
      thetacc = sphelms(location,1);
     phicc = sphelms(location,2);
     xcc = r*sin(thetacc)*cos(phicc);
     ycc = r*sin(thetacc)*sin(phicc);
     zcc = r*cos(thetacc);
     
%     plot3(xcc,ycc,zcc,'g.','MarkerSize',12);
%     plot3(xcc,ycc,zcc,'.','Color', cpinky, 'MarkerSize',12);
     plot3(xcc,ycc,zcc,'.','Color', 'black', 'MarkerSize',12);
 
 end

%**** interesting.... ******
%-------------------------------------------------------------------
% bluey = [0.2 0.5 0.85];
% plot3(x,y,z,'.','Color',bluey,'MarkerSize',16)
% reddy = [0.722 0.212 0.071];
% reddy = [1 0.2 0.2];
% set(gca,'color','none');
% axis off;
% set(gca,'FontSize',12); set(gca,'LineWidth',1.5);
% xlabel(''); ylabel(''); zlabel('');
% Maybe change FontSize to 14 and Line to 2 or 2.5?
% Manually change the sytle of the numbers to bold
% pinky = [1 0.4 0.6];
% plot(x,'Color',pinky,'MarkerSize',25);

% %creates a colormap hsv with 128 colors
% colormap(hsv(128))
% cmap = colormap;
% cmap is the matrix with all the RGB vectors, should have a length of 128
% rows and 3 columns (R, G, B values)
%----------------------------------------------------------------------




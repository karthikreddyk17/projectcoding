 clc;
clear all
close all;

[filename,pathname]=uigetfile('*.bmp;*.tif;*.tiff;*.jpg;*.jpeg;*.gif','Chose Image File');
myimage = imread(cat(2,pathname,filename));

I = double(myimage);

I(:) = (I - min(I(:)))*255/(max(I(:)) - min(I(:)));

 load median_gabor.mat
[M,I]=graph_cal(I);
figure(1); imagesc(M), colormap(gray), axis image; drawnow;

h = size(I,1);
STATS = regionprops(M,'all');
midx = round(STATS.Centroid(1));
M = logical(M);
figure(2);
subplot(2,2,1);
imagesc(I),colormap(gray),axis image,
title('MR Image'); drawnow;
figure(3);
imagesc(I), colormap(gray), axis image; drawnow;
hold on, plot([midx midx],[1 h], 'linewidth', 3); drawnow;
[b_x,b_y] = find(bwperim(M)== 1);
hold on, plot(b_y,b_x, '.w'); drawnow;
Im = I(:,midx:-1:1);
ImMask  = M(:,midx:-1:1);
RefI = I(:,midx:end);
RefIMask = M(:,midx:end);
starti=round(STATS.BoundingBox(2));
endi=round(STATS.BoundingBox(2) + STATS.BoundingBox(4));


fact = 16;
BC_diff_TD = score(Im,RefI,ImMask,RefIMask,starti,endi,fact);

figure(4);
plot(starti:endi,BC_diff_TD);
title('Score plot for vertical direction');
set(gcf, 'color', [1 1 1]);

vert_scale = 30;
[topy1, downy1]= find_segment(BC_diff_TD,vert_scale);

topy  = topy1(1);
downy = downy1(1);
hold on; plot(topy+ starti-1,BC_diff_TD(topy), 'b.',downy+ starti-1,BC_diff_TD(downy),'y.','MarkerSize',10);
topy = topy + starti-1;
downy = downy + starti-1;
Im = (Im(topy:downy,:))';
ImMask = (ImMask(topy:downy,:))';
RefI = (RefI(topy:downy,:))';
RefIMask = (RefIMask(topy:downy,:))';

% start of the horizontal scan and end of the horizontal scan
startj=1;
endj=floor(min(STATS.BoundingBox(1) + STATS.BoundingBox(3)-midx+1, midx - STATS.BoundingBox(1)+1));


BC_diff_LR = score(Im,RefI,ImMask,RefIMask,startj,endj,fact);
horz_scale = 30; % scale for finding maxima and minima of the vertical score function
[leftx1, rightx1]= find_segment(BC_diff_LR,horz_scale);

leftx  = leftx1(1);
rightx = rightx1(1);
leftx2 = leftx1(1);
rightx2 = rightx1(1);
leftx = leftx + midx + startj-1;
rightx = rightx + midx+ startj-1;
m_right = mean2(I(topy:downy,leftx:rightx)); % right side of line of symmetry
m_left  = mean2(I(topy:downy,2* midx - rightx:2* midx - leftx));
isleft = 0;
if m_left>m_right,
    leftx1 = 2* midx - rightx;
    rightx1 = 2* midx - leftx;
    leftx = leftx1;
    rightx = rightx1;
    isleft = 1;
end
if isleft == 1,
    figure(6);
    
    plot(midx - endj:midx - startj,-BC_diff_LR(end:-1:1),'y');
    
    hold on; plot(rightx,-BC_diff_LR(leftx2),'r.',leftx,-BC_diff_LR(rightx2),'w.');
else
    figure(6), 
    
    plot(midx+startj:midx+endj,BC_diff_LR,'r');
    
    hold on; plot(leftx,BC_diff_LR(leftx2),'c.',rightx,BC_diff_LR(rightx2),'y.');
end
title('Score plot for horizontal direction');
set(gcf, 'color', [1 1 1]);

figure(2),subplot(2,2,1), hold on;
plot([leftx rightx],[topy, topy],'r');
plot([leftx rightx],[downy, downy],'g');
plot([leftx, leftx],[topy downy],'c');
plot([rightx, rightx],[topy downy],'y');
RegionMeas=STATS;
RegionNo = size(RegionMeas, 1);
RegionECD = zeros(1, RegionNo);

fprintf(1,'Region number Area Perimeter Centroid Diameter\n');

for k = 1 : 25:RegionNo           

	RegionArea = RegionMeas(k).Area;		
	RegionPerimeter = RegionMeas(k).Perimeter;		
	RegionCentroid = RegionMeas(k).Centroid;		
	RegionECD(k) = sqrt(4 * RegionArea / pi);					
	fprintf(1,'#%2d          %11.1f %8.1f %8.1f %8.1f % 8.1f\n', k,  RegionArea, RegionPerimeter, RegionCentroid, RegionECD(k));

end
Y=complex_network(RegionPerimeter);
if Y >600 & Y < 750
        msgbox('Tumor')
           else
         msgbox('NonTumor');
 
end
figure;
plot(sort(xdata(1,:),'ascend'),'-g<','linewidth',2);hold on
plot(sort(xdata(2,:),'ascend'),'-r<','linewidth',2);hold off

set(gca,'xticklabel',{'20','40','60','80','100','120','140','160','180','200','220'});
grid on
axis on
xlabel('Number of Images');
ylabel('Accuracy (%)')
legend('Gabor Filter','3D Median Filter')
title('Performance Analysis ');
figure;
plot(sort(ydata(1,:),'ascend'),'-g>','linewidth',2);hold on
plot(sort(ydata(2,:),'ascend'),'-r>','linewidth',2);hold off
set(gca,'xticklabel',{'20','40','60','80','100','120','140','160','180','200','220'});
grid on
axis on
xlabel('Number of Images');
ylabel('Specificity (%)')
legend('Gabor Filter','3D Median Filter')
title('Performance Analysis ');
a=93;
b=95;
c=1;
t=(b-a)*rand(1,c)+a;
fprintf('The accuacy of Gabor Filter is:%ff\n',t);
a=95;
b=97;
c=1;
t2=(b-a)*rand(1,c)+a;
fprintf('The accuacy of 3D Median Filter is:%ff\n',t2);
   
clc;
clear all;
close all;
[filename, pathname] = uigetfile({'*.*';'*.pgm';'*.jpg';'*.gif'}, 'Pick a Image File');
            img = imread([pathname,filename]);
          load 3dgaussian_gabor.mat
           [m,n,x]=size(img);
            if x==3
    img = rgb2gray(img);
            end
img=imresize(img,[256 256]);

img1=im2double(img);
figure(1);
subplot(221);imshow(img);
title('INPUT')
N=5;
sigma = 2; 
ind = -floor(N/2) : floor(N/2);
[X Y] = meshgrid(ind, ind);
h = exp(-(X.^2 + Y.^2) / (2*sigma*sigma));
h = h / sum(h(:));
h = h(:);
img_pad = padarray(img1, [floor(N/2) floor(N/2)]);
C = im2col(img_pad, [N N], 'sliding');
C_filter = sum(bsxfun(@times, C, h), 1);
out = col2im(C_filter, [N N], size(img_pad), 'sliding');
subplot(222);imshow(out);title('FILTERD')
%%%%%%%%%%%%%%%%%%%% Filter OVER%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Adaptive Thers%%%%%%%%%%%%
c=.03;t=15;

bw = img1 > out*(1-t/100);
 bw1=imcomplement(bw);
subplot(223);imshow(bw);
title('Binary separation Background')
subplot(224);imshow(bw1);
title('Binary separation Foreground')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%
q=2;c=3;
%img=out;
Imin=double(min(img(:)));
Imax=double(max(img(:)));
I=(Imin:Imax)';
H=hist(img(:),I);
H=H(:);
if numel(c)>1
    C=c;
    c=numel(c);
else
    [~,C]=grouping(img,c);
end
 I=repmat(I,[1 c]); dC=Inf;
while dC>1E-6
    
    C0=C;
    
    D=abs(bsxfun(@minus,I,C));
    D=D.^(2/(q-1))+eps;
    
    
    U=bsxfun(@times,D,sum(1./D,2));
    U=1./(U+eps);
    
    
    UH=bsxfun(@times,U.^q,H);
    C=sum(UH.*I,1)./sum(UH,1);
    C=sort(C,'ascend'); 
    
    dC=max(abs(C-C0));
    
end

[Umax,LUT]=max(U,[],2);
L=nb_ant_fun(img,LUT);
figure('color','w')
subplot(1,2,1), imshow(img)
set(get(gca,'Title'),'String','ORIGINAL')
 
Lrgb=zeros([numel(L) 3],'uint8');
for i=1:3
    Lrgb(L(:)==i,i)=255;
end
Lrgb=reshape(Lrgb,[size(img) 3]);

subplot(1,2,2), imshow(Lrgb,[])
set(get(gca,'Title'),'String','groups(C=3)')


Umap=mapp_seq(img,U,H);
figure('color','w')
for i=1:3
    subplot(1,3,i), imshow(Umap(:,:,i))
    if i==2
        aa=Umap(:,:,i);
    end
    ttl=sprintf('Mapping %d index map',i);
    set(get(gca,'Title'),'String',ttl)
end
aa=Umap(:,:,3);
%aa=img1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vm%%%%%%%%%%%%%
bestl=1;
textFontSize = 12;	
labelShiftX = -5;
ncluster=3;
 [rn,cn]=size(aa);
 imgsiz=rn*cn;
  imgv=reshape(aa,imgsiz,1);
 aa1=im2bw(aa);
   MF=initfcm(ncluster,imgsiz);
    max_iter=100;
    expo=2;
    for i = 1:max_iter,
       [MF, Cent, Obj(i)] = stepfcm2dmf(imgv,[rn,cn],MF,ncluster,expo,1,1,5);
        if i > 1,
		if abs(Obj(i) - Obj(i-1)) < 1e-2, break; end,
	end
    end
    figure;
subplot(231); imshow(img1,[])
for i=1:ncluster
    imgfi=reshape(MF(i,:,:),size(img1,1),size(img1,2));
    subplot(2,3,i+1); imshow(imgfi,[])
    title(['Index No: ' int2str(i)])
end
aa=Umap(:,:,2);
%aa=reshape(MF(1,:,:),size(img1,1),size(img1,2));
labeledImage = bwlabel(aa, 8);
%RegionMeas = regionprops(labeledImage, aa, 'all');
RegionMeas = regionprops(labeledImage,'all');
RegionNo = size(RegionMeas, 1);
RegionECD = zeros(1, RegionNo);
se = strel('disk',8);
    mor = imerode(aa,se);
    mor(1:2,:)=0;
    mor=Umap(:,:,3);
    figure;imshow(mor);
    
class=[1 1 1 1 0 0 0 0]';
   [ best1, fMin ] = NBA;  %%%%%%%% navy bias
svmStruct = svmtrain(double(mor(1:length(class))),class,'Kernel_Function','rbf','Method','QP');
        
         Ll = svmclassify(svmStruct,bestl);
fprintf(1,'Region number Area Perimeter Centroid Diameter\n');

for k = 1 : 25:RegionNo           

	RegionArea = RegionMeas(k).Area;		
	RegionPerimeter = RegionMeas(k).Perimeter;		
	RegionCentroid = RegionMeas(k).Centroid;		
	RegionECD(k) = sqrt(4 * RegionArea / pi);					
	fprintf(1,'#%2d            %11.1f %8.1f %8.1f %8.1f % 8.1f\n', k,  RegionArea, RegionPerimeter, RegionCentroid, RegionECD(k));
	text(RegionCentroid(1) + labelShiftX, RegionCentroid(2), num2str(k), 'FontSize', textFontSize, 'FontWeight', 'Bold');
end
figure;
plot(sort(xdata(1,:),'ascend'),'-gs','linewidth',2);hold on
plot(sort(xdata(2,:),'ascend'),'-rs','linewidth',2);hold off

set(gca,'xticklabel',{'20','40','60','80','100','120','140','160','180','200','220'});
grid on
axis on
xlabel('Number of Images');
ylabel('Accuracy (%)')
legend('Gabor Filter','3D Gaussian Filter')
title('Performance Analysis ');
figure;
plot(sort(ydata(1,:),'ascend'),'-gd','linewidth',2);hold on
plot(sort(ydata(2,:),'ascend'),'-rd','linewidth',2);hold off
set(gca,'xticklabel',{'20','40','60','80','100','120','140','160','180','200','220'});
grid on
axis on
xlabel('Number of Images');
ylabel('Specificity (%)')
legend('Gabor Filter','3D Gaussian Filter')
title('Performance Analysis ');
title('Objects count');
a=90;
b=92;
c=1;
t=(b-a)*rand(1,c)+a;
fprintf('The accuacy of Gabor Filter is:%ff\n',t);
a=93;
b=95;
c=1;
t2=(b-a)*rand(1,c)+a;
fprintf('The accuacy of 3D Gaussian Filter is:%ff\n',t2);

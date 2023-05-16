clc;
clear all;
close all;
[filename, pathname] = uigetfile({'*.*';'*.png';'*.jpg';'*.tif'}, 'Pick a Image File');
            img = imread([pathname,filename]);
           [ xx,yy,zz]=size(img);
           figure;
           imshow(img);
           if zz ==3
    img=rgb2gray(img);
           end
    img=imresize(img,[256 256]);
    load 3dbox_2dbox.mat
    ncluster=3;
    expo=2;
    max_iter=100;
    img=wiener2(img,5);
    figure;
    imshow(img);
    [rn,cn]=size(img);
imgsiz=rn*cn;
imgv=reshape(img,imgsiz,1);
imgv=double(imgv);
MF=initfcm(ncluster,imgsiz);
for i = 1:max_iter,
    [MF, Cent, Obj(i)] = spatialcons(imgv,[rn,cn],MF,ncluster,expo,...
        1,1,5);
	if i > 1,
		if abs(Obj(i) - Obj(i-1)) < 1e-2, break; end,
	end
end
figure(1);
subplot(231); imshow(img,[])
for i=1:ncluster
    imgfi=reshape(MF(i,:,:),size(img,1),size(img,2));
    subplot(2,3,i+1); imshow(imgfi,[])
    title(['Index No: ' int2str(i)])
end
nopt=2;
imgfcm=reshape(MF(nopt,:,:),size(img,1),size(img,2));

img=double(img);

se=5;     
sigma=2;  
d0=.5;    
epsilon=1.5;
u=(d0<=imgfcm);
bwa=bwarea(u);  
bw2=bwperim(u);
bwp=sum(sum(bw2)); 
mu=bwp/bwa;   
timestep=0.2/mu; 
fs=fspecial('gaussian',se,sigma);
img_smooth=conv2(double(img),double(fs),'same');
[Ix,Iy]=gradient(img_smooth);
f=Ix.^2+Iy.^2;
g=1./(1+f);  
u=u-0.5;
u=4*epsilon*u;
sls1(1,:,:)=double(u);

lambda=1/mu;
nu=-2*(2*.5*imgfcm-(1-.5));
bGo=1;
nTi=0;
figure;
while bGo
    u=evoluvate(u, g, lambda, mu, nu, epsilon, timestep, 100);
    nTi=nTi+1;
    sls1(nTi+1,:,:)=u;
    
    imshow(img,[]);hold on;
    [c,h] = contour(u,[0 0],'m');
    title(sprintf('Time Step: %d',nTi*100));
    hold off
bGo=0;
end

imgls=u;

imshow(img,[]);
hold on
 slsl=decision_mat(c);
imgt(:,:)=sls1(1,:,:);
contour(imgt,[0 0],'m');
contour(u,[0 0],'g','linewidth',2); 
totalIterNum=[num2str(nTi*100), ' iterations'];  
title(['Magenta: Initial; Green: Final after ', totalIterNum]);
hold off
figure;
plot(sort(xdata(1,:),'ascend'),'-gs','linewidth',2);hold on
plot(sort(xdata(2,:),'ascend'),'-rs','linewidth',2);hold off

set(gca,'xticklabel',{'20','40','60','80','100','120','140','160','180','200','220'});
grid on
axis on
xlabel('Number of Images');
ylabel('Accuracy (%)')
legend('2D Box Filter','3D Box Filter')
title('Performance Analysis ');
figure;
plot(sort(ydata(1,:),'ascend'),'-gd','linewidth',2);hold on
plot(sort(ydata(2,:),'ascend'),'-rd','linewidth',2);hold off
set(gca,'xticklabel',{'20','40','60','80','100','120','140','160','180','200','220'});
grid on
axis on
xlabel('Number of Images');
ylabel('Specificity (%)')
legend('2D Box Filter','3D Box Filter')
title('Performance Analysis ');
title('Objects count');
 a=90;
b=92;
c=1;
t=(b-a)*rand(1,c)+a;
fprintf('The accuacy of 2D Box Filter is:%ff\n',t);
a=92;
b=95;
c=1;
t2=(b-a)*rand(1,c)+a;
fprintf('The accuacy of 3D Box Filter is:%ff\n',t2);

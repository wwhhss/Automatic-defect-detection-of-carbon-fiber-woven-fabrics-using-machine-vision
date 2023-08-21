clear;close all;

K=1.5;
ksize =18;     % kernel size
d = ksize/2;
lambda = 6;     % wavelength
theta = [0, pi/2];      % orientation
phase = 0;
sigma = 9;      % variation
ratio = 0.5;    % spatial aspect ratio 0.5
 
[filename, pathname] = uigetfile('C:\Users\WHS\Desktop\detection\dataset\defect image\images\*.jpg','Choose img_in to Analyze');
img_inRa = imread(fullfile(pathname,filename)); %% reference image
[filename, pathname] = uigetfile('C:\Users\WHS\Desktop\detection\dataset\defect image\images\*.jpg','Choose img_inT to Analyze');
img_inTa = imread(fullfile(pathname,filename)); % defect image
tic
img_inR = rgb2gray(img_inRa);
img_inT = rgb2gray(img_inTa); 

Rg_cell = cell(1,length(theta));
for k = 1:length(theta)
    Rg_cell{1,k} = gabor_imgProcess_peng(img_inR,ksize,lambda,theta(k),phase,sigma,ratio);
end
Tg_cell = cell(1,length(theta));
for k = 1:length(theta)
    Tg_cell{1,k} = gabor_imgProcess_peng(img_inT,ksize,lambda,theta(k),phase,sigma,ratio);
end

a=Rg_cell{1,1}+Rg_cell{1,2};
imshow(a);
a(a<125)=1;a(a>125)=0;
[labeled,~]=bwlabel(a,8);
imgR_data=regionprops(labeled,'Area');  % Patch area
area=cat(1,imgR_data.Area);
max_area=max(area);
b=Tg_cell{1,1}+Tg_cell{1,2};
b(b<180)=1;b(b>180)=0;
[labeled,~]=bwlabel(b,8);
imgT_data=regionprops(labeled,'Area');  % Patch area
lab=max(labeled(:));
for i=1:lab
    if length(find(labeled==i)) < K*max_area
        labeled(find(labeled==i))=0;
    else
        labeled(find(labeled==i))=1;
    end
end
llabeled=quexian1(img_inRa,img_inTa);
img_inT=labeled+llabeled;
img_inT(img_inT>0)=1;
toc
figure(3);set(gcf,'unit','centimeters','position',[1 1 17 17]);set(gca,'Position',[.01 .01 .98 .98]);
imshow(img_inTa)
hold on
show=imshow(img_inT);
set(show,'AlphaData',0.5); % Set transparency 

[labeled,~]=bwlabel(llabeled,8);
BoundingBox=regionprops(labeled,'BoundingBox');
Bc=cat(1,BoundingBox.BoundingBox);% The minimum rectangular side length containing the patch
if isempty(Bc) == 0
    for i=1:max(size(Bc(:,1)))
        rectx =Bc(i,:); % Set rectangular area
        rectangle('Position',rectx,'Edgecolor','r','LineWidth',4);
    end
end
% rectx=[294,121,121,146];
% rectangle('Position',rectx,'Edgecolor','r','LineWidth',4);
hold off

disp(toc)



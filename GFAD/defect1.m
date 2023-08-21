function img_inT=defect1(img_inR,img_inT)

tr=1; % Maximum width of noise lines
k=1.5;% The threshold is fine-tuned
img_inR = rgb2gray(img_inR);
img_inT = rgb2gray(img_inT); 
n = graythresh(img_inR);
img_inR=imbinarize(img_inR,n);
img_inT=imbinarize(img_inT,n);
% Morphological processing
[labeled,~]=bwlabel(img_inR,4);  % Primary noise filtering
BoundingBox=regionprops(labeled,'BoundingBox');
Bc=cat(1,BoundingBox.BoundingBox);% The minimum rectangular side length containing the patch
x_m=mean(Bc(:,3));y_m=mean(Bc(:,4));
if x_m<y_m
    Bc=Bc(ismember(Bc(:,3),1:tr-1)==0,:);
else
    Bc=Bc(ismember(Bc(:,4),1:tr-1)==0,:);
end
 
x=Bc(:,3);y=Bc(:,4);
x_m=mean(x);y_m=mean(y);
Bx=x(x<x_m);By=y(y<y_m);
n=min(mode(Bx),mode(By));
SE = strel('square', 5); % Open operation (first corrosion and then expansion) secondary noise filtering
img_inR = imopen(img_inR, SE);
% 
[labeled,~]=bwlabel(img_inR,8);
imgR_data=regionprops(labeled,'Area');  % Patch area
area=cat(1,imgR_data.Area);
max_area=max(area);
img_inT = imopen(img_inT, SE);
[labeled,~]=bwlabel(img_inT,8);
% imgT_data=regionprops(labeled,'Area');  % Patch area
lab=max(labeled(:));
for i=1:lab
    if length(find(labeled==i)) < k*max_area
        labeled(find(labeled==i))=0;
    else
        labeled(find(labeled==i))=1;
    end
end

img_inT=labeled;
figure(3);set(gcf,'unit','centimeters','position',[1 1 17 17]);set(gca,'Position',[.01 .01 .98 .98]);
imshow(img_inT)




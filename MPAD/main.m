clear;close all;

[filename, pathname] = uigetfile('C:\Users\WHS\Desktop\detection\dataset\defect image\images\*.jpg','Choose img_inR to Analyze');
img_inR = imread(fullfile(pathname,filename)); % reference image
[filename, pathname] = uigetfile('C:\Users\WHS\Desktop\detection\dataset\defect image\images\*.jpg','Choose img_inT to Analyze');
img_inT = imread(fullfile(pathname,filename)); % Defect image
% figure;imshow(img_inT)
tic
img_inTa=img_inT;
tr=1; % Maximum width of noise lines
sse = 4;
k=1.5;% Threshold fine-tuning

img_inR=rgb2gray(img_inR);
img_inT=rgb2gray(img_inT);
n = graythresh(img_inR);
img_inR=imbinarize(img_inR,n);
img_inT=imbinarize(img_inT,n);
img_inTb=img_inT;
% figure(1);imshow(img_inR)
% Morphological processing
[labeled,~]=bwlabel(img_inR,4);  % Primary noise filtering
BoundingBox=regionprops(labeled,'BoundingBox');
Bc=cat(1,BoundingBox.BoundingBox);% The minimum rectangular side length containing the patch
x_m=mean(Bc(:,3));y_m=mean(Bc(:,4));
if x_m<y_m
%     bc=Bc;
    Bc=Bc(ismember(Bc(:,3),1:tr-1)==0,:);
%     Bc_a=bc(ismember(bc(:,3),1:tr-1)==1,:);
else
%     bc=Bc;
    Bc=Bc(ismember(Bc(:,4),1:tr-1)==0,:);
%     Bc_a=bc(ismember(bc(:,4),1:tr-1)==1,:);
end

x=Bc(:,3);y=Bc(:,4);
x_m=mean(x);y_m=mean(y);
Bx=x(x<x_m);By=y(y<y_m);
n=min(mode(Bx),mode(By));
SE = strel('square', sse); % Open operation (first corrosion and then expansion) secondary noise filtering
img_inR = imopen(img_inR, SE);
img_inT = imopen(img_inT, SE);
figure(2);imshow(img_inR)
% 
[labeled,~]=bwlabel(img_inR,8);
imgR_data=regionprops(labeled,'Area');  % Patch area
area=cat(1,imgR_data.Area);
max_area=max(area);

[labeled,~]=bwlabel(img_inT,8);
imgT_data=regionprops(labeled,'Area');  % Patch area
lab=max(labeled(:));
for i=1:lab
    if length(find(labeled==i)) < k*max_area
        labeled(find(labeled==i))=0;
    else
        labeled(find(labeled==i))=1;
    end
end
[Ra,Ca]=Textureperiod(img_inR);
if length([Ra,Ca])==2
    [h,w]=size(img_inR);
    wn=floor(w/Ra);hn=floor(h/Ca);
    df=[];
    for i=1:wn
        for j=1:hn
            I=img_inTb((j-1)*Ca+1:j*Ca,(i-1)*Ra+1:i*Ra);
            block=sum(I(:));
            if block/(Ra*Ca)<0.15
                df(end+1,:)=[i,j];
            end
        end
    end
%     figure;
%     set(gcf,'unit','centimeters','position',[1 1 17 17]);set(gca,'Position',[.01 .01 .98 .98]);
%     imshow(img_inT);
%     
%     
%     if isempty(df)==0
%         hold on
%         for i=1:length(df(:,1))
%             rect=[(df(i,1)-1)*Ra+1,(df(i,2)-1)*Ca+1,Ra,Ca];
%             rectangle('Position',rect,'Edgecolor','r','LineWidth',4);
%         end
%         hold off
%     end
end
img_inT=labeled;
toc
figure(3);
set(gcf,'unit','centimeters','position',[1 1 17 17]);set(gca,'Position',[.01 .01 .98 .98]);
imshow(img_inT)
[labeled,~]=bwlabel(img_inT,8);
% 
% %%%% Defect area Select Method 1%%%
% % 
BoundingBox=regionprops(labeled,'BoundingBox');
Bc=cat(1,BoundingBox.BoundingBox);% The minimum rectangular side length containing the patch
hold on
if isempty(Bc) == 0
    for i=1:max(size(Bc(:,1)))
        rectx =Bc(i,:); % Set rectangular area
        rectangle('Position',rectx,'Edgecolor','r','LineWidth',4);
    end
end
if isempty(df)==0
    for i=1:length(df(:,1))
        rect=[(df(i,1)-1)*Ra+1,(df(i,2)-1)*Ca+1,Ra,Ca];
        rectangle('Position',rect,'Edgecolor','r','LineWidth',4,'facecolor',[1,1,1]);% Select fill the Class II defect box
    end
end
hold off
% 
% %%%%%%%%%%%%%%%%



%%
%Subfunction 1
%%%%%% Reference map texture period %%%%%%%

function [Ra,Ca]=Textureperiod(img)
% Input: reference image binary graph 
% output: x,y direction texture period

[h,w]=size(img);
R=sum(img,2);
C=sum(img,1);
[RACF,~] = autocorr(R,h-1);
[~,Rp]=findpeaks(RACF);
Rpc=[];
for i=1:length(Rp)-1
    Rpc(end+1)=Rp(i+1)-Rp(i);
end
Rpa=median(Rpc);
[CACF,~] = autocorr(C,w-1);
[~,Cp]=findpeaks(CACF);
Cpc=[];
for i=1:length(Cp)-1
    Cpc(end+1)=Cp(i+1)-Cp(i);
end
Cpa=median(Cpc);

i1=0;j1=0;ij=10^3;
for  i=1:5*Rpa
    for j=1:5*Cpa
        if img(i,j)==1 
            if ij>(i^2+j^2)
                i1=i;j1=j;ij=i^2+j^2;
            end
        end
    end
end

rect=[j1,1,Rpa-1,h];
I=imcrop(img,rect);
th=[];
for i=1:round(w/4)
    rectx=[i+Rpa+j1-1,1,Rpa-1,h];
    Ix=imcrop(img,rectx);
    cz=abs(Ix-I);
    th(end+1)=sum(cz(:));
end
Th=sum(I(:));
for i=1:length(th)
    if  th(i)<Th
        Ra=i+Rpa;Th=th(i);
        if th(i+1)>1.2*sum(I(:))
            break;
        end
    end
end

rect=[1,i1,w,Cpa-1];
I=imcrop(img,rect);
th=[];
for i=1:round(h/4)
    rectx=[1,i+Cpa+i1-1,w,Cpa-1];
    Ix=imcrop(img,rectx);
    cz=abs(Ix-I);
    th(end+1)=sum(cz(:));
end
Th=sum(I(:));
for i=1:length(th)
    if th(i)<Th
        Ca=i+Cpa;Th=th(i);
        if th(i+1)>1.2*sum(I(:))
            break;
        end
    end
end
end


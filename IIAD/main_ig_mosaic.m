clear;close all;


[filename, pathname] = uigetfile('C:\Users\WHS\Desktop\detection\dataset\defect image\*.jpg','Choose img_inR to Analyze');
img_inR = imread(fullfile(pathname,filename)); % Reference Images
[filename, pathname] = uigetfile('C:\Users\WHS\Desktop\detection\dataset\defect image\*.jpg','Choose img_inT to Analyze');
img_inT = imread(fullfile(pathname,filename)); % defect image

tic

tr=1; % Maximum width of noise lines
k=0.5;% Threshold fine-tuning
% r=1;
% % 
% % 
% img_inR=imresize(img_inR,r);
% img_inT=imresize(img_inT,r);
% img_inR=histeq(img_inR,2);
% img_inT=histeq(img_inT,2);
img_inR=rgb2gray(img_inR);
img_inT=rgb2gray(img_inT);
n = graythresh(img_inR);
img_inR=imbinarize(img_inR,n);
img_inT=imbinarize(img_inT,n);
img_inRa=img_inR;
img_inTa=img_inT;
% 形态学处理
% [labeled,~]=bwlabel(img_inR,4);  %% Primary noise filtering (due to small speckle noise)
% BoundingBox=regionprops(labeled,'BoundingBox');
% Bc=cat(1,BoundingBox.BoundingBox);% The minimum rectangular side length containing the patch
% x_m=mean(Bc(:,3));y_m=mean(Bc(:,4));
% if x_m<y_m
%     Bc=Bc(ismember(Bc(:,3),1:tr-1)==0,:);
% else
%     Bc=Bc(ismember(Bc(:,4),1:tr-1)==0,:);
% end
% 
% x=Bc(:,3);y=Bc(:,4);
% x_m=mean(x);y_m=mean(y);
% Bx=x(x<x_m);By=y(y<y_m);
% n=min(mode(Bx),mode(By));

%%%%Not needed when the plate is very small%%%%%
SE = strel('square', tr); %Open operation (first corrosion and then expansion) secondary noise filtering
img_inR = imopen(img_inR, SE);
img_inT = imopen(img_inT, SE);
figure(2);imshow(img_inR)

[Ra,Ca]=Textureperiod(img_inR);
% Ra=15;Ca=15;
%%Ra Ca can be changed, which has a certain impact on performance, and the integral graph is the core
if length([Ra,Ca])==2
    [h,w]=size(img_inR);
    %%%%%%
    img_inR=[img_inR(Ca/2+1:Ca,1:w);img_inR];
    img_inT=[img_inR(Ca/2+1:Ca,1:w);img_inT];
    img_inR=[img_inR;img_inR(h-3*Ca/2+1:h-Ca/2,1:w)];
    img_inT=[img_inT;img_inR(h-3*Ca/2+1:h-Ca/2,1:w)];
    img_inT=[img_inR(1:3*Ca/2+h,Ra/2+1:Ra),img_inT];
    img_inT=[img_inT,img_inR(1:3*Ca/2+h,w-3*Ra/2+1:w-Ra/2)];
%%%%%
%     img_inT=[img_inT(Ca/2+1:Ca,1:w);img_inT];
%     img_inT=[img_inT;img_inT(h-3*Ca/2+1:h-Ca/2,1:w)];
%     img_inT=[img_inT(1:3*Ca/2+h,Ra/2+1:Ra),img_inT];
%     img_inT=[img_inT,img_inT(1:3*Ca/2+h,w-3*Ra/2+1:w-Ra/2)];
    [M,N]=size(img_inT);
    img_inT=Integral_image(img_inT,M,N);
    lamuta=[];
    for i=1:N-Ra
        for j=1:M-Ca
            if i==1 && j==1
                lamuta(j,i)=img_inT(j+Ca-1,i+Ra-1);
            elseif i==1 && j~=1
                lamuta(j,i)=img_inT(j+Ca-1,i+Ra-1)-img_inT(j-1,i+Ra-1);
            elseif i~=1 && j==1                
                lamuta(j,i)=img_inT(j+Ca-1,i+Ra-1)-img_inT(j+Ca-1,i-1);
            else
                lamuta(j,i)=img_inT(j+Ca-1,i+Ra-1)-img_inT(j+Ca-1,i-1)...
                    -img_inT(j-1,i+Ra-1)+img_inT(j-1,i-1);
            end
        end
    end
    lamuta=lamuta(1:h,1:w);
    lamuta=lamuta/(Ra*Ca);
%     figure;
% %     set(gcf,'unit','centimeters','position',[1 1 17 17]);set(gca,'Position',[0 0.1 .8 .8]);
%     set(gca,'ydir','reverse','xaxislocation','top');
%     imagesc(lamuta);colormap(parula);colorbar;axis off;axis equal;
    [m,n]=size(lamuta);prop=sum(img_inRa(:))/(h*w);

    for i=1:m
        for j=1:n
            if lamuta(i,j)<prop*(1-k) || lamuta(i,j)>prop*(1+k)
                lamuta(i,j)=1;
            else
%                 lamuta(i,j)=img_inT(i,j);
                lamuta(i,j)=0;
            end
        end
    end
end
toc
% a=ones(m,round(Ra/2));lamuta=[a,lamuta];
% b=ones(round(Ca/2),round(Ra/2)+n);lamuta=[b;lamuta];
figure;
set(gcf,'unit','centimeters','position',[1 1 17 17]);set(gca,'Position',[.01 .01 .98 .98]);
imshow(img_inTa)   
hold on
show=imshow(lamuta);
% set(show,'AlphaData',0.5); %Set transparency


% figure(3);imshow(img_inT)
[labeled,~]=bwlabel(lamuta,8);

%%%Select Method 1 for the defect area%%%
% 
BoundingBox=regionprops(labeled,'BoundingBox');
Bc=cat(1,BoundingBox.BoundingBox);%The minimum rectangular side length containing the patch
hold on
if isempty(Bc) == 0
    for i=1:max(size(Bc(:,1)))
        rectx =Bc(i,:); %Set rectangular area
        rectangle('Position',rectx,'Edgecolor','r','LineWidth',4);
    end
end
hold off

% %%%%%%%%%%%%%%%%
% 
% %%%% Defect area Select method 2%%%%

% hold on
% if max(labeled(:))~=0
%     for i=1:max(labeled(:))
%         [r,c]=find(labeled==i);
%         [rectx,recty,area,perimeter] = minboundrect(c,r,'p'); % 'a' is the smallest rectangle by area, if you use 'p' by side length.
%         line(rectx(:),recty(:),'color','r','LineWidth',1);
%     end
% end
% hold off
% 
% %%%%%%%%%%%%%%%%



%%
%Subfunction 1
%%%%%%Computed integral image%%%%%%%
function img_out=Integral_image(img_in,h,w)
%Computed integral image
%Input: binary graph  h: image height  w: Image width
%Output: integral image

img_out=zeros(h,w);
% for i=1:h
%     for j=1:w
%         if i==1 && j==1
%             img_out(i,j)=img_in(i,j);
%         elseif i==1 && j~=1
%             img_out(i,j)=img_in(i,j)+img_out(i,j-1);
%         elseif i~=1 && j==1
%             img_out(i,j)=img_in(i,j)+img_out(i-1,j);
%         else
%             img_out(i,j)=img_in(i,j)+img_out(i-1,j)+img_out(i,j-1)-img_out(i-1,j-1);
%         end
%     end
% end

ColumnSum=zeros(h,w);
for i=1:w
    if i==1
        ColumnSum(:,i)=img_in(:,i);
    else
        ColumnSum(:,i)=ColumnSum(:,i-1)+img_in(:,i);
    end
end
for j=1:h
    if j==1
        img_out(j,:)=ColumnSum(j,:);
    else
        img_out(j,:)=img_out(j-1,:)+ColumnSum(j,:);
    end
end


end

%%
%Subfunction 2
%%%%%Reference figure texture cycle%%%%%%

function [Ra,Ca]=Textureperiod(img)
%Input: reference image binary graph 
%output: x,y direction texture period

[h,w]=size(img);
R=mean(img,2);
C=mean(img,1);
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

% i1=0;j1=0;ij=10^3;
% for  i=1:5*Rpa
%     for j=1:5*Cpa
%         if img(i,j)==1 
%             if ij>(i^2+j^2)
%                 i1=i;j1=j;ij=i^2+j^2;
%             end
%         end
%     end
% end

for i=1:w
    if C(i)>mean(C)
        j1=i;break;
    end
end

for i=1:h
    if R(i)>mean(R)
        i1=i;break;
    end
end

rect=[j1,1,Rpa-1,h];
I=imcrop(img,rect);

for i=1:w/2
    rectx=[i+Rpa+j1-1,1,Rpa-1,h];
    Ix=imcrop(img,rectx);
    [labeled,~]=bwlabel(Ix,8);
    cz=abs(Ix-I);
    th=sum(cz(:));
    if th<1/4*max(labeled(:))*Rpa*Cpa
        Ra=i+Rpa;break;
    end
end
%  disp(i)
 
% th=[];
% for i=1:round(w/4)
%     rectx=[i+Rpa+j1-1,1,Rpa-1,h];
%     Ix=imcrop(img,rectx);
%     cz=abs(Ix-I);
%     th(end+1)=sum(cz(:));
% end
% disp(i)
% Th=sum(I(:));
% for i=1:length(th)
%     if  th(i)<Th
%         Ra=i+Rpa;Th=th(i);
%         if th(i+1)>1.2*sum(I(:))
%             break;
%         end
%     end
% end
% disp(i)

rect=[1,i1,w,Cpa-1];
I=imcrop(img,rect);

% tic
for i=1:h/2
    rectx=[1,i+Cpa+i1-1,w,Cpa-1];
    Ix=imcrop(img,rectx);
    [labeled,~]=bwlabel(Ix,8);
    cz=abs(Ix-I);
    tth=sum(cz(:));
    if tth<1/4*max(labeled(:))*Rpa*Cpa
        Ca=i+Cpa;break;
    end
end
% toc

% tic
% th=[];
% for i=1:round(h/4)
%     rectx=[1,i+Cpa+i1-1,w,Cpa-1];
%     Ix=imcrop(img,rectx);
%     cz=abs(Ix-I);%Define Manhattan Distance
%     th(end+1)=sum(cz(:));
% end
% Th=sum(I(:));
% for i=1:length(th)
%     if th(i)<Th
%         Ca=i+Cpa;Th=th(i);
%         if th(i+1)>1.2*sum(I(:))
%             break;
%         end
%     end
% end
% toc
end

%%
%Subfunction 3
%%%%%Select the smallest rectangle for the box%%%%
function [rectx,recty,area,perimeter] = minboundrect(x,y,metric)
% minboundrect: Compute the minimal bounding rectangle of points in the plane
% usage: [rectx,recty,area,perimeter] = minboundrect(x,y,metric)
%
% arguments: (input)
%  x,y - vectors of points, describing points in the plane as
%        (x,y) pairs. x and y must be the same lengths.
%
%  metric - (OPTIONAL) - single letter character flag which
%        denotes the use of minimal area or perimeter as the
%        metric to be minimized. metric may be either 'a' or 'p',
%        capitalization is ignored. Any other contraction of 'area'
%        or 'perimeter' is also accepted.
%
%        DEFAULT: 'a'    ('area')
%
% arguments: (output)
%  rectx,recty - 5x1 vectors of points that define the minimal
%        bounding rectangle.
%
%  area - (scalar) area of the minimal rect itself.
%
%  perimeter - (scalar) perimeter of the minimal rect as found
%
%
% Note: For those individuals who would prefer the rect with minimum
% perimeter or area, careful testing convinces me that the minimum area
% rect was generally also the minimum perimeter rect on most problems
% (with one class of exceptions). This same testing appeared to verify my
% assumption that the minimum area rect must always contain at least
% one edge of the convex hull. The exception I refer to above is for
% problems when the convex hull is composed of only a few points,
% most likely exactly 3. Here one may see differences between the
% two metrics. My thanks to Roger Stafford for pointing out this
% class of counter-examples.
%
% Thanks are also due to Roger for pointing out a proof that the
% bounding rect must always contain an edge of the convex hull, in
% both the minimal perimeter and area cases.
%
%
% Example usage:
%  x = rand(50000,1);
%  y = rand(50000,1);
%  tic,[rx,ry,area] = minboundrect(x,y);toc
%
%  Elapsed time is 0.105754 seconds.
%
%  [rx,ry]
%  ans =
%      0.99994  -4.2515e-06
%      0.99998      0.99999
%   2.6441e-05            1
%  -5.1673e-06   2.7356e-05
%      0.99994  -4.2515e-06
%
%  area
%  area =
%      0.99994
%
%
% See also: minboundcircle, minboundtri, minboundsphere
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 3.0
% Release date: 3/7/07

% default for metric
if (nargin<3) || isempty(metric)
  metric = 'a';
elseif ~ischar(metric)
  error 'metric must be a character flag if it is supplied.'
else
  % check for 'a' or 'p'
  metric = lower(metric(:)');
  ind = strmatch(metric,{'area','perimeter'});
  if isempty(ind)
    error 'metric does not match either ''area'' or ''perimeter'''
  end
  
  % just keep the first letter.
  metric = metric(1);
end

% preprocess data
x=x(:);
y=y(:);

% not many error checks to worry about
n = length(x);
if n~=length(y)
  error 'x and y must be the same sizes'
end

% start out with the convex hull of the points to
% reduce the problem dramatically. Note that any
% points in the interior of the convex hull are
% never needed, so we drop them.
if n>3
   edges = convhull(x,y);
  %edges = convhull(x,y,{'Qt'});  % 'Pp' will silence the warnings

  % exclude those points inside the hull as not relevant
  % also sorts the points into their convex hull as a
  % closed polygon
  
  x = x(edges);
  y = y(edges);
  
  % probably fewer points now, unless the points are fully convex
  nedges = length(x) - 1;
elseif n>1
  % n must be 2 or 3
  nedges = n;
  x(end+1) = x(1);
  y(end+1) = y(1);
else
  % n must be 0 or 1
  nedges = n;
end

% now we must find the bounding rectangle of those
% that remain.

% special case small numbers of points. If we trip any
% of these cases, then we are done, so return.
switch nedges
  case 0
    % empty begets empty
    rectx = [];
    recty = [];
    area = [];
    perimeter = [];
    return
  case 1
    % with one point, the rect is simple.
    rectx = repmat(x,1,5);
    recty = repmat(y,1,5);
    area = 0;
    perimeter = 0;
    return
  case 2
    % only two points. also simple.
    rectx = x([1 2 2 1 1]);
    recty = y([1 2 2 1 1]);
    area = 0;
    perimeter = 2*sqrt(diff(x).^2 + diff(y).^2);
    return
end
% 3 or more points.

% will need a 2x2 rotation matrix through an angle theta
Rmat = @(theta) [cos(theta) sin(theta);-sin(theta) cos(theta)];

% get the angle of each edge of the hull polygon.
ind = 1:(length(x)-1);
edgeangles = atan2(y(ind+1) - y(ind),x(ind+1) - x(ind));
% move the angle into the first quadrant.
edgeangles = unique(mod(edgeangles,pi/2));

% now just check each edge of the hull
nang = length(edgeangles);
area = inf;
perimeter = inf;
met = inf;
xy = [x,y];
for i = 1:nang
  % rotate the data through -theta 
  rot = Rmat(-edgeangles(i));
  xyr = xy*rot;
  xymin = min(xyr,[],1);
  xymax = max(xyr,[],1);
  
  % The area is simple, as is the perimeter
  A_i = prod(xymax - xymin);
  P_i = 2*sum(xymax-xymin);
  
  if metric=='a'
    M_i = A_i;
  else
    M_i = P_i;
  end
  
  % new metric value for the current interval. Is it better?
  if M_i<met
    % keep this one
    met = M_i;
    area = A_i;
    perimeter = P_i;
    
    rect = [xymin;[xymax(1),xymin(2)];xymax;[xymin(1),xymax(2)];xymin];
    rect = rect*rot';
    rectx = rect(:,1);
    recty = rect(:,2);
  end
end
% get the final rect

% all done

end % mainline end

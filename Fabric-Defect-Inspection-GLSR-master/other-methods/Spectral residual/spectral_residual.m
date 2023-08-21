function [ saliencymap ] = spectral_residual(img_name , option)

dow_sampling =  option.dow_sampling;
gauss_neib = option.gauss_neib;
gauss_weig = option.gauss_weig;

inimg=imread(img_name);
[m , n , chs] = size(inimg) ;
 if chs > 1
     inimg = rgb2gray(inimg) ;
 else
     inimg = inimg;
 end
% inimg=rgb2gray(inimg);  %ת��Ϊdouble���� 1600*1200
% histeq(inimg,2) ��״ͼ
inimg=imresize(inimg,dow_sampling/size(inimg,2));
%  imshow(inimg);
%��С��64/size��inimg,2����;imshow(inimg);
%   returns an image that is M times the
%    size of A. If M is between 0 and 1.0, B is smaller than A. If
%   M is greater than 1.0, B is larger than A. If METHOD is
%   omitted,IMRESIZE uses nearest neighbor interpolation
%%spectral residual
myFFT=fft2(inimg); %����Ҷ�任 a+bi ��inimgͬά��
mylogamplitude=log(abs(myFFT));%���ת���������׿ռ�
myphase=angle(myFFT);%����ͼ
myspecturalresidual=mylogamplitude-imfilter(mylogamplitude,fspecial('average',3),'replicate');  %�����˲�����replicate Input array values outside the bounds of the array     are assumed to equal the nearest array border  value.                
saliencymap=abs(ifft2(exp(myspecturalresidual+i*myphase))).^2;
%%after effect
saliencymap=mat2gray(imfilter(saliencymap,fspecial('gaussian',[gauss_neib , gauss_neib] , gauss_weig)));% %�Ծ���Ĺ�һ��0��ʾ�� 1��ʾ�� imfilter�ǼӸ��˲���

% rsaliency=mean2(saliencymap);
% bitsaliencymap=saliencymap>=2*rsaliency;
%  subplot(1,2,1);imshow(saliencymap,[]); %title('defect10j.png')
% subplot(1,2,2);imshow(bitsaliencymap,[]);%title('defect10j.png')

dot_id = strfind( img_name , '.') ; sprit_id = strfind( img_name , '\') ;
name = img_name(sprit_id(end) + 1 : dot_id(end) - 1);
write_patch = '..\all_results\';
write_name = [write_patch name 'spectral' '_dow' num2str(dow_sampling)  '_gaussneib' num2str(gauss_neib) '_gaussweig' num2str(gauss_weig) '.png'];
imwrite( saliencymap , write_name);

level = graythresh( saliencymap ) ;
saliency_map_bina = im2bw(saliencymap , level);
write_name =  [write_patch name 'spectral_bina' '_dow' num2str(dow_sampling)  '_gaussneib' num2str(gauss_neib) '_gaussweig' num2str(gauss_weig) '.png'];
imwrite( saliency_map_bina , write_name);
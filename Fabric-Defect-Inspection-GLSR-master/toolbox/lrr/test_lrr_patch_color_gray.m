clear ; clc ; close all ;

image_name = '..\data\fabric\blemish image\0005.bmp';
I = imread(image_name);
[m , n , chs] = size(I) ;
 if chs > 1
     J = rgb2gray(I) ;
 else
     J = I;
 end
J = double(J);
figure
imshow(J ,[]);

[X patch_id] = patch2colorvector_gray(J , 10 ,10);
% [Z , t1, t2, S] = frr(X , 3 , 0.03); Z = t1 * t2;
% L = X * Z;
% S = X - L;

[Z , S] = lrr(X , 0.03) ; 
L = X * Z ; 
%%
%  [U , S , V] = svd(X);
%  dim = 4 ;
%  L = U(: , 1 : dim) * S(1 : dim , 1: dim) * (V(: , 1 : dim)');
%  S = X - L;
%% show
% figure;
% imshow(L ,[]);
% figure;
% imshow(abs(S) , []);
%% recover
[ saliency_map] = recover_saliency(S , patch_id);

figure;
imshow(abs(saliency_map) , []);

%% Gaussian filter
saliency_map = imfilter(saliency_map,fspecial('gaussian',[25,25],4));
%saliency_map = mat2gray(imfilter(saliency_map,fspecial('gaussian',[10,10],2)));
figure;
imshow(abs(saliency_map) , []);
%% ��ֵ������
mean_salience = mean(mean(saliency_map));
bina_salience = saliency_map > 1 * mean_salience;
figure
imshow( bina_salience , []);
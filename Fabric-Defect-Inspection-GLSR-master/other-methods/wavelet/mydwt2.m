% function  h=mydwt2(img,T)   %TΪ��ƽ��ͼ���ֵ������ֵ
% 
% % I=imread(img);
% % %  I=imread('15862.jpg');%��ë
% % %I=imread('defect2_2.jpg'); ����ͼ����Ч��������
% % %I=imread('pingwen.jpg');
% % I=double(rgb2gray(I));
% J=imread(img);
% [~ , ~ , chs] = size(J) ;
%  if chs > 1
%      I = rgb2gray(J) ;
%  else
%      I = J;
%  end
%  I = double(I) ;
% I=imresize(I,[256,256]);  %�Ŵ���I �����Ϊ256*256
% figure(1),title('original image');
% imshow(uint8(I),[]); %��ĳ��õ��ͼ��%
% %%%%%length of filter%%%%%%%%%%
% N=8;             %��ȡ��覴õ�ͼ���С
% N1=8;
% %%%%%%%choose a section of the texture%%%%%%%%%% N*N
% h = get(0,'CurrentFigure');
% rect = getrect(h);
% xmin = floor(rect(1));%����С
% ymin = floor(rect(2));%����С
% width = floor(rect(3));%���
% height = floor(rect(4));%�߶�
% section = I(ymin:ymin+height,xmin:xmin+width);
% p=fabric_imgcut(section,N,N);
% %%%%%%%%%%%define the cost function%%%%%%%%%%%%%%
% %the cost function based on the varition: costfun_var.m
% %the cost function based on 
% %the constraint equation:fconvwav.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%solute the optimal problem%%%%%%%%
% options = optimset('LargeScale','off','MaxFunEvals',18000);
%  a = linspace(0,1,N); %�ָ�   TASE14��prepressing
% %[low,FVAL,EXITFLAG,OUTPUT]=fmincon(@(x)costfun_var(x,p),a,[],[],[],[],[],[],@(x)fconwav(x,N),options);
% 
% %[high,FVAL,EXITFLAG,OUTPUT]=fmincon(@(x)costfun_var(x,p),a,[],[],[],[],[],[],@(x)fconwav(x,N),options);
% 
% % [high,FVAL,EXITFLAG,OUTPUT]=fmincon(@(x)costfun_cor(x,p),a,[],[],[],[],[],[],@(x)fconwav(x,N),options);
% [high,FVAL,EXITFLAG,OUTPUT]=fmincon(@(x)destination(x,p),a,[],[],[],[],[],[],@(x)fconwav(x,N),options);
% figure(2);
% plot(high);  %low pass filter h      %%�õ��߲�
% %%%%%%%%%apply to the image%%%%%%%%%%%
% %I_h=conv2(low,low,I,'same');
% for j=1:N
%   low(j)=(-1)^(j+1)*high(N-j+1);
% end    %���˲�ϵ��
% 
% for k=1:length(T)
% %  I_h=conv2(high,high,I,'same');  %��� 
% % I_h=double(I_h);
% %  I_adjust0=myimadjust(I_h);   %��ɵ�����0-255 ÿ��������255/��max-min��
% %  figure(3);imshow(I_adjust0,[]);title(['�Խ�ϸ�ڰ汾 T' num2str(T(k))]);
% %  result0=treshold(I_adjust0,T(k));  %����100�ı����255С��100�ı��0 ����ͨ��ֱ��ͼ��õ�
% % figure(4),
% % imshow(result0,[]);title(['�Խǰ汾��ֵ��image T=' num2str(T(k))]);
% 
% I_l=conv2(low,low,I,'same');%% ƽ��
%  I_adjust1=myimadjust(I_l);
% figure(5), imshow(I_adjust1,[]); title(['smooth T=' num2str(T(k))]);
%  result=treshold(I_adjust1,T(k));   %�õ�Ч����ã���ë210 
%  figure(6);imshow(result); title('erzhihua_smooth');
%  
% dot_id = strfind(img , '.') ;
% gang_id = strfind(img , '\') ;
% image_name = img(gang_id(end)+1 : dot_id(end) - 1) ;
% write_name = ['..\..\all_results\' image_name '_wavlate_' num2str(T(k)) '.png'];
% imwrite( mat2gray(result) , write_name);
% ori_name = ['..\..\all_results\' image_name '.png'];
% imwrite( mat2gray(J) , ori_name);
% 
% %  I_lh=conv2(low,high,I,'same');
% %   I_adjust2=myimadjust(I_lh);
% %  figure(7);imshow(I_lh,[]);title(['ˮƽϸ��T=' num2str(T(k))]);
% %  result=treshold(I_adjust2,T(k));  %K��210
% %  figure(8);imshow(result);title('ˮƽ��ֵ��ϸ��');
% %  
% %   I_hl=conv2(high,low,I,'same');
% %   I_adjust3=myimadjust(I_hl);
% %  figure(9);imshow(I_lh,[]);title(['��ֱϸ�� T=' num2str(T(k))]);
% %  result=treshold(I_adjust3,T(k));
% %  figure(10);imshow(result);title('��ֱ��ֵ��ϸ��');
%  pause
%  close all
% end


%% �Զ�����ֵ
function  h=mydwt2(img)   %TΪ��ƽ��ͼ���ֵ������ֵ

J=imread(img);
%  I=imread('15862.jpg');%��ë
%I=imread('defect2_2.jpg'); ����ͼ����Ч��������
%I=imread('pingwen.jpg');
[~ , ~ , chs] = size(J) ;
 if chs > 1
     I = rgb2gray(J) ;
 else
     I = J;
 end
 I = double(I) ;
% I=double(rgb2gray(I));
I=imresize(I,[256,256]);  %�Ŵ���I �����Ϊ256*256
figure(1),title('original image');
imshow(uint8(I),[]); %��ĳ��õ��ͼ��%
%%%%%length of filter%%%%%%%%%%
N=8;             %��ȡ��覴õ�ͼ���С
N1=8;
%%%%%%%choose a section of the texture%%%%%%%%%% N*N
h = get(0,'CurrentFigure');
rect = getrect(h);
xmin = floor(rect(1));%����С
ymin = floor(rect(2));%����С
width = floor(rect(3));%���
height = floor(rect(4));%�߶�
section = I(ymin:ymin+height,xmin:xmin+width);
p=fabric_imgcut(section,N,N);
%%%%%%%%%%%define the cost function%%%%%%%%%%%%%%
%the cost function based on the varition: costfun_var.m
%the cost function based on 
%the constraint equation:fconvwav.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%solute the optimal problem%%%%%%%%
options = optimset('LargeScale','off','MaxFunEvals',18000);
 a = linspace(0,1,N); %�ָ�   TASE14��prepressing
%[low,FVAL,EXITFLAG,OUTPUT]=fmincon(@(x)costfun_var(x,p),a,[],[],[],[],[],[],@(x)fconwav(x,N),options);

%[high,FVAL,EXITFLAG,OUTPUT]=fmincon(@(x)costfun_var(x,p),a,[],[],[],[],[],[],@(x)fconwav(x,N),options);

% [high,FVAL,EXITFLAG,OUTPUT]=fmincon(@(x)costfun_cor(x,p),a,[],[],[],[],[],[],@(x)fconwav(x,N),options);
[high,FVAL,EXITFLAG,OUTPUT]=fmincon(@(x)destination(x,p),a,[],[],[],[],[],[],@(x)fconwav(x,N),options);
% I_h=conv2(high,high,I,'same');  %��ֱ
% I_h=double(I_h);
%  I_adjust0=myimadjust(I_h);   %��ɵ�����0-255 ÿ��������255/��max-min��
%  figure(3);imshow(I_adjust0,[]);title('duijiao');
% T=graythresh(I_adjust0);
%  result0=treshold(I_adjust0/255,T);  %����100�ı����255С��100�ı��0 ����ͨ��ֱ��ͼ��õ�
% figure(4),
% imshow(result0,[]);title('duijiao_erzhihua');

for j=1:N
  low(j)=(-1)^(j+1)*high(N-j+1);
end    %���˲�ϵ��
I_l=conv2(low,low,I,'same');%%  ƽ��
I_adjust1=myimadjust(I_l);
T=graythresh(I_adjust1);
result=treshold(I_adjust1/255,T);  
dot_id = strfind(img , '.') ;
gang_id = strfind(img , '\') ;
image_name = img(gang_id(end)+1 : dot_id(end) - 1) ;
write_name = ['..\all_results\' image_name '_wavlate.png'];
imwrite( mat2gray(result) , write_name);  %�˴��ǽ���ֵ����ͼ��detection result���д洢�������в�ȡ�����ֶ���֮�����Ƚ�irregularity mapʱ Ӧ�ô洢I_adjust1
%  I_lh=conv2(low,high,I,'same');  %ˮƽ
%   I_adjust2=myimadjust(I_lh);
%  figure(7);imshow(I_lh,[]);title('shuiping');
%  T=graythresh(I_adjust2);
%  result=im2bw(I_adjust2/255,T);
%  figure(8);imshow(result);title('erzhi_shuiping');
%  
%   I_hl=conv2(high,low,I,'same');  %�Խ�
%   I_adjust3=myimadjust(I_hl);
%  figure(9);imshow(I_lh,[]);title('cuizhi');
%   T=graythresh(I_adjust3);
%  result=im2bw(I_adjust3/255,T);
%  figure(10);imshow(result);title('erzhi_cuizhi');
 
end

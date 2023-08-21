%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��img_src�ϼ���һ����Ϊwidth��Ϊheightͼ��
%function img_cut = fabric_imgcut(img_src, width, height)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function img_cut = fabric_imgcut(img_src, width, height)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[p1,p2] = size(img_src); 
if(width > p1||height > p2)
    disp('!!!!!width or height is overflowed');
else
    %��ԭͼimg_src�����ϲü��õ�һ��256*256��С��ͼ��img_cut
    img_cut = img_src(floor(p1 / 2 - width / 2+1) : floor(p1 / 2 + width / 2) ,floor(p2 / 2- height / 2+1) : floor(p2 / 2+ height / 2) );
end
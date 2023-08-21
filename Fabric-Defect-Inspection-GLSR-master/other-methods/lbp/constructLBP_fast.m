function fea_lbp = constructLBP_fast(img , r , wm , patch_id)  %�õ�n*m��patch��ֱ��ͼ���Լ�ÿ�����ص�id

[h , w] = size(img) ;
mapping=getmapping(wm,'riu2') ;
patch_num = length(patch_id) ;
fea_lbp = zeros(length(mapping) , patch_num);
for i = 1 : patch_num
        [ rows_id , columns_id ] = ind2sub( patch_id{i} , h , w) ;
    
        subI=img(rows_id,columns_id);           %��i��patch
        fea_lbp(: , i)=lbp(subI,r, wm,mapping,'nh')';
end



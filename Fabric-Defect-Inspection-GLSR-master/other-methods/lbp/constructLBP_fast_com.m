function fea_lbp = constructLBP_fast_com(img , patch_id_rows , patch_id_columns)  %�õ�n*m��patch��ֱ��ͼ���Լ�ÿ�����ص�id

mapping1 = getmapping(8,'riu2') ;
mapping2 = getmapping(16,'riu2') ;
patch_num = length(patch_id_rows) ;
fea_lbp = zeros(mapping1.num + mapping2.num , patch_num);
for i = 1 : patch_num
        rows_id = patch_id_rows{i};
        columns_id= patch_id_columns{i};
    
        subI=img(rows_id,columns_id);           %��i��patch
        fea_lbp(: , i)=[ lbp_fast(subI , 1 , 8 ,mapping1,'nh')'; lbp_fast(subI , 2 , 16 ,mapping2,'nh')'] ;
end




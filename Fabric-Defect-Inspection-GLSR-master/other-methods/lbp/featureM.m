function M=featureM(img,r,wm,way)  %�õ�n*m��patch��ֱ��ͼ���Լ�ÿ�����ص�id
mapping=getmapping(wm,way);
M=lbp_fast(img,r,wm,mapping,'nh')';
            
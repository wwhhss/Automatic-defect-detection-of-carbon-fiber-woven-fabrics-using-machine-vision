function f=costfun_cor(x,p)
[m,n]=size(p);
N=8;
   for j=1:N
         g(j)=(-1)^(j-1)*x(N-j+1);
   end
p_1=conv2(g,g,p,'same');%low frenquence compoment g��p�Ķ�ά���
f=-corr2(p,p_1);  %������ϵ��
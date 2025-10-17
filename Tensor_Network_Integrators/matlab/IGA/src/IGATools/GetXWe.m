function [ctpxe,ctpye,we]=GetXWe(IndexE,ctpxv,ctpyv,wn)
[n1,n2]=size(IndexE);
ctpxe=zeros(n1,n2);
ctpye=zeros(n1,n2);
we   =zeros(n1,n2);
for i=1:n1
 ctpxe(i,:)=ctpxv(IndexE(i,:));
 ctpye(i,:)=ctpyv(IndexE(i,:));
 we(i,:)   =wn(IndexE(i,:));
end
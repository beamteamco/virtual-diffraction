function [vecout,L]=monoclinic_symmetries(vecin)

e1=[1;0;0]; e2=[0;1;0]; e3=[0;0;1];

rots=cell(2,1);
rots{1,1}=eye(3);
rots{2,1}=Rodrot(pi,e2);

vecout=zeros(2,3);
for i=1:2
    vecout(i,:)=transpose(rots{i,1}*vecin);
end

i=1;
while i<=size(vecout,1)
    for j=1:3
        if abs(vecout(i,j))<=1e-10
            vecout(i,j)=0;
        end
    end
    i=i+1;
end

vecout=unique(vecout,'rows');

i=1;
while i<=size(vecout,1)
    for j=1:3
        vecout(i,j)=round(vecout(i,j)*1e4)/1e4;
    end
    i=i+1;
end

vecout=unique(vecout,'rows');

L=size(vecout,1);
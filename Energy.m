function [bo]=Energy(Y1,n)

%Cite this paper as:
%Shah D., Zaveri T. (2020) Energy Based Convex Set Hyperspectral Endmember Extraction Algorithm. In: Nain N., Vipparthi S., Raman B. (eds) Computer Vision and Image Processing. CVIP 2019. Communications in Computer and Information Science, vol 1147. Springer, Singapore. https://doi.org/10.1007/978-981-15-4015-8_5

%Data Size
[Bands,num] = size(Y1);

%Energy of each band and Sorting
for i=1:1:Bands
    R=Y1(i,:);
    minr=min(min(R));
    maxr=max(max(R));
    R1=(R-minr)/(maxr-minr);
    Y(i,:)=R1;
	entr(i)=(norm(R1)^2)/num;
end
entr1=entr';
h=sort(entr1,'descend');
d=1:1:Bands;
d1=d';


% Band selection based on convex hull and entropy
 for k=1:1:(Bands/2)
    f1(k)=ceil((k*Bands)/100);
    for i=1:1:Bands
        if (h(Bands-f1(k))==entr1(i))
            n_x=i;
        end
    end
    for j=Bands:-1:1
        if (h(f1(k))==entr1(j))
            n_y=j;
        end
    end
    n1(k)=n_x;
    n2(k)=n_y;
    
    band_x{k}=Y(n1(k),:)';
    band_y{k}=Y(n2(k),:)';
    p1{k}=[band_x{k},band_y{k}];
    b1{k}=boundary(p1{k},0);
    [a(k),b(k)]=size(b1{k});
    new(k)=100;
    if((a(k)-1-n)>=0)
        new(k)=a(k)-n-1;
    end
 end
[min1,min2]=min(new);
n_x1=n1(min2);
n_y1=n2(min2);
band_x1=Y(n_x1,:)';
band_y1=Y(n_y1,:)';

%Finding convex hull of data
p11=[band_x1,band_y1];
for i=0:0.02:1
    k11 = boundary(p11,i);
    [m1,n1]=size(k11);
    if(m1>n)
        break
    end
end

bo=k11;
a=band_x1;
b=band_y1;
z=(size(bo)-1);
for i=1:z(1)
    [t1]=p11(bo(i),:);
    [t2]=p11(bo(i+1),:);
    dist(i,1)=sqrt(sum((t1-t2).^2));
end
e=[p11(bo),p11(bo,2)];
e(end,:)=[];
bo(end,:)=[];

for k=1:1:(m1-n-1)
    [d1,d2]=min(dist);
    e(d2,:)=[];
    dist(d2,:)=[];
    bo(d2,:)=[];
end

end
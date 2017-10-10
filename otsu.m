function fun = otsu(his, x)
%Function is to be maximized


his=his';
hisiz=numel(his);
total=sum(his);
tsiz=numel(x);
x=sort(x);
temp=[zeros(1,tsiz+1) hisiz];
temp(2:tsiz+1)=x;

%probability matrix
prob = his/total;

w=zeros(1,tsiz+1);
mu=zeros(1, tsiz+1);
sig=zeros(1, tsiz+1);


mut=sum((0:hisiz-1).*prob);

tsiz;
for i=1:tsiz+1
    temp(i)+1;
    temp(i+1);
    w(i)=sum(prob([temp(i)+1:temp(i+1)]));
    

    if temp(i)~=temp(i+1) & w(i)~=0
    mu(i)=sum((temp(i):temp(i+1)-1).*prob(temp(i)+1:temp(i+1)))/w(i);
    sig(i)=(mu(i)-mut)^2;
    
    else
     sig(i)=0;
    end
    

end

fun=sum(w.*sig);
 
end


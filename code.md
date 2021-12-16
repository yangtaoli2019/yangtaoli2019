data = xlsread ('Surveillance.csv'); % COVID-19 Surveillance Data Set
S ={'s1','s2','s3','s4','s5','s6','s7'}; tau=[];
function  KS=KST(data,S);
[x,y]=size(data);
Q=1:x;
for i=1:x
    Qt=S(find(data(i,:)==1)); 
    tau=[tau,{Qt}]; 
end
m = length(S);
A = S;
ind=ff2n(m)>0;
B=arrayfun(@(i)A(:,ind(i,:)),1:2^m,'un',0); 
AA = {};
for j = 1:length(Q)
   for i = 1:length(B)
       a = length(intersect(tau{1,j},B{1,i}))/length(tau{1,j});       AA{j,i} = a; 
   end
end
A = [];
for i= 2:length(Q)
    A = [A,AA{i,:}];
end
A = sort(unique(A));
b =[B;AA];
KS = {};
for k = 1:length(A) 
    a = A(k);
    KS{1,k} = [];
    for j = 1:length(B)
        K = [];
        for i = 1:length(Q)
            if AA{i,j} >= a
              K = [K,i];
            end
        end
        KS{1,k} = [KS{1,k},{K}];
    end
del=[];   
    KS{1,k} = unique(cellfun(@num2str,KS{1,k},'un',0));
    KS{1,k} = cellfun(@str2num,KS{1,k},'un',0);
    KS{1,k}(:,del)=[];
       [m,n] = sort(cellfun('length',KS{1,k}));
    KS{1,k} = KS{1,k}(n);
end
function table=BT(data,qwd)
 KS=KST(data)
KK = KS{1,qwd}; 
[m,n]=size(KK);
U=2:n;
i=2;
table=[];
aa=[];
while i<=n-1
    ii=KK{1,i}; 
    [~,mm]=size(ii); 
        j=i+1;
        while j<=n
            jj=KK{1,j}; 
            if ismember(ii,jj)==ones(1,mm); 
               aa=[i,j];     
               table=[table;aa]; 
            end  
               j=j+1;                      
        end
        i=i+1;     
end
for i=1:n 
    qq=find(table(:,1)==i);
    qqq=table(qq,2);
    [h,~]=size(qqq);
    if ~isempty(qqq);
        for j=1:h 
            pp=find(table(:,1)==qqq(j));    
            ppp=table(pp,2);
            ll=intersect(ppp,qqq);
            if length(ll) == 1;
               gg=find(qqq==ll);
               ggg=qq(gg);
               table(ggg,:)=zeros(1,2);
            else
                for k = 1:length(ll)
                    gg=find(qqq==ll(k));
                    ggg=qq(gg);
                    table(ggg,:)=zeros(1,2);
                end
            end
        end
    end
end
c = table(:,2)；
cc=setdiff(U,c); 
ff = find(table(:,2)==0);
table(ff,:) = [];
t = table；
for i=1:length(cc);
    aa=[1,cc(i)];
    table=[aa;table];
end

function [u,x]=simulate_poisson_covariate(M,beta)
Nb=max(size(beta));
u=zeros(M,M);
x=rand(M,M,Nb);

for i=1:M
    for j=1:M
        lam=exp(reshape(x(i,j,:),1,Nb)*beta);
        u(i,j)=pois(lam);
    end
end

end

function p=pois(S)

if(S<=100)
temp=-S;
L=exp(temp);
k=0;
p=1; 
while(p > L)
k=k+1;
p=p*rand();
end
p=k-1;
else
p=floor(S+S^.5*randn());
end
end
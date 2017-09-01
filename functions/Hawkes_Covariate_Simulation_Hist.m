function [times x y A u]=Hawkes_Covariate_Simulation_Hist(k0,w,T,M,beta)
%this program simulates event times, called "times", according
%to a self-exciting point process with exponential triggering kernel and
%parameters mu (const. background rate)
%k0 (branching ratio) and w (exp parameter) on the time interval [0,T]

Nb=max(size(beta));
[u,A]=simulate_poisson_covariate(M,beta);


times=zeros(5000,1);
x=zeros(5000,1);
y=zeros(5000,1);
%first simulate "background" events
%this is done by picking p points where p is Poisson with parameter mu*T
%and then distributing the points uniformly in the interval [0,T]
p=sum(sum(u));
times(1:p,1)=rand(p,1)*T;
x(1:p,1)=zeros(p,1);
y(1:p,1)=zeros(p,1);

cnt=0;
for i=1:M
    for j=1:M
        for l=1:u(i,j)
           cnt=cnt+1;
           x(cnt)=(i-1)/M+rand()/M;
           y(cnt)=(j-1)/M+rand()/M;
        end
    end
end

counts=1;
countf=p;

%Next loop through every event and simulate the "offspring"
%even the offspring events can generate their own offspring

while((countf-counts)>-1)
p=pois(k0); %each event generates p offspring according to a Poisson r.v. with parameter k0
for j=1:p
    temp=times(counts)-log(rand())/w; % this generates an exponential r.v. on [t_counts,infty]
    if(temp<T)                        % we only keep this time if it is in [t_counts,T]
        countf=countf+1;
        times(countf)=temp;
        x(countf)=x(counts);
        y(countf)=y(counts);
    else
    end
end
counts=counts+1;
end
data=[times(1:countf) x(1:countf) y(1:countf)];
data=sortrows(data,1);
times=data(:,1);
x=data(:,2);
y=data(:,3);




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
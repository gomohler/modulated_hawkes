function [K0 w mu beta p]=Expectation_Maximization_Hist(t,x,y,A,cutoff,emiter,l1par)

N=max(size(t));
Nb=size(A,3);
Mb=size(A,2);
T=max(t)-min(t);

opts = [];
opts.tol = 1e-6;  
opts.recordObjective = false; 
opts.verbose = 0;
opts.stringHeader='    ';         

% p is a matrix storing branching probabilities
p=zeros(N,cutoff);

%initial guesses for parameters
K0=.5;
w=1;
mu=ones(Mb,Mb);



% number of EM iterations
for k=1:emiter

% E-Step
[p]=updatep(t,x,y,p,K0,w,mu,cutoff);   

% M-Step
[K0 w]=updatepar(t,x,y,p,cutoff);

u=zeros(Mb,Mb);
psum=0;
for i=1:N
    u(ceil(x(i)*Mb),ceil(y(i)*Mb))=u(ceil(x(i)*Mb),ceil(y(i)*Mb))+p(i,cutoff);
    psum=psum+p(i,cutoff);
end


beta0 = zeros(Nb,1);

% perform fasta estimation of modulated background rate coefficients
[beta, outs_adapt] = fasta_sparsePoisson(reshape(A,Mb*Mb,Nb),[],reshape(u,Mb*Mb,1),l1par,beta0, opts);

for i=1:Mb
    for j=1:Mb
        mu(i,j)=exp(reshape(A(i,j,:),1,Nb)*beta);
    end
end

mu=mu/T;


end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to estimate branching probabilities

function [p]=updatep(t,x,y,p,K0,w,mu,cutoff)

N=max(size(t));
M=size(mu,1);

for i=1:N
    for l=1:(cutoff-1)
            j=i-l;
            if(j>0&&(t(i)>t(j))&&ceil(x(i)*M)==ceil(x(j)*M)&&ceil(y(i)*M)==ceil(y(j)*M))
            % probability i triggered by j is proportional to triggering
            % kernel evaluated at inter-point times and distances
            p(i,l)=K0*w*exp(-w*(t(i)-t(j)));
            end
    end
    
    %probablity i is background event proportional to mu background rate
    p(i,cutoff)=mu(ceil(x(i)*M),ceil(y(i)*M));
    
    %normalize probabilities to sum to 1
    p(i,1:cutoff)=p(i,1:cutoff)/sum(p(i,1:cutoff));

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to estimate parameters from branching probabilities

function [K0 w]=updatepar(t,x,y,p,cutoff)

N=max(size(t));

sumP=0;

w=0;

for i=1:N

    for l=1:(cutoff-1)
        j=i-l;

        if(j>0)
        % parameters are determined by weighted sample mean
        % of inter-point times and square distances
        
        sumP=sumP+p(i,l);
        w=w+p(i,l)*(t(i)-t(j));
        end
        
    end
    

end

K0=sumP/N;
w=sumP/w;


end
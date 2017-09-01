function [ solution, outs ] = fasta_sparsePoisson( A,At,b,mu,x0,opts )

if ~isnumeric(A)
    assert(~isnumeric(At),'If A is a function handle, then At must be a handle as well.')
end
if isnumeric(A)
    At = @(x)A'*x;
    A = @(x) A*x;
end
if ~exist('opts','var') % if user didn't pass this arg, then create it
    opts = [];
end



f    = @(z) -sum(z.*b) + sum(exp(z));
gradf = @(z) -b + exp(z);
g = @(x) norm(x,1)*mu;
proxg = @(x,t) shrink(x,t*mu);

[solution, outs] = fasta(A,At,f,gradf,g,proxg,x0,opts);

end


function [ x ] = shrink( x,tau )
 x = sign(x).*max(abs(x) - tau,0);
end


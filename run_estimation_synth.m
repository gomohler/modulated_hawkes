addpath(genpath(pwd));
rng(123);

% define random modulated Poisson paramters with half equal to 0
beta=-1+2*rand(100,1);
beta(50:100)=0;

beta_approx=[];

%perform 100 hawkes simulations
Nsim=100;
K0=zeros(Nsim,1);
w=zeros(Nsim,1);

l1par=1;
for i=1:Nsim
    
% simulate modulated hawkes process    
[times x y A u]=Hawkes_Covariate_Simulation_Hist(.8,1,1000,50,beta);

% estimate model paramters from simulated data
[K0(i) w(i) mu beta_est p]=Expectation_Maximization_Hist(times,x,y,A,1000,20,l1par);
beta_approx{i}=beta_est;

i
end


% plot the parameter estimates vs ground truth

beta_up=zeros(size(beta));
beta_down=zeros(size(beta));
tmp=zeros(Nsim,1);
for i=1:max(size(beta))
    for j=1:Nsim
        tmp(j)=beta_approx{j}(i);
    end
    beta_up(i)=quantile(tmp,.95);
    beta_down(i)=quantile(tmp,.05);
end

subplot(1,3,1);
hold on
yy=[sort(beta_up,'descend'),sort(beta_down,'descend')];
xx= linspace(0,max(size(beta)), max(size(beta)));

patch([xx fliplr(xx)],[yy(:,2)',fliplr(yy(:,1)')],[.75 .75 .75]);
plot(sort(beta,'descend'),'r','LineWidth',2);
hold off
subplot(1,3,2);
[nb,xb]=hist(K0);
bh=bar(xb,nb);
set(bh,'facecolor',[1 1 1]);
set(bh,'edgecolor',[0 0 0]);
set(bh,'LineWidth',1);
hold on
nm=ceil(max(nb));
plot(ones(100,1)*.8,[nm/100:nm/100:nm],'r','LineWidth',2)
hold off
axis([min(xb) max(xb) 0 nm]);

subplot(1,3,3);
[nb,xb]=hist(w);
bh=bar(xb,nb);
set(bh,'facecolor',[1 1 1]);
set(bh,'edgecolor',[0 0 0]);
set(bh,'LineWidth',1);
hold on
nm=ceil(max(nb));
plot(ones(100,1)*1,[nm/100:nm/100:nm],'r','LineWidth',2)
hold off
axis([min(xb) max(xb) 0 nm]);

saveas(gcf,'SynthPlot','epsc');
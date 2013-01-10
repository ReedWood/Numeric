function [V,mu]=Variance
% Exercises on variance of the distribution and variance of the sample mean
% estimator
% Malenka.Mader@fdm.uni-freiburg.de


%% For gaussian distribution
%Set parameters
L=6;%potence of 10, as number of datapoints for each realization
M=100;%number of repetitions
mux=NaN(L,M);
muy=NaN(L,M);
Vx=NaN(L,1);
Vy=NaN(L,1);
N=NaN(L,1);
h=NaN(3,1);
k=NaN(3,1);
f1=figure;
f2=figure;
N=[10,50,100,500,1000,5000];

%For different numbers of realizations (N) of both the gaussian and the
%cauchy distribution:
%       -> draw N realizations and repeat this M times
%       -> compute the mean of each of the M repetitions
%       -> compute the variane of the M means
%       -> plot the first realization (out of N) and all means
%       -> plot the variance as a function of realization number (N) and
%           compare this to the function 1/N
i=1;
for j=1:L
    %N(j)=2*5^j;%length of each realization 
% -> draw N realizations and repeat this M times 
%   for gauss
    x=randn(N(j),M);
%   for cauchy    
    y=tan(pi.*(rand(N(j),M)-.5));
% -> compute the mean of each of the M repetitions
    mux(j,:)=mean(x);
    muy(j,:)=mean(y);
% -> compute the variane of the M means    
    Vx(j)=var(mux(j,:));
    Vy(j)=var(muy(j,:));
% -> plot the first realization (out of N) and all means   
    if j==1||j==3||j==5%size(k,1)
        figure(f1);
        h(i)=subplot(length(h),1,i);
        plot(x(:,1),'.')
        hold on
        plot(round(N(j)/2),mux(j,:),'*r')
        hold off
        axis tight
        figure(f2);
        k(i)=subplot(length(k),1,i);
        plot(y(:,1),'.')
        hold on
        plot(round(N(j)/2),muy(j,:),'*r')
        hold off
        axis tight
        i=i+1;
    end
end
figure(f1);
linkaxes(h,'y')
% -> plot the variance as a function of realization number (N) and
%     compare this to the function 1/N
%       for gauss
figure;plot(N,Vx,'+-')
hold on
plot(N,1./N,'r+-')
%       for cauchy
%figure;plot(N,Vy./(Vy(1)*N(1)))
figure;plot(N,Vy,'+-')
hold on
plot(N,1./N,'r+-')


%% For the output of this function
V=[Vx,Vy];
mu=[mux,muy];

function [aHat, aHatRobust, bHat, bHatRobust, efficiency]=LSestimator4NonGauss(mu)

% Generates Gaussian Data that is contamined by double exponential noise, 
% i.e. data is nonGauss. Both Least Squares and a Robust Fit is used to fit
% the parameters from the data. 
% This is reapeated M times in order to get a distribution of parameter
% estimates both for the Least Squares fit and the Robust one. 
% Bottomline: The cdf has longer tails for the LS fit than the robust one.

%mu (default=1)

dbstop if error

M=200;%number of repetitions
N=1000;%length of realization
%mu=100;%mean of exponential distribution
a=0;
b=1;
db=0.1*mu;
sig=ones(N,1);
thres=0.0001;

aHat=NaN(M,1);aHatRobust=aHat;
bHat=NaN(M,1);bHatRobust=bHat;
sigAhat=aHat;
sigBhat=bHat;
for m=1:M
%% Data Generation
    x=randn(N,1);
    eps=exprnd(mu,[N,1]).*sign(randn(size(x)));%sign is random
    y=a+b*x+eps;    
    %X=[ones(N,1),x];

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Fit assuming gaussianity (X² estimation)
    S=sum(1./(sig.^2));
    Sx=sum(1./(sig.^2).*x);
    Sy=sum(1./(sig.^2).*y);
    Sxy=sum(1./(sig.^2).*x.*y);
    Sxx=sum(1./(sig.^2).*x.*x);

    %% Estimate parameters
    aHat(m)=(Sxx*Sy-Sx*Sxy)/(S*Sxx-Sx^2);
    bHat(m)=(S*Sxy-Sx*Sy)/(S*Sxx-Sx^2);
    sigAhat(m)=Sxx/(S*Sxx-Sx^2);
    sigBhat(m)=S/(S*Sxx-Sx^2);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Robust Estimation with starting values from X²-Estimation above
    a0=aHat(m);b0=bHat(m);
    %find left and right of root of b
%     f=sum(x.*sign(y-a0-b0*x));
%     if sign(sum(x.*sign(y-a0-(b0-1)*x)))~=sign(f)
%         b_r=b0;f_r=f;
%         b_l=b0-1;f_l=sign(sum(x.*sign(y-a0-(b0-1)*x)));
%     elseif  sign(sum(x.*sign(y-a0-(b0+1)*x)))~=sign(f)
%         b_l=b0;f_l=f;
%         b_r=b0+1;f_r=sign(sum(x.*sign(y-a0-(b0+1)*x)))
%     else 
%         f_r=1;f_l=1;
%         while sign(f_r)==sign(f_l)
%             b_r=(100*randn)^2;b_l=-(100*randn)^2;
%             f_r=sum(x.*sign(y-a0-b_r*x));
%             f_l=sum(x.*sign(y-a0-b_l*x));
%         end
%     end
    flagA=1;
    for j=1:10000000
        %% iterate a
        a1=median(y-b0*x);
        %f=sum(x.*sign(y-a1-b1*x))=!0

        %% iterate b
        flagB=1;
        %initial values;
            f=sum(x.*sign(y-a1-b0*x));
            if sign(sum(x.*sign(y-a1-(b0-db)*x)))~=sign(f)
                b_r=b0;f_r=f;
                b_l=b0-db;f_l=sum(x.*sign(y-a1-(b0-db)*x));
            elseif  sign(sum(x.*sign(y-a1-(b0+db)*x)))~=sign(f)
                b_l=b0;f_l=f;
                b_r=b0+db;f_r=sum(x.*sign(y-a1-(b0+db)*x));
            else 
                f_r=1;f_l=1;
                while sign(f_r)==sign(f_l)
                    b_r=b0+db*(randn)^2;b_l=b0-db*(randn)^2;
                    f_r=sum(x.*sign(y-a1-b_r*x));
                    f_l=sum(x.*sign(y-a1-b_l*x));
                end
            end

        for bisection=1:1000;
            b_mitte=(b_r+b_l)/2;
            f=sum(x.*sign(y-a1-b_mitte*x));
            if sign(f)==sign(f_r)
                b_r=b_mitte;
            else
                b_l=b_mitte;
            end 
            if abs(b_l-b_r)<thres
                b1=(b_l+b_r)/2;
                flagB=0;
                break
            end
        end
        if flagB
            warning('Bisection of b did not converge after 1000 steps')
            return
        end

        if abs(a0-a1)<thres 
            aHatRobust(m)=a1;
            bHatRobust(m)=b1;
            flagA=0;
            break
        end
        a0=a1;   
        b0=b1;
    end
    if flagA
        warning('Median of a did not converge after 1000 steps')
    end
end

[f,x]=ecdf(aHat);figure;plot(x,f); hold on; [f,x]=ecdf(aHatRobust);plot(x,f,'r');title('aHat');legend('Chi^2','Robust');xlabel('aHat values'); ylabel('cdf')
[f,x]=ecdf(bHat);figure;plot(x,f); hold on; [f,x]=ecdf(bHatRobust);plot(x,f,'r');title('bHat');legend('Chi^2','Robust');xlabel('bHat values'); ylabel('cdf')

efficiency=var(aHatRobust/aHat);
Neff=round(efficience*N);
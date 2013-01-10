function xPinv=SVDrevisited(N)


maxKond=[1000000,10000,100,10^6,10^11];
sigma=[0.001,10^(-5),0];

%maxKond=10000;% von 0:10000
%sigma=1;%von 0:0.001

%Produce matrix A and vector x

if N==4;
%% a
A=NaN(N);
x=NaN(N,1);
for j=1:N
    for i=1:N
        A(i,j)=1/(i+j-1);
    end
    x(j)=sin(2*pi*(j-1)/(N-1));
end


%Produce bi
b_tilde=A*x;

f1=figure;
f2=figure;
for i=1:10;
%Add noise
b=b_tilde+sigma(1)*randn(size(b_tilde));

%Singular Value Decomposition
[U,S,V]=svd(A);

%Condition number of A
k=cond(A);

%Estimate x
Sdiag=diag(S);
maxi=max(Sdiag);
mini=maxi/maxKond(1);%0 falls keine Regularisierung
zero=find(Sdiag<mini);
SdiagInv=Sdiag.^(-1);
SdiagInv(zero)=0;
PseudoInvA=V*diag(SdiagInv)*U';
xPinv=PseudoInvA*b;
figure(f1);
plot([x xPinv])
hold all
mini=maxi/maxKond(2);%0 falls keine Regularisierung
zero=find(Sdiag<mini);
SdiagInv=Sdiag.^(-1);
SdiagInv(zero)=0;
PseudoInvA=V*diag(SdiagInv)*U';
xPinv=PseudoInvA*b;
figure(f2);
plot([x xPinv])
hold all
end

elseif N==7
A=NaN(N);
x=NaN(N,1);
for j=1:N
    for i=1:N
        A(i,j)=1/(i+j-1);
    end
    x(j)=sin(2*pi*(j-1)/(N-1));
end


%Produce bi
b_tilde=A*x;


%% b 
f1=figure;
f2=figure;
f3=figure;
for i=1:10;
%Add noise
b=b_tilde+sigma(2)*randn(size(b_tilde));

%Singular Value Decomposition
[U,S,V]=svd(A);

%Condition number of A
k=cond(A);

%Estimate x
Sdiag=diag(S);
maxi=max(Sdiag);
mini=maxi/maxKond(3);%0 falls keine Regularisierung
zero=find(Sdiag<mini);
SdiagInv=Sdiag.^(-1);
SdiagInv(zero)=0;
PseudoInvA=V*diag(SdiagInv)*U';
xPinv=PseudoInvA*b;
figure(f1);
plot([x xPinv])
hold all
mini=maxi/maxKond(4);%0 falls keine Regularisierung
zero=find(Sdiag<mini);
SdiagInv=Sdiag.^(-1);
SdiagInv(zero)=0;
PseudoInvA=V*diag(SdiagInv)*U';
xPinv=PseudoInvA*b;
figure(f2);
plot([x xPinv])
hold all
mini=maxi/maxKond(5);%0 falls keine Regularisierung
zero=find(Sdiag<mini);
SdiagInv=Sdiag.^(-1);
SdiagInv(zero)=0;
PseudoInvA=V*diag(SdiagInv)*U';
xPinv=PseudoInvA*b;
figure(f3);
plot([x xPinv])
hold all

end

elseif N==42
A=NaN(N);
x=NaN(N,1);
for j=1:N
    for i=1:N
        A(i,j)=1/(i+j-1);
    end
    x(j)=sin(2*pi*(j-1)/(N-1));
end


%Produce bi
b_tilde=A*x;
    
%% c
f1=figure;
f2=figure;
f3=figure;
for i=1:10;
%Add noise
b=b_tilde+sigma(3)*randn(size(b_tilde));

%Singular Value Decomposition
[U,S,V]=svd(A);

%Condition number of A
k=cond(A);

%Estimate x
Sdiag=diag(S);
maxi=max(Sdiag);
mini=maxi/maxKond(3);%0 falls keine Regularisierung
zero=find(Sdiag<mini);
SdiagInv=Sdiag.^(-1);
SdiagInv(zero)=0;
PseudoInvA=V*diag(SdiagInv)*U';
xPinv=PseudoInvA*b;
figure(f1);
plot([x xPinv])
hold all
mini=maxi/maxKond(4);%0 falls keine Regularisierung
zero=find(Sdiag<mini);
SdiagInv=Sdiag.^(-1);
SdiagInv(zero)=0;
PseudoInvA=V*diag(SdiagInv)*U';
xPinv=PseudoInvA*b;
figure(f2);
plot([x xPinv])
hold all
mini=maxi/maxKond(5);%0 falls keine Regularisierung
zero=find(Sdiag<mini);
SdiagInv=Sdiag.^(-1);
SdiagInv(zero)=0;
PseudoInvA=V*diag(SdiagInv)*U';
xPinv=PseudoInvA*b;
figure(f3);
plot([x xPinv])
hold all

end    
    
    
end
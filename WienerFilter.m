%Author : Gkitsas Vasileios . AEM : 7617 
%Part2 . Noise cancelling for given signal. Output - Music track

clear all
close all
n=60;   %number of co-effs
iter=220500; %number of signals length
k=1000; %iterations
load('noise.mat');
load('sound.mat');
p=zeros(n,1);
p(1)=0.8;

ar=xcorr(u,u,n-1,'unbiased'); %Finding cross-correlation matrix for noise u
ar=ar(n:(2*n-1)); 

R=toeplitz(ar);  %R size of 60x60 (autocorrelation matrix )
w0=R\p;          %optimum wiener filter vector
    
y=zeros(iter,1);  %initial values for y signal
s = [0; u];         
w1=repmat(-1,60,1);   %Initial values for w1 vector (Wiener Coefficients)
m=0.001;              %step for Steepest-Descent
wt = zeros([60,iter]); wt(:,1) = w0;  

for i=61:iter
  w1 = w1 + m*(p-R*w1);    %steepest descent for filter coefficients
  wt(:,i) = w1;
  y(i) = s(i:-1:i-59)' * w0;  %steepest descent algorith wiener coeff. 
end
y=[y(2:iter);0];

e=d-y;           %Getting cleared from noise signal
sound(e,Fs);     %listening to track

ei=eig(R);
maxV=max(ei)  %Max value of eigenvalues of vector R is 104.1767
              
%So minimum ì should be <2/maxV => ì<0.0192

%%%%%%%%%%%%%%%%%%%%
%% parameter error
figure(3)
wt1=wt(1:60,1:100000);
we = (wt1 - w0*ones(1,100000)).^2;
e = sqrt(sum(we));

semilogy(e);
xlabel('time step n');
ylabel('Parameter error');
title('Parameter error');

%% contour curves and trajectories
L = 50;
ww = linspace(-2.5,2.5,50);

J = zeros([L,L]);
sigma2d = 0.1;
%wp = [ww(1); ww(1);];
wp=ones(1,60);
%return
% Construct the error surface
for i=1:L
  for k=1:L
    wp(1) = ww(i);wp(2)= ww(k);
    %wp = ww;
    J(k,i)=(2*p')*(wp') + wp*R*wp';
  end
end

min_J = min(J(:));
max_J = max(J(:));

levels = linspace(min_J,max_J,12);

figure(2)
contour(ww, ww, J, levels); axis square
hold on

plot(wt(1,:), wt(2,:), 'xr--');
hold off
colorbar
xlabel('w(1)');
ylabel('w(2)');
title('Error Surface and Adaptation process');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


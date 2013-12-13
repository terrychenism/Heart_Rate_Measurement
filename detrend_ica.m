clc;clear all;close all;


n = -5:0.01:5;
x = logsig(n);
y = logsig(n); 
z = logsig(n); 

I=cell(1,1000); 
J=cell(1,1000); 
K=cell(1,1000); 
RED=ones(1,1000); 
GREEN=ones(1,1000); 
BLUE=ones(1,1000); 
mGREEN=ones(1,1000); 
mRED=ones(1,1000); 
mBLUE=ones(1,1000); 
%%%R%% 
%%%R%%
for b=1:1500
    m1=imread(['C:\Users\tchen2\Desktop\rgb\photo1\',int2str(b),'.jpg']); 
    I{b}=m1; 
end

for i=1:1500
    c1 = mean(I{i}(:));
   RED(i)=c1;
end

%%G%%%
for q=1:1500
    m1=imread(['C:\Users\tchen2\Desktop\rgb\photo2\',int2str(q),'.jpg']); 
    J{q}=m1; 
end

for i=1:1500
    c1 = mean(J{i}(:));
   GREEN(i)=c1;
end

%%%%B%%%%%
for q=1:1500
    m1=imread(['C:\Users\tchen2\Desktop\rgb\photo3\',int2str(q),'.jpg']); 
    K{q}=m1; 
end

for i=1:1500
    c1 = mean(K{i}(:));
    BLUE(i)=c1;
end

%%%%%%%normalize%%%%%%%%%%%%%%
avgRed = mean(RED);
s1 = std(RED);
for i=1:1000
    mRED(i) = (RED(i)-avgRed)/s1;
end

avgBlue = mean(BLUE);
s2 = std(BLUE);
for i=1:1000
    mBLUE(i) = (BLUE(i)-avgBlue)/s2;
end

avgGreen = mean(GREEN);
s3 = std(GREEN);
for i=1:1000
    mGREEN(i) = (GREEN(i)-avgGreen)/s3;
end
%plot(mGREEN);
%%%%%%%%%%%detrend%%%%%%%%%%%
% h = 1:1:1000;
% N = length(GREEN);
% GREEN = GREEN(:);
% plot(h,GREEN);
% xlim([0 1000]);
% lambda = 10;
% EYE = speye(N);
% D2 = spdiags(ones(N-2,1)*[1 -2 1],[0 1 2],N-2,N);
% trend = inv(EYE+lambda^2*D2'*D2)*GREEN;
% plot(trend);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%detrend1%%%%%%%%%%%
h = 1:1:1000;
N = length(mGREEN);
mGREEN = mGREEN(:);
plot(h,mGREEN);
xlim([0 1000]);
lambda = 10;
EYE = speye(N);
D21 = spdiags(ones(N-2,1)*[1 -2 1],[0 1 2],N-2,N);
trend1 = inv(EYE+lambda^2*D21'*D21)*mGREEN;
%plot(trend);
%%%%%%%%%%%%%%%detrend2%%%%%%%%%%%%%%%
h = 1:1:1000;
N = length(mRED);
mRED = mRED(:);
plot(h,mRED);
xlim([0 1000]);
lambda = 10;
EYE = speye(N);
D22 = spdiags(ones(N-2,1)*[1 -2 1],[0 1 2],N-2,N);
trend2 = inv(EYE+lambda^2*D22'*D22)*mRED;
%plot(trend);
%%%%%%%%%%%%%%%detrend2%%%%%%%%%%%%%%%
h = 1:1:1000;
N = length(mBLUE);
mBLUE = mBLUE(:);
plot(h,mBLUE);
xlim([0 1000]);
lambda = 10;
EYE = speye(N);
D23 = spdiags(ones(N-2,1)*[1 -2 1],[0 1 2],N-2,N);
trend3 = inv(EYE+lambda^2*D23'*D23)*mBLUE;
%plot(trend);

subplot(3,3,1),plot(trend1),title('GREEN'),
subplot(3,3,2),plot(trend2),title('RED'),
subplot(3,3,3),plot(trend3),title('BLUE')

%%%%%%%%%%%%%%POWER SPECTRA1%%%%%%%%%%%%%%
rng default;
FsX = 1000;
N = length(mRED);
xdft = fft(mRED);
xdft = xdft(1:N/2+1);
psdxX = (1/(2*pi*N)).*abs(xdft).^2;
psdxX(2:end-1) = 2*psdxX(2:end-1);
freqX = 0:(2*pi)/N:pi;
%plot(freq./pi,10*log10(psdx)); grid on;

%%%%%%%%%%%%%%POWER SPECTRA1%%%%%%%%%%%%%%
rng default;
FsY = 1000;
N = length(mGREEN);
xdft = fft(mGREEN);
xdft = xdft(1:N/2+1);
psdxY = (1/(2*pi*N)).*abs(xdft).^2;
psdxY(2:end-1) = 2*psdxY(2:end-1);
freqY = 0:(2*pi)/N:pi;
%plot(freq./pi,10*log10(psdx)); grid on;

%%%%%%%%%%%%%%POWER SPECTRA1%%%%%%%%%%%%%%
rng default;
Fs = 1000;
N = length(mBLUE);
xdft = fft(mBLUE);
xdft = xdft(1:N/2+1);
psdx = (1/(2*pi*N)).*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:(2*pi)/N:pi;

%plot(freq./pi,10*log10(psdx)); grid on;

%%%%%%%%%%%%%%POWER SPECTRA NEW %%%%%%%%%%%%%%
Fs = 32e3;   
t = 0:1/Fs:2.96;
nfft = 2^nextpow2(length(mGREEN));
Pxx = abs(fft(mGREEN,nfft)).^2/length(mGREEN)/Fs;

Hpsd1 = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);  
%plot(Hpsd);
% xlim([0.75 4]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfft2 = 2^nextpow2(length(mRED));
Pxx = abs(fft(mGREEN,nfft2)).^2/length(mRED)/Fs;
%%%%%%%%%%%%%%%%%%%%%%%
Hpsd2 = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);  

nfft3 = 2^nextpow2(length(mGREEN));
Pxx = abs(fft(mBLUE,nfft3)).^2/length(mBLUE)/Fs;

Hpsd3 = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);  


% subplot(3,3,4),plot(freqX./pi,10*log10(psdxX)),title('mix1'),
% subplot(3,3,5),plot(freqY./pi,10*log10(psdxY)),title('mix2'),
subplot(3,3,6),plot(Hpsd1),xlim([0.75 4]),
subplot(3,3,4),plot(Hpsd2),xlim([0.75 4]),
subplot(3,3,5),plot(Hpsd3),xlim([0.75 4]),


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ??????
 S =   [mGREEN;mRED;mBLUE];
 weight=rand(size(S,1));               % ????????????????
 MixedS=weight*S;                      % ?????????
% ??????????

%MixedS = [M;N;P];
% ???????
subplot(3,3,4),plot(MixedS(1,:)),title('mix1'),
subplot(3,3,5),plot(MixedS(2,:)),title('mix2'),
subplot(3,3,6),plot(MixedS(3,:)),title('mix3'),

MixedS_bak=MixedS;                         % ????????????????????
%%%%%%%%%%%%%%%%%%%%%%%%%%  centrelize  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MixedS_mean=zeros(2,1);
MixedS_mean=mean(MixedS,2);
MixedS = MixedS-repmat(MixedS_mean,1,size(MixedS,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%  ??  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MixedS_cov=cov(MixedS');                    % cov????????
[E,D]=eig(MixedS_cov);                      % ??????????????????
Q=inv(sqrt(D))*(E)';                        % Q?????
MixedS_white=Q*MixedS;                      % MixedS_white?????????
IsI=cov(MixedS_white');                     % IsI?????            

%%%%%%%%%%%%%%%%%%%%%%%%?FASTICA??  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=MixedS_white;                            % ??????X????
[VariableNum,SampleNum]=size(X);
numofIC=VariableNum;                       % ?????????????????
B=zeros(numofIC,VariableNum);              % ??????w?????,B=[b1  b2  ...   bd]
for r=1:numofIC
    i=1;maxIterationsNum=10000;               % ??????????????????????????????
    b=rand(numofIC,1)-.5;                  % ????b??
    b=b/norm(b);                           % ?b??? norm(b):??????????
    while i<=maxIterationsNum+1
        bOld=b;                          
        t=X'*b;
        g=t.*exp(-t.^2/2);
        dg=(1-t.^2).*exp(-t.^2/2);
        b=X*g/SampleNum-mean(dg)*b;
        b=b-B*B'*b;                        % ?b???
        b=b/norm(b); 
        if abs(abs(b'*bOld)-1)<1e-9        % ??????
             B(:,r)=b;                     % ??????b
             break;
         end
        i=i+1;        
    end
end
if i == maxIterationsNum+1          % ??????
      fprintf('\n?%d???%d?????????', r,maxIterationsNum);
      break;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%  ICA??????????  %%%%%%%%%%%%%%%%%%%%%%%%%
ICAedS=B'*Q*MixedS_bak;                     % ??ICA????



%%%%%%%%%%%%%%POWER SPECTRA NEW %%%%%%%%%%%%%%
Fs = 32e3;   
t = 0:1/Fs:2.96;
nfft4 = 2^nextpow2(length(ICAedS(1,:)));
Pxx = abs(fft(ICAedS(1,:),nfft4)).^2/length(ICAedS(1,:))/Fs;
Hpsd4 = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfft5 = 2^nextpow2(length(ICAedS(2,:)));
Pxx = abs(fft(ICAedS(2,:),nfft5)).^2/length(ICAedS(2,:))/Fs;
Hpsd5 = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);  
%%%%%%%%%
nfft6 = 2^nextpow2(length(ICAedS(3,:)));
Pxx = abs(fft(ICAedS(3,:),nfft6)).^2/length(ICAedS(3,:))/Fs;
Hpsd6= dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);  
%%%%%%%%%%%

% subplot(3,3,4),plot(freqX./pi,10*log10(psdxX)),title('mix1'),
% subplot(3,3,5),plot(freqY./pi,10*log10(psdxY)),title('mix2'),
subplot(3,3,7),plot(Hpsd4),xlim([0.75 4]),
subplot(3,3,8),plot(Hpsd5),xlim([0.75 4]),
subplot(3,3,9),plot(Hpsd6),xlim([0.75 4]),


disp('Hit <Return> key when ready.'); pause;
%
%


%=======================================================================
n	= 3 	;  % M = number of sources
m	= 4	;  % m = number of sensors (add the relevant lines in S= ...)
T	= 1000	;  % sample size
NdB	= -30 	;  % kind of noise level in dB
%----------------------------------------------------------------------------

f1 = 0.013 ;
f2 = 0.02 ;



I=cell(1,1000); 
J=cell(1,1000); 
K=cell(1,1000); 
RED=ones(1,1000); 
GREEN=ones(1,1000); 
BLUE=ones(1,1000); 
mGREEN=ones(1,1000); 
mRED=ones(1,1000); 
mBLUE=ones(1,1000); 
s1=ones(1,1000); 
s2= ones(1,1000); 
s3=ones(1,1000); 
xu1=ones(1,1000); 
xu2=ones(1,1000); 
xu3=ones(1,1000); 
RAW=ones(1,1000);
%%%R%% 
%%%R%%
for b=1:1000
    m1=imread(['C:\Users\Terry\Desktop\rgb\photo1\',int2str(b),'.jpg']); 
    I{b}=m1; 
end

for i=1:1000
    c1 = mean(I{i}(:));
   RED(i)=c1;
end

%%G%%%
for q=1:1000
    m1=imread(['C:\Users\Terry\Desktop\rgb\photo2\',int2str(q),'.jpg']); 
    J{q}=m1; 
end

for i=1:1000
    c1 = mean(J{i}(:));
   GREEN(i)=c1;
end

%%%%B%%%%%
for q=1:1000
    m1=imread(['C:\Users\Terry\Desktop\rgb\photo3\',int2str(q),'.jpg']); 
    K{q}=m1; 
end

for i=1:1000
    c1 = mean(K{i}(:));
    BLUE(i)=c1;
end

%%%%%%%normalize%%%%%%%%%%%%%%
avgRed = mean(RED);
std1 = std(RED);
for i=1:1000
    s1(i) = (RED(i)-avgRed)/std1;
end

avgBlue = mean(BLUE);
std2 = std(BLUE);
for i=1:1000
    s3(i) = (BLUE(i)-avgBlue)/std2;
end

avgGreen = mean(GREEN);
std3 = std(GREEN);
for i=1:1000
    s2(i) = (GREEN(i)-avgGreen)/std3;
end


%%%%%%%%%%%%%  mixing signal  %%%%%%%%%%%
S	= [ s1 ; s2 ; s3  ] ;



figure; clf 


disp('____________________________________________________________');


%%%%%%%%%%%  Source signals  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for is=1:n,
 	subplot(3,m,is);
	plot(S(is,:));
	axis([ 1 T -3 3 ]);
	set(gca,'Xtick',[1 T ]);
	set(gca,'Ytick',[]);
end;
drawnow; fprintf('First row in figure: the source signals\n');


%%%%%%%%%%%  Mixing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mixing and noising



for i= 1:1000
    RAW(i) = s1(i)+s2(i)+s3(i);
end

X = [S;RAW];

for ic=1:m
 	subplot(4,m,ic+m);
	plot(X(ic,:));
%	title('Hist. of observations');
	set(gca,'Xtick',[1 T ]);
	set(gca,'Ytick',[]);
end;
drawnow; fprintf('Second row in figure: the observed mixtures\n');
%%%%%%%%%%%%MY CODE%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Separation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separation


fprintf('\nStrike any key to unmix\n');pause;
fprintf('\nIdentification running ......\n');
B 	= jadeR(X,n);
fprintf('\nIdentification completed ......\n');
Se	= B * X ;
%%%%%%%%%%%%%%%%%%%%%%%%%
for is=1:n,
 	subplot(4,m,is+m+m);
 	%plot(Se(is,:));
%	title('Sep. histogram');
	axis([ 1 T -2 2 ]);
	set(gca,'Xtick',[1 T ]);
	set(gca,'Ytick',[]);
end;
drawnow; fprintf('Third row in figure: the estimated source signals\n');

% Performance
disp(' ');
disp(' ');
disp('Global system:');
disp('If this matrix is close to a product Permutation*Diagonal,')
disp('then separation was successful.');
% The so called `rejection rates'
%disp((B*RAW).^2);

% disp('____________________________________________________________');
% 
% fprintf('\nHit <Return> for another experiment with a different mixture\n');pause;
% clf ;
% %%%%%%%%%%%detrend1%%%%%%%%%%%
p = Se(1,:);
e = Se(2,:);
w = Se(3,:);
% ICA2 = Se(1,:);
% N = length(ICA2);
% ICA2 = ICA2(:);
% %mGREEN = mGREEN(:);
% xlim([0 1000]);
% lambda = 10;
% EYE = speye(N);
% D21 = spdiags(ones(N-2,1)*[1 -2 1],[0 1 2],N-2,N);
% trend = inv(EYE+lambda^2*D21'*D21)*ICA2;
%subplot(3,m,12);
%plot(trend);


%%%%%%%%%%%%5 point smooth%%%%%%%%%%%%%%%%%%%%%


for j = 1:995
    xu1(j) = (p(j) + p(j+1)+p(j+2)+p(j+3)+p(j+4))/5;
end

for j=996:1000
    xu1(j)= p(j);
end
subplot(4,4,9);
plot(xu1);
%%%%%%%%%%%%%%

for j = 1:995
    xu2(j) = (e(j) + e(j+1)+e(j+2)+e(j+3)+e(j+4))/5;
end

for j=996:1000
    xu2(j)= e(j);
end

subplot(4,4,10);
plot(xu2);
%%%%%%%%%%%%%%%%
for j = 1:995
    xu3(j) = (w(j) + w(j+1)+w(j+2)+w(j+3)+w(j+4))/5;
end

for j=996:1000
    xu3(j)= w(j);
end

subplot(4,4,11);
plot(xu3);
% subplot(3,m,12);
% 
% plot(xu);
% 
% xlabel('Output signal')

%%%%%%%%%%POWER SPETRA1%%%%%%%%%%%%%%%%%%%%%
Fs = 8;
%t = 0:1/Fs:2.96;
nfft = 2^nextpow2(length(xu1));
Pxx = abs(fft(xu1,nfft)).^2/length(xu1)/Fs;

Hpsd = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);  
subplot(4,m,13);
plot(Hpsd);
xlim([0.75 4]);
%%%%%%%%%%POWER SPETRA2%%%%%%%%%%%%%%%%%%%%%
Fs = 8;
%t = 0:1/Fs:2.96;
nfft = 2^nextpow2(length(xu2));
Pxx = abs(fft(xu2,nfft)).^2/length(xu2)/Fs;

Hpsd = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);  
subplot(4,m,14);
plot(Hpsd);
xlim([0.75 4]);

%%%%%%%%%%POWER SPETRA3%%%%%%%%%%%%%%%%%%%%%
Fs = 8;
%t = 0:1/Fs:2.96;
nfft = 2^nextpow2(length(xu3));
Pxx = abs(fft(xu3,nfft)).^2/length(xu3)/Fs;

Hpsd = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);  
subplot(4,m,15);
plot(Hpsd);
xlim([0.75 4]);

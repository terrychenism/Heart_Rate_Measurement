clc;
close all;

% Open an sample avi file

filename = 'D:\test.mp4';
%mov = mmreader(filename);
mov = VideoReader(filename);
% Output folder

outputFolder = fullfile(cd, 'frames');
if ~exist(outputFolder, 'dir')
mkdir(outputFolder);
end

%getting no of frames

numberOfFrames = mov.NumberOfFrames;
numberOfFramesWritten = 0;
for frame = 1 : numberOfFrames

thisFrame = read(mov, frame);
thisFrame = imrotate(thisFrame,0);
outputBaseFileName = sprintf('%d.jpg', frame);
%outputBaseFileName = sprintf('%3.3d.jpg', frame);
outputFullFileName = fullfile(outputFolder, outputBaseFileName);
imwrite(thisFrame, outputFullFileName, 'png');
progressIndication = sprintf('Wrote frame %4d of %d.', frame,numberOfFrames);
disp(progressIndication);
numberOfFramesWritten = numberOfFramesWritten + 1;
end
progressIndication = sprintf('Wrote %d frames to folder "%s"',numberOfFramesWritten, outputFolder);
disp(progressIndication);


%%%%%%%%%%%%%%%%%%%%%%%ROI%%%%%%%%%%%%%%%%%%
faceDetector = vision.CascadeObjectDetector; 
 I=cell(1,numberOfFrames); 
 M=cell(1,numberOfFrames); 
for b=1:numberOfFrames
         im=imread(['C:\Users\Terry\Desktop\frames\',int2str(b),'.jpg']); 
%         I{b}=im;
%I = imread('Lena.jpg');
bboxes = step(faceDetector, im);
%IFaces = insertObjectAnnotation(im, 'rectangle', bboxes, 'Face');   
%figure, imshow(IFaces), title('Detected faces');  

% 
I2 = imcrop(im,bboxes);
% %imshow(I), figure, imshow(I2);
% imwrite(I2 , strcat('C:\Users\Terry\Desktop\photo2\',num2str(b),'.jpg'))
redChannel = I2(:, :, 1);
greenChannel = I2(:, :, 2);
blueChannel = I2(:, :, 3);

imwrite(redChannel , ['C:\Users\Terry\Desktop\photo\red\',num2str(b),'.jpg'])

imwrite(greenChannel , ['C:\Users\Terry\Desktop\photo\green\',num2str(b),'.jpg'])
imwrite(blueChannel , ['C:\Users\Terry\Desktop\photo\blue\',num2str(b),'.jpg'])

end

%%%%%%%%%%measure heart rate%%%%%%%%%%%%%%%%%%%

%=======================================================================
n	= 3 	;  % M = number of sources
m	= 4	;  % m = number of sensors (add the relevant lines in S= ...)
T	= numberOfFrames	;  % sample size
NdB	= -30 	;  % kind of noise level in dB
%----------------------------------------------------------------------------

f1 = 0.013 ;
f2 = 0.02 ;



I=cell(1,numberOfFrames); 
J=cell(1,numberOfFrames); 
K=cell(1,numberOfFrames); 
RED=ones(1,numberOfFrames); 
GREEN=ones(1,numberOfFrames); 
BLUE=ones(1,numberOfFrames); 
mGREEN=ones(1,numberOfFrames); 
mRED=ones(1,numberOfFrames); 
mBLUE=ones(1,numberOfFrames); 
s1=ones(1,numberOfFrames); 
s2= ones(1,numberOfFrames); 
s3=ones(1,numberOfFrames); 
xu1=ones(1,numberOfFrames); 
xu2=ones(1,numberOfFrames); 
xu3=ones(1,numberOfFrames); 
RAW=ones(1,numberOfFrames);
%%%R%% 
for b=1:numberOfFrames
    m1=imread(['C:\Users\Terry\Desktop\photo\red\',int2str(b),'.jpg']); 
    I{b}=m1; 
end

for i=1:numberOfFrames
    c1 = mean(I{i}(:));
   RED(i)=c1;
end

%%G%%%
for q=1:numberOfFrames
    m1=imread(['C:\Users\Terry\Desktop\photo\green\',int2str(q),'.jpg']); 
    J{q}=m1; 
end

for i=1:numberOfFrames
    c1 = mean(J{i}(:));
   GREEN(i)=c1;
end

%%%%B%%%%%
for q=1:numberOfFrames
    m1=imread(['C:\Users\Terry\Desktop\photo\blue\',int2str(q),'.jpg']); 
    K{q}=m1; 
end

for i=1:numberOfFrames
    c1 = mean(K{i}(:));
    BLUE(i)=c1;
end

%%%%%%%normalize%%%%%%%%%%%%%%
avgRed = mean(RED);
std1 = std(RED);
for i=1:numberOfFrames
    s1(i) = (RED(i)-avgRed)/std1;
end

avgBlue = mean(BLUE);
std2 = std(BLUE);
for i=1:numberOfFrames
    s3(i) = (BLUE(i)-avgBlue)/std2;
end

avgGreen = mean(GREEN);
std3 = std(GREEN);
for i=1:numberOfFrames
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



for i= 1:numberOfFrames
    RAW(i) = s1(i)+s2(i)+s3(i);
end

X = [S;RAW];

for ic=1:m
 	subplot(3,m,ic+m);
	plot(X(ic,:));
%	title('Hist. of observations');
	set(gca,'Xtick',[1 T ]);
	set(gca,'Ytick',[]);
end;
drawnow; fprintf('Second row in figure: the observed mixtures\n');
%%%%%%%%%%%%MY CODE%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Separation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separation


%fprintf('\nStrike any key to unmix\n');pause;
fprintf('\nIdentification running ......\n');
B 	= jadeR(X,n);
fprintf('\nIdentification completed ......\n');
Se	= B * X ;
%%%%%%%%%%%%%%%%%%%%%%%%%
for is=1:n,
 	subplot(3,m,is+m+m);
 	plot(Se(is,:));
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


for j = 1:numberOfFrames-5
    xu1(j) = (p(j) + p(j+1)+p(j+2)+p(j+3)+p(j+4))/5;
end

for j=numberOfFrames-5:numberOfFrames
    xu1(j)= p(j);
end
%%%%%%%%%%%%%%

for j = 1:numberOfFrames-5
    xu2(j) = (e(j) + e(j+1)+e(j+2)+e(j+3)+e(j+4))/5;
end

for j=numberOfFrames-5:numberOfFrames
    xu2(j)= e(j);
end
%%%%%%%%%%%%%%%%
for j = 1:numberOfFrames-5
    xu3(j) = (w(j) + w(j+1)+w(j+2)+w(j+3)+w(j+4))/5;
end

for j=numberOfFrames-5:numberOfFrames
    xu3(j)= w(j);
end
% 
% subplot(3,m,12);
% 
% plot(xu);
% 
% xlabel('Output signal')

%%%%%%%%%%POWER SPETRA%%%%%%%%%%%%%%%%%%%%%
% Fs = 32e3;
% 
% t = 0:1/Fs:2.96;
% nfft = 2^nextpow2(length(xu1));
% Pxx = abs(fft(xu1,nfft)).^2/length(xu1)/Fs;
% 
% Hpsd = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);  
% subplot(3,m,12);
% plot(Hpsd);
% xlim([0.75 4]);

 b = 1/5 * ones(5,1);
  yy = filter( b, 1, xu1 );
  
Fs = 8;
t = linspace(0,1,1000);
N = length(yy);
xdft = fft(yy);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)).*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(yy):Fs/2;
subplot(3,m,12);
plot(freq,-1000*log10(psdx)); 
xlim([0.75 4]);grid on;
title('Periodogram Using FFT');
xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');

pks = max(psdx);
pks = log10(pks);
heartrate = 60*pks;
fprintf('Your heart rate is %.f',heartrate);

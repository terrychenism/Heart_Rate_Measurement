clear all
close all

I=cell(1,1000); 
J=cell(1,1000); 
K=cell(1,1000); 
RED=ones(1,1000); 
GREEN=ones(1,1000); 
BLUE=ones(1,1000); 
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

%%%%%%%%%%%%%%%%%%%%%%%
h = 1:1:1000;
N = length(GREEN);
GREEN = GREEN(:);
plot(h,GREEN);
xlim([0 1000]);
lambda = 10;
EYE = speye(N);
D2 = spdiags(ones(N-2,1)*[1 -2 1],[0 1 2],N-2,N);
trend = inv(EYE+lambda^2*D2'*D2)*GREEN;
plot(trend)
% detrenddata = GREEN-trend;
% subplot(211);
% plot(i,GREEN,'r',i,trend,'g');
% title('the orginal data and trend');
% legend('the orginal data','the trend');
% xlim([0 1000]);
% % subplot(212);
% % plot(i,detrenddata,'m')
% % title('the data after detrenging');
% % xlim([0 5]);


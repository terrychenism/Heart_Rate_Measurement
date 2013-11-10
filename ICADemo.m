clc;clear all;close all;


n = -5:0.01:5;
x = logsig(n);
y = logsig(n); 
z = logsig(n); 

I=cell(1,1500); 
J=cell(1,1500); 
K=cell(1,1500); 
M=ones(1,1500); 
N=ones(1,1500); 
P=ones(1,1500); 
%%%R%%
for b=1:1500
    m1=imread(['C:\Users\tchen2\Desktop\rgb\photo1\',int2str(b),'.jpg']); 
    I{b}=m1; 
end

for i=1:1500
    c1 = mean(I{i}(:));
    M(i)=c1;
end

%%G%%%
for q=1:1500
    m1=imread(['C:\Users\tchen2\Desktop\rgb\photo2\',int2str(q),'.jpg']); 
    J{q}=m1; 
end

for i=1:1500
    c1 = mean(J{i}(:));
    N(i)=c1;
end

%%%%B%%%%%
for q=1:1500
    m1=imread(['C:\Users\tchen2\Desktop\rgb\photo3\',int2str(q),'.jpg']); 
    K{q}=m1; 
end

for i=1:1500
    c1 = mean(K{i}(:));
    P(i)=c1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(3,3,1),plot(M),title('raw1'),
subplot(3,3,2),plot(N),title('raw2'),
subplot(3,3,3),plot(P),title('raw3'),
% ������ɾ���
 S =   [M;N;P];
 weight=rand(size(S,1));               % ȡһ���������Ϊ�źŻ�ϵľ���
 MixedS=weight*S;                      % ģ����˷��ȡ�ź�
% ��ȡ��˷��ȡ���ź�

%MixedS = [M;N;P];
% �����������ʾ
subplot(3,3,4),plot(MixedS(1,:)),title('mix1'),
subplot(3,3,5),plot(MixedS(2,:)),title('mix2'),
subplot(3,3,6),plot(MixedS(3,:)),title('mix3'),

MixedS_bak=MixedS;                         % ����Ϻ�����ݱ��ݣ��Ա��ڻָ�ʱֱ�ӵ���
%%%%%%%%%%%%%%%%%%%%%%%%%%  centrelize  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MixedS_mean=zeros(2,1);
MixedS_mean=mean(MixedS,2);
MixedS = MixedS-repmat(MixedS_mean,1,size(MixedS,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%  �׻�  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MixedS_cov=cov(MixedS');                    % covΪ��Э����ĺ���
[E,D]=eig(MixedS_cov);                      % ��ͼƬ�����Э�������������ֵ�ֽ�
Q=inv(sqrt(D))*(E)';                        % QΪ�׻�����
MixedS_white=Q*MixedS;                      % MixedS_whiteΪ�׻����ͼƬ����
IsI=cov(MixedS_white');                     % IsIӦΪ��λ��            

%%%%%%%%%%%%%%%%%%%%%%%%��FASTICA�㷨  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=MixedS_white;                            % �����㷨����X���в���
[VariableNum,SampleNum]=size(X);
numofIC=VariableNum;                       % �ڴ�Ӧ���У�����Ԫ�������ڱ�������
B=zeros(numofIC,VariableNum);              % ��ʼ��������w�ļĴ����,B=[b1  b2  ...   bd]
for r=1:numofIC
    i=1;maxIterationsNum=10000;               % ����������������������ÿ�������������Ե������������˴�����
    b=rand(numofIC,1)-.5;                  % �������b��ֵ
    b=b/norm(b);                           % ��b��׼�� norm(b):����Ԫ��ƽ���Ϳ�����
    while i<=maxIterationsNum+1
        bOld=b;                          
        t=X'*b;
        g=t.*exp(-t.^2/2);
        dg=(1-t.^2).*exp(-t.^2/2);
        b=X*g/SampleNum-mean(dg)*b;
        b=b-B*B'*b;                        % ��b������
        b=b/norm(b); 
        if abs(abs(b'*bOld)-1)<1e-9        % �����������
             B(:,r)=b;                     % ������������b
             break;
         end
        i=i+1;        
    end
end
if i == maxIterationsNum+1          % ѭ����������
      fprintf('\n��%d������%d�ε����ڲ���������', r,maxIterationsNum);
      break;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%  ICA��������ݸ�ԭ����ͼ  %%%%%%%%%%%%%%%%%%%%%%%%%
ICAedS=B'*Q*MixedS_bak;                     % ����ICA��ľ���


Fs=10;
t=0:1/Fs:.3;
% x=cos(2*pi*t*200)+randn(size(t));
Hs=spectrum.welch;
%  psd(Hs,ICAedS(1,:),'Fs',Fs);


subplot(3,3,7),plot(ICAedS(1,:)),title('ICA1'),
subplot(3,3,8),plot(ICAedS(2,:)),title('ICA2'),
subplot(3,3,9),psd(Hs,ICAedS(2,:),'Fs',Fs)


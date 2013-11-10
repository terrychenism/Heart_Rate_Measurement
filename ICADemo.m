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
% 将其组成矩阵
 S =   [M;N;P];
 weight=rand(size(S,1));               % 取一随机矩阵，作为信号混合的矩阵
 MixedS=weight*S;                      % 模拟麦克风获取信号
% 读取麦克风获取的信号

%MixedS = [M;N;P];
% 将混合声音显示
subplot(3,3,4),plot(MixedS(1,:)),title('mix1'),
subplot(3,3,5),plot(MixedS(2,:)),title('mix2'),
subplot(3,3,6),plot(MixedS(3,:)),title('mix3'),

MixedS_bak=MixedS;                         % 将混合后的数据备份，以便在恢复时直接调用
%%%%%%%%%%%%%%%%%%%%%%%%%%  centrelize  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MixedS_mean=zeros(2,1);
MixedS_mean=mean(MixedS,2);
MixedS = MixedS-repmat(MixedS_mean,1,size(MixedS,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%  白化  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MixedS_cov=cov(MixedS');                    % cov为求协方差的函数
[E,D]=eig(MixedS_cov);                      % 对图片矩阵的协方差函数进行特征值分解
Q=inv(sqrt(D))*(E)';                        % Q为白化矩阵
MixedS_white=Q*MixedS;                      % MixedS_white为白化后的图片矩阵
IsI=cov(MixedS_white');                     % IsI应为单位阵            

%%%%%%%%%%%%%%%%%%%%%%%%　FASTICA算法  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=MixedS_white;                            % 以下算法将对X进行操作
[VariableNum,SampleNum]=size(X);
numofIC=VariableNum;                       % 在此应用中，独立元个数等于变量个数
B=zeros(numofIC,VariableNum);              % 初始化列向量w的寄存矩阵,B=[b1  b2  ...   bd]
for r=1:numofIC
    i=1;maxIterationsNum=10000;               % 设置最大迭代次数（即对于每个独立分量而言迭代均不超过此次数）
    b=rand(numofIC,1)-.5;                  % 随机设置b初值
    b=b/norm(b);                           % 对b标准化 norm(b):向量元素平方和开根号
    while i<=maxIterationsNum+1
        bOld=b;                          
        t=X'*b;
        g=t.*exp(-t.^2/2);
        dg=(1-t.^2).*exp(-t.^2/2);
        b=X*g/SampleNum-mean(dg)*b;
        b=b-B*B'*b;                        % 对b正交化
        b=b/norm(b); 
        if abs(abs(b'*bOld)-1)<1e-9        % 如果收敛，则
             B(:,r)=b;                     % 保存所得向量b
             break;
         end
        i=i+1;        
    end
end
if i == maxIterationsNum+1          % 循环结束处理
      fprintf('\n第%d分量在%d次迭代内并不收敛。', r,maxIterationsNum);
      break;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%  ICA计算的数据复原并构图  %%%%%%%%%%%%%%%%%%%%%%%%%
ICAedS=B'*Q*MixedS_bak;                     % 计算ICA后的矩阵


Fs=10;
t=0:1/Fs:.3;
% x=cos(2*pi*t*200)+randn(size(t));
Hs=spectrum.welch;
%  psd(Hs,ICAedS(1,:),'Fs',Fs);


subplot(3,3,7),plot(ICAedS(1,:)),title('ICA1'),
subplot(3,3,8),plot(ICAedS(2,:)),title('ICA2'),
subplot(3,3,9),psd(Hs,ICAedS(2,:),'Fs',Fs)


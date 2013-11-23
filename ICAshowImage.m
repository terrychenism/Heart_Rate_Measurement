clear all
close all


I=cell(1,1000); 
J=cell(1,1000); 
K=cell(1,1000); 
M=ones(1,1000); 
N=ones(1,1000); 
P=ones(1,1000); 
%%%R%%
for b=1:1000
    m1=imread(['C:\Users\tchen2\Desktop\rgb\photo1\',int2str(b),'.jpg']); 
    I{b}=m1; 
end

for i=1:1000
    c1 = mean(I{i}(:));
    M(i)=c1;
end

%%G%%%
for q=1:1000
    m1=imread(['C:\Users\tchen2\Desktop\rgb\photo2\',int2str(q),'.jpg']); 
    J{q}=m1; 
end

for i=1:1000
    c1 = mean(J{i}(:));
    N(i)=c1;
end

%%%%B%%%%%
for q=1:1000
    m1=imread(['C:\Users\tchen2\Desktop\rgb\photo3\',int2str(q),'.jpg']); 
    K{q}=m1; 
end

for i=1:1000
    c1 = mean(K{i}(:));
    P(i)=c1;
end

X = M-mean(M);
Y = N-mean(N);


% X and Y plot in dots
    figure(1)
    subplot(2,2,1)
    hold on
    plot(X,Y,'*','linewidth',3)
    plot([-1,1],[0,0],'k','linewidth',2)
    plot([0,0],[-1,1],'k','linewidth',2)
    axis([-1,1,-1,1])
    % Plot histogram of X and Y
    subplot(2,2,2)
    hold on
    bar(hist(X),'r');
    bar(hist(Y),'b');
    % X and Y plot in line
    subplot(2,1,2)
    hold on
    plot(X(1:100),'r','linewidth',2)
    plot(Y(1:100),'k:','linewidth',2)
    plot([1,100],[0,0],'k','linewidth',2)
    axis([1,100,-0.6,0.6])
    % Covariance and Kurtosis of X and Y
    CXY=X*Y'/1000;
    kX=(sum(X.^4)/length(X))/(sum(X.^2)/length(X))^2-3;
    kY=(sum(Y.^4)/length(Y))/(sum(Y.^2)/length(Y))^2-3;
    disp('     CXY      Kurt(X)      Kurt(Y)')
    disp([CXY,kX,kY])
    
A=[0.7,0.3;0.4,0.6]
input('Anykey to translate X and Y by A : ')
% Translation X and Y into X1 and Y1, or in matrix format Z1
Z1=A*[X;Y];
X1=Z1(1,:);
Y1=Z1(2,:);
clear Z1


% x and y plot in dots
    figure(2)
    subplot(2,2,1)
    hold on
    plot(X1,Y1,'*','linewidth',3)
    plot([-1,1],[0,0],'k','linewidth',2)
    plot([0,0],[-1,1],'k','linewidth',2)
    axis([-1,1,-1,1])
    % Plot histogram of X and Y
    subplot(2,2,2)
    hold on
    bar(hist(X1),'r');
    bar(hist(Y1),'b');    
    % x and y plot in line
    subplot(2,1,2)
    hold on
    plot(X1(1:100),'r','linewidth',2)
    plot(Y1(1:100),'k:','linewidth',2)
    plot([1,100],[0,0],'k','linewidth',2)
    axis([1,100,-0.6,0.6])
    % Covariance and Kurtosis of x and y
    CXY1=X1*Y1'/1000;
    kX1=(sum(X1.^4)/length(X1))/(sum(X1.^2)/length(X1))^2-3;
    kY1=(sum(Y1.^4)/length(Y1))/(sum(Y1.^2)/length(Y1))^2-3;
    disp('    CXY1      Kurt(X1)      Kurt(Y1)')
    disp([CXY1,kX1,kY1])
    
    
input('Anykey to do PCA in X1 and Y1 :')
% PCA get X2 and Y2
Z1(1,:)=X1-mean(X1);
Z1(2,:)=Y1-mean(Y1);
CZ1=Z1*Z1'/1000;
[D,L]=eig(CZ1);
D=flipud(D)
L=fliplr(flipud(L))
Z2=D*Z1;
X2=Z2(1,:);
Y2=Z2(2,:);
clear Z1;


% X2x and Y2 plot in dots
    figure(3)
    subplot(2,2,1)
    hold on
    plot(X2,Y2,'*','linewidth',3)
    plot([-1,1],[0,0],'k','linewidth',2)
    plot([0,0],[-1,1],'k','linewidth',2)
    axis([-1,1,-1,1])
    % Plot histogram of X and Y
    subplot(2,2,2)
    hold on
    bar(hist(X2),'r');
    bar(hist(Y2),'b');    
    % x and y plot in line
    subplot(2,1,2)
    hold on
    plot(X2(1:100),'r','linewidth',2)
    plot(Y2(1:100),'k:','linewidth',2)
    plot([1,100],[0,0],'k','linewidth',2)
    axis([1,100,-0.6,0.6])

    % Covariance and Kurtosis of x and y
    CXY2=X2*Y2'/1000;
    kX2=(sum(X2.^4)/length(X2))/(sum(X2.^2)/length(X2))^2-3;
    kY2=(sum(Y2.^4)/length(Y2))/(sum(Y2.^2)/length(Y2))^2-3;
    disp('     CXY2      Kurt(X2)      Kurt(Y2)')
    disp([CXY2,kX2,kY2])

input('Anykey to recover X and Y from X2 and Y2 by inverse(D*A) : ')
% Recover X and Y from X2 and Y2
disp('Z2=D*(A*Z)=D*A*Z,   Z1 denote [X,Y] and Z2 denote [X2,Y2]')
disp('inverse(D*A)=');
R=inv(D*A);
disp(R);
Z3=R*Z2;
X3=Z3(1,:);
Y3=Z3(2,:);


    % X3 and Y2 plot in dots
    figure(4)
    subplot(2,2,1)
    hold on
    plot(X3,Y3,'*','linewidth',3)
    plot([-1,1],[0,0],'k','linewidth',2)
    plot([0,0],[-1,1],'k','linewidth',2)
    axis([-1,1,-1,1])
    % Plot histogram of X and Y
    subplot(2,2,2)
    hold on
    bar(hist(X3),'r');
    bar(hist(Y3),'b');    
    % X3 and Y3 plot in line
    subplot(2,1,2)
    hold on
    plot(X3(1:100),'r','linewidth',2)
    plot(Y3(1:100),'k:','linewidth',2)
    plot([1,100],[0,0],'k','linewidth',2)
    axis([1,100,-0.6,0.6])

    % Covariance and Kurtosis of x and y
    CXY3=X3*Y3'/1000;
    kX3=(sum(X3.^4)/length(X3))/(sum(X3.^2)/length(X3))^2-3;
    kY3=(sum(Y3.^4)/length(Y3))/(sum(Y3.^2)/length(Y3))^2-3;
    disp('     CXY2      Kurt(X2)      Kurt(Y2)')
    disp([CXY3,kX3,kY3])

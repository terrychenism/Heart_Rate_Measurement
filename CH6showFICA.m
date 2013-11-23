%ShowFastICA
close all
clear all

N=1000;

X=rand(2,N);

% Cancel mean
        figure(1)
        subplot(2,1,1)
        hold on
            plot(X(1,1:50),'r','linewidth',3)
            plot(X(2,1:50),'b-.','linewidth',3)
        subplot(2,3,4)
            [H,A]=hist(X(1,:));
            bar(A,H,'r')
            dA=(A(2)-A(1))/1.9;
            axis([min(A)-dA,max(A)+dA,0,max(H)+10])
        subplot(2,3,5)
            [H,A]=hist(X(2,:));
            bar(A,H,'b')
            dA=(A(2)-A(1))/1.9;
            axis([min(A)-dA,max(A)+dA,0,max(H)+10])
        subplot(2,3,6)
            plot(X(1,:),X(2,:),'*','linewidth',3)
            axis([min(X(1,:))-0.2,max(X(1,:))+0.2,min(X(2,:))-0.2,max(X(2,:))+0.2])
            
% Mixture
A=[0.6,0.4;0.3,0.7];
X=A*X;
        figure(2)
        subplot(2,1,1)
        hold on
            plot(X(1,1:50),'r','linewidth',3)
            plot(X(2,1:50),'b-.','linewidth',3)
        subplot(2,3,4)
            [H,A]=hist(X(1,:));
            bar(A,H,'r')
            dA=(A(2)-A(1))/1.9;
            axis([min(A)-dA,max(A)+dA,0,max(H)+10])
        subplot(2,3,5)
            [H,A]=hist(X(2,:));
            bar(A,H,'b')
            dA=(A(2)-A(1))/1.9;
            axis([min(A)-dA,max(A)+dA,0,max(H)+10]);
        subplot(2,3,6)
            plot(X(1,:),X(2,:),'*','linewidth',3)
            axis([min(X(1,:))-0.1,max(X(1,:))+0.1,min(X(2,:))-0.1,max(X(2,:))+0.1])
           

% Cancel and Whitening

    X(1,:)=X(1,:)-mean(X(1,:));
    X(2,:)=X(2,:)-mean(X(2,:));
    CXX=X*X'/N;
    [D,L]=eig(CXX);
    X=sqrt(inv(L))*D*X;
        figure(3)
        subplot(2,1,1)
        hold on
            plot(X(1,1:50),'r','linewidth',3)
            plot(X(2,1:50),'b-.','linewidth',3)
        subplot(2,3,4)
            [H,A]=hist(X(1,:));
            bar(A,H,'r')
            dA=(A(2)-A(1))/1.9;
            axis([min(A)-dA,max(A)+dA,0,max(H)+10])
        subplot(2,3,5)
            [H,A]=hist(X(2,:));
            bar(A,H,'b')
            dA=(A(2)-A(1))/1.9;
            axis([min(A)-dA,max(A)+dA,0,max(H)+10])
        subplot(2,3,6)
            plot(X(1,:),X(2,:),'*','linewidth',3)
            axis([min(X(1,:))-0.2,max(X(1,:))+0.2,min(X(2,:))-0.2,max(X(2,:))+0.2])
% FastIAC
    w1=rand(1,2);
    w1=w1/sqrt(sum(w1.^2));
    for k=1:10,
        y1=w1*X;
        y1=y1-mean(y1);
        k4=sum(y1.^4)/length(y1);
        k2=sum(y1.^2)/length(y1);
        kt1(k)=k4/k2-3;
        dw=X*diag((w1*X).^3);
        w1=sum(dw')/length(dw)-3*w1;
        w1=w1/sqrt(sum(w1.^2));
    end
    
    w2=w1;
    while (w2-w1)*(w2-w1)'<0.1,
        w2=rand(1,2);
        w2=w2/sqrt(sum(w2.^2));
        for k=1:10,
            y2=w2*X;
            y2=y2-mean(y2);
            k4=sum(y2.^4)/length(y2);
            k2=sum(y2.^2)/length(y2);
            kt2(k)=k4/k2-3;
            dw=X*diag((w2*X).^3);
            w2=sum(dw')/length(dw)-3*w2;
            w2=w2/sqrt(sum(w2.^2));
        end
    end
    y1=(y1-mean(y1))/sqrt(sum(y1.^2)/length(y1));
    y2=(y2-mean(y2))/sqrt(sum(y2.^2)/length(y2));
    
        figure(4)
        subplot(2,1,1)
        hold on
            plot(y1(1:50),'r','linewidth',3)
            plot(y2(1:50),'b-.','linewidth',3)
        subplot(2,3,4)
            [H,A]=hist(y1);
            bar(A,H,'r')
            dA=(A(2)-A(1))/1.9;
            axis([min(A)-dA,max(A)+dA,0,max(H)+10])
         subplot(2,3,5)
            [H,A]=hist(y2);
            bar(A,H,'b')
            dA=(A(2)-A(1))/1.9;
            axis([min(A)-dA,max(A)+dA,0,max(H)+10])
         subplot(2,3,6)
            plot(y1,y2,'*','linewidth',3)
            axis([min(y1)-0.2,max(y1)+0.2,min(y2)-0.2,max(y2)+0.2])
            
        figure(5)
        hold on
            plot(kt1,'r','linewidth',3)
            plot(kt2,'b-.','linewidth',3)
            
       W=[w1;w2]

    
    

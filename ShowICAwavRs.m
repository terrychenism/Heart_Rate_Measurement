% ��Ԫ���룬��׼ICA���򿪶�������ļ�����ϣ��÷������뷨�ֳ����ж���Դ
% 
close all
clear all

% ��ȡ�����ź�,X-ԭʼ�ź�(�У��źţ��У�����),f-����Ƶ��, b-���ݾ���
[X(1,:),f,b]=wavread('s01.wav');               % Ů���質
[X(2,:),f,b]=wavread('s02.wav');               % �����������С�������
[X(3,:),f,b]=wavread('s03.wav');               % ��������ͯ��
%[X(4,:),f,b]=wavread('s04.wav');               % �����
%[X(5,:),f,b]=wavread('s05.wav');               % ��������
%[X(6,:),f,b]=wavread('s06.wav');               % �����質
%[X(7,:),f,b]=wavread('s07.wav');               % Ů������Ԥ��
%[X(8,:),f,b]=wavread('s08.wav');               % ͯ���質

%>>STEP 1 : ԭʼ�źŹ�һ��;
            % 1-1.��һ����ȡ��ֵ�������λ����
            % 1-2.��ʾÿ���źŲ��������Ͷ�
            % 1-3 ѡ�񲥷�
            
    %>1-1.��һ����ȥ��ֵ�������λ��
        M=size(X);              % ԭʼ���ݹ�ģ
        for m=1:M(1)            % ��һ��
            X(m,:)=(X(m,:)-mean(X(m,:)))/std(X(m,:));
        end

    %>1-2.��ʾԭʼ�ź�(fig.1)�������Ͷ� 
        figure(1)
        set(1,'position',[10,10,800,600]);
        whitebg(1,'k');
        for m=1:M(1)
            subplot(M(1),1,m)
            plot(X(m,:),'r','linewidth',1);
            if m==1,
                title('ԭʼ�ź�')
            end
            axis([1,M(2),min(X(m,:)),max(X(m,:))]);
            kurtXo(m)=(sum(X(m,:).^4)/M(2))/(sum(X(m,:).^2)/M(2))^2-3;
        end
        
    %>1-3.ѡ�񲥷�
        disp('  ');
        while (m~=0)&&(m<=M(1)),
            m=input('��������Բ���ԭʼ����, ���� 0 ��ת :');
            if isempty(m),
                break
            else
                wavplay(X(m,:)/3,f);
            end
        end
        %Xorg=X;  % remain the original signals to Xorg
%>> Conclusion after step 1:  ԭʼ�źŲ�����. X �� Xorg Ϊ��һ������. M(1) �� M(2) Ϊ��һ�����ݹ�ģ;
    

%>>STEP 2: ���ԭʼ���ݵõ��۲�����
            % 2-1.���ɻ�Ͼ���;
            % 2-2.��ʾ�������;
            % 2-3.ѡ�񲥷� 
            
    %>2-1.���ɻ�Ͼ���: A,��Ͼ���; X,��ϵõ��Ĺ۲����; N,����źŵĹ�ģ.        
        A=1/7+(6/7)*rand(M(1),M(1));                % ��Ͼ���, �����С���ȱ�Ϊ 7
        X=A*X;                                      % X ����źţ������۲��ź�
        M=size(X);                                  % �۲��źŵĹ�ģ
        
    %>2-2.��ʾ�۲��ź�               
        figure(2)
        set(2,'position',[20,20,800,600]);
        whitebg(2,'k');
        for m=1:M(1)
            subplot(M(1),1,m)
            plot(X(m,:),'r','linewidth',1)
            axis([1,M(2),min(X(m,:)),max(X(m,:))]);
            kurtXm(m)=(sum(X(m,:).^4)/M(2))/(sum(X(m,:).^2)/M(2))^2-3;
            if m==1,
                title('�۲�(���)�ź�')
            end
        end
        
    %>2-3.ѡ�񲥷Ż���ź� 
        disp('  ');
        while (m~=0)&&(m<=M(1)),
            m=input('��������Բ��Ź۲��ź�, ���� 0 ��ת :');
            if isempty(m),
                break
            else
                wavplay(X(m,:)/3,f);
            end
        end
%>> Conclusion after step 2:  ԭʼ�ź��ޱ���. X Ϊ�۲��ź�. N(1) & N(2) Ϊ���ģ�� A Ϊ��Ͼ���;
        
        
%>>STEP 3:  �������г�ȡ�����������������м�ȥ����������ѹ��һ��ά���õ�������, ��������N(2)-1�Σ����ʣ��һ����������
        %3-1.��X�г�ȡһ����������Y. X ��ǰ�μ�ȥ����������ʣ������ݣ� 
        %3-2.�ӵ�ǰ����X�м�ȥ��ǰ���������õ���X;
        %3-3.�Լ�ȥ��������ʣ���X��PCA����һά��
        %3-4.�Խ�ά���Xȥ��ֵ�������һ��,ת3-1;
        %3-5.��ʾ����Դ fig3.
        %3-6.ѡ�񲥷�   
        
        Ma=M;
        for m=1:M(1)-1
        
            %>3-1. ����ICA�ӳ��򣬴� X �г�ȡһ���������� Y(m,:);
            [y,B,BB,ktY]=SubICA(X,35);
            Y(m,:)=y;                   % �� m ����������
            kurtY(m)=max(ktY);          % �� m �������������Ͷ�
            clear B BB ktY
            
            %>3-2.�ӵ�ǰ���� X �м�ȥ�յõ��Ķ���������������
            %       �Ż����� a ʹĿ�꺯�� sum([X(ma,:)-a*y]^2) �ﵽ��С��
            %       ���� a �ĵ����㷨��a=a+[X(ma,:)-a*y]*y'
            for ma=1:Ma(1)
                a=1;
                ymax=max(y);
                ymin=min(y);
                for k=1:1000
                    da=(X(ma,:)-a*y)*y'/Ma(2);
                    a=a+da;
                    if da<(abs(da)<=ymin/(10000*ymax)),
                        break
                    end
                end
                X(ma,:)=X(ma,:)-a*y;
            end
            clear y
        
            %>3-3.��PCA����ȥ��������������ݽ���һά.
            [D,L]=eig(X*X'/Ma(2));
            D=D';
            X=D(2:Ma(1),:)*X;
            clear D L
            
            %>3-4.��ά������ݹ�һ��(��ȥ����Դ�󷽲�ٵ�λ��)
            Ma=size(X);
            for ma=1:Ma(1)
                X(ma,:)=(X(ma,:)-mean(X(ma,:)))/std(X(ma,:));
            end
            
        end
        % �����ݽ������һ��ʱ����Ϊ���һ������Դ
        kurtY(m+1)=(sum(X.^4)/Ma(2))/(sum(X.^2)/Ma(2))^2-3;
        Y(m+1,:)=X;    
        
        %>3-5.��ʾ��������
        figure(3)
        set(3,'position',[30,30,800,600]);
        whitebg(3,'k');
        for m=1:M(1)
            subplot(M(1),1,m)
            plot(Y(m,:),'r','linewidth',1)
            axis([1,M(2),min(Y(m,:)),max(Y(m,:))]);
            if m==1,
                title('��������')
            end                
        end
        
        %>3-6.ѡ�񲥷Ŷ�������
        disp(' ');
        while (m~=0)&&(m<=M(1)),
            m=input('��������Բ��Ŷ������������� 0 ��ת :');
            if isempty(m),
                break
            else
                wavplay(Y(m,:)/3,f);
            end
        end
        
        clear X n m D L Na a b k
        
disp(' ')       
disp('ԭʼ�źŵ��Ͷ� :')
disp(kurtXo)
disp('  ')
disp('����źŵ��Ͷ� :')
disp(kurtXm)
disp('  ')
disp('����Դ���Ͷ� :')
disp(kurtY)
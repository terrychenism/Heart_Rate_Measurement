% 多元分离，标准ICA，打开多个声音文件，混合，用反复分离法分出所有独立源
% 
close all
clear all

% 读取声音信号,X-原始信号(行：信号，列：样本),f-采样频率, b-数据精度
[X(1,:),f,b]=wavread('s01.wav');               % 女声歌唱
[X(2,:),f,b]=wavread('s02.wav');               % 男声配乐朗诵——大理
[X(3,:),f,b]=wavread('s03.wav');               % 男声配乐童话
%[X(4,:),f,b]=wavread('s04.wav');               % 打击乐
%[X(5,:),f,b]=wavread('s05.wav');               % 男声朗诵
%[X(6,:),f,b]=wavread('s06.wav');               % 男声歌唱
%[X(7,:),f,b]=wavread('s07.wav');               % 女声天气预报
%[X(8,:),f,b]=wavread('s08.wav');               % 童声歌唱

%>>STEP 1 : 原始信号归一化;
            % 1-1.归一化：取均值、均方差单位化；
            % 1-2.显示每个信号并计算其峭度
            % 1-3 选择播放
            
    %>1-1.归一化：去均值、均方差单位化
        M=size(X);              % 原始数据规模
        for m=1:M(1)            % 归一化
            X(m,:)=(X(m,:)-mean(X(m,:)))/std(X(m,:));
        end

    %>1-2.显示原始信号(fig.1)并计算峭度 
        figure(1)
        set(1,'position',[10,10,800,600]);
        whitebg(1,'k');
        for m=1:M(1)
            subplot(M(1),1,m)
            plot(X(m,:),'r','linewidth',1);
            if m==1,
                title('原始信号')
            end
            axis([1,M(2),min(X(m,:)),max(X(m,:))]);
            kurtXo(m)=(sum(X(m,:).^4)/M(2))/(sum(X(m,:).^2)/M(2))^2-3;
        end
        
    %>1-3.选择播放
        disp('  ');
        while (m~=0)&&(m<=M(1)),
            m=input('输入序号以播放原始声音, 输入 0 跳转 :');
            if isempty(m),
                break
            else
                wavplay(X(m,:)/3,f);
            end
        end
        %Xorg=X;  % remain the original signals to Xorg
%>> Conclusion after step 1:  原始信号不保留. X 和 Xorg 为归一化数据. M(1) 和 M(2) 为归一化数据规模;
    

%>>STEP 2: 混合原始数据得到观测数据
            % 2-1.生成混合矩阵;
            % 2-2.显示混合数据;
            % 2-3.选择播放 
            
    %>2-1.生成混合矩阵: A,混合矩阵; X,混合得到的观测矩阵; N,混合信号的规模.        
        A=1/7+(6/7)*rand(M(1),M(1));                % 混合矩阵, 最大最小幅度比为 7
        X=A*X;                                      % X 混合信号，即，观测信号
        M=size(X);                                  % 观测信号的规模
        
    %>2-2.显示观测信号               
        figure(2)
        set(2,'position',[20,20,800,600]);
        whitebg(2,'k');
        for m=1:M(1)
            subplot(M(1),1,m)
            plot(X(m,:),'r','linewidth',1)
            axis([1,M(2),min(X(m,:)),max(X(m,:))]);
            kurtXm(m)=(sum(X(m,:).^4)/M(2))/(sum(X(m,:).^2)/M(2))^2-3;
            if m==1,
                title('观测(混合)信号')
            end
        end
        
    %>2-3.选择播放混合信号 
        disp('  ');
        while (m~=0)&&(m<=M(1)),
            m=input('输入序号以播放观测信号, 输入 0 跳转 :');
            if isempty(m),
                break
            else
                wavplay(X(m,:)/3,f);
            end
        end
%>> Conclusion after step 2:  原始信号无保留. X 为观测信号. N(1) & N(2) 为其规模， A 为混合矩阵;
        
        
%>>STEP 3:  从数据中抽取独立分量，从数据中减去独立分量并压缩一次维数得到新数据, 反复进行N(2)-1次，最后剩余一个独立分量
        %3-1.从X中抽取一个独立分量Y. X 是前次减去独立分量后剩余的数据； 
        %3-2.从当前数据X中减去当前独立分量得到新X;
        %3-3.对减去独立分量剩余的X用PCA降低一维；
        %3-4.对降维后的X去均值、方差归一化,转3-1;
        %3-5.显示独立源 fig3.
        %3-6.选择播放   
        
        Ma=M;
        for m=1:M(1)-1
        
            %>3-1. 调用ICA子程序，从 X 中抽取一个独立分量 Y(m,:);
            [y,B,BB,ktY]=SubICA(X,35);
            Y(m,:)=y;                   % 第 m 个独立分量
            kurtY(m)=max(ktY);          % 第 m 个独立分量的峭度
            clear B BB ktY
            
            %>3-2.从当前数据 X 中减去刚得到的独立分量，方法：
            %       优化参数 a 使目标函数 sum([X(ma,:)-a*y]^2) 达到最小；
            %       参数 a 的迭代算法：a=a+[X(ma,:)-a*y]*y'
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
        
            %>3-3.用PCA将减去独立分量后的数据降低一维.
            [D,L]=eig(X*X'/Ma(2));
            D=D';
            X=D(2:Ma(1),:)*X;
            clear D L
            
            %>3-4.降维后的数据归一化(减去独立源后方差不再单位化)
            Ma=size(X);
            for ma=1:Ma(1)
                X(ma,:)=(X(ma,:)-mean(X(ma,:)))/std(X(ma,:));
            end
            
        end
        % 当数据降至最后一行时，此为最后一个独立源
        kurtY(m+1)=(sum(X.^4)/Ma(2))/(sum(X.^2)/Ma(2))^2-3;
        Y(m+1,:)=X;    
        
        %>3-5.显示独立分量
        figure(3)
        set(3,'position',[30,30,800,600]);
        whitebg(3,'k');
        for m=1:M(1)
            subplot(M(1),1,m)
            plot(Y(m,:),'r','linewidth',1)
            axis([1,M(2),min(Y(m,:)),max(Y(m,:))]);
            if m==1,
                title('独立分量')
            end                
        end
        
        %>3-6.选择播放独立分量
        disp(' ');
        while (m~=0)&&(m<=M(1)),
            m=input('输入序号以播放独立分量，输入 0 跳转 :');
            if isempty(m),
                break
            else
                wavplay(Y(m,:)/3,f);
            end
        end
        
        clear X n m D L Na a b k
        
disp(' ')       
disp('原始信号的峭度 :')
disp(kurtXo)
disp('  ')
disp('混合信号的峭度 :')
disp(kurtXm)
disp('  ')
disp('独立源的峭度 :')
disp(kurtY)
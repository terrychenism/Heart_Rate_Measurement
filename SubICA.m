% Sub-routing ICA with grade method
    % input data:
        % X observed signals, rows be variables and collums be samples;
        % n loop times when finding optimal B, a vector as y=B*X is independent component;
    % output data
        % y the found independent compnent;
        % B the vector in the orthogonal data space, the independent component is the map of X on B;
        % DB the increcement of B when evolued in optimizing procedure
        % BB the process series of B evolued;
        % ktY the process series of kurtosis evolued
        % X the orthogonalized input data
        
    % procedures
        %1.orthogonalize signals by PCA
        %2.normalize the orthogonal signals
        %3.generate a normalized vector B
        %4.calculing dB and naximum kurt(y) by B, get a B by maximum kurt(y)
    
    function [y,B,BB,ktY]=SubICA(X,loops)
    
        N=size(X);
        
        %1. orthogonal the signals by PCA
        CXX=X*X'/N(2);          % covariance of observed(mixed) signal
        [D,L]=eig(CXX);         % eigen value decompose
        X=D'*X;                 % orthogonal transform
        clear CXX D L
        
        %2.normalizing the orthogonal the signals
        for n=1:N(1)
            X(n,:)=(X(n,:)-mean(X(n,:)))/std(X(n,:));
        end

            
        %3.generate a normalized vector B
        B=randn(1,N(1));            %
        B=B/sqrt(sum(B.^2));        % Generalizing
        y=B*X;
        
        %4.calculing dB and naximum kurt(y) by B, get a B with maximum kurt(y)
        k=0;
        while k<loops
            k=k+1;
            kty=sum(y.^4)/N(2)-3;
            ktY(k)=kty;
            % calculate dB=E{X(y^3)}    
                dB=X*(y.^3)'/N(2);      
                B=B+0.1*kty*dB';
            % normalizing B be unity vector
            B=B/sqrt(sum(B.^2));
            % record dB and B series for drawn
            BB(k,:)=B;
            y=B*X;
            if (k>3)&&(abs(ktY(k)-ktY(k-1))<10^-8)
                break
            end

        end

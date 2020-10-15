function [A,L,Ac,xy] = get_network( options)
%This function generates a network and returns the following parameters:

%         A:            The matrix that defines the nodes connection
%         L:            Laplacian matrix
%         An:           Algebraic connectivity
%         xy:           Nodes coordinate
%
% Inputs:
%         options.N:                    Number of nodes                     (default=10)
%         options.type:                 Topology type: binomial or radial   (default='radial')
%         options.parameter:            Probablility for binomila topology and radius for (default=0.2)
%                                       radial topology
%                                       radial topology2
%         options.minAllowableDistance  (default=0.5)
%         options.sigma_net             Position variance (default=1)



%% Network generation
options=initoptions('get_network',options);

N=options.N;
type=options.type;
parameter=options.parameter;

A=zeros(N);       % topology matrix
L=zeros(N);       % Laplacian matrix






Ac=0;
while(Ac<1e-4)
    
    % Generating random coordinate
    
    
    %x=round(sqrt(options.sigma_net)*randn(100*N,1)) ;
    %y=round(sqrt(options.sigma_net)*randn(100*N,1)) ;
    
    minAllowableDistance = options.minAllowableDistance ;
    
    % Initialize first point.
    keeperX = sqrt(options.sigma_net)*randn(1);
    keeperY = sqrt(options.sigma_net)*randn(1);
    % Try dropping down more points.
    counter = 2;
    %k=2;
    while counter<=N
        
        % Get a trial point.
        thisX = round(sqrt(options.sigma_net)*randn(1));
        thisY = round(sqrt(options.sigma_net)*randn(1));
        % See how far is away from existing keeper points.
        distances = sqrt((thisX-keeperX).^2 + (thisY - keeperY).^2);
        minDistance = min(distances);
        if minDistance >= minAllowableDistance
            keeperX(counter) = thisX;
            keeperY(counter) = thisY;
            counter = counter + 1;
        end
        %k=k+1;
    end
    
    
    x=keeperX';
    y=keeperY';
    xy=[x y];
    
    
    switch type
        case 'binomial'
            p=parameter;
            % Linking the nodes randomly
            
            for i=1:length(x)
                A(i,i)=1;
                for j=i+1:length(y)
                    pn=rand(1,1);
                    if (pn>p)
                        A(i,j)=1;
                        A(j,i)=1;
                    end
                end
            end
            
            
        case 'radial'
            r=parameter;
            for i=1:length(x)
                
                for j=1:length(y)
                    rn=norm([x(j)-x(i) y(j)-y(i)],2);
                    A(i,j)=(rn<=r);
                end
            end
            
        case 'gaussian'
            
            for i=1:length(x)
                A(i,i)=1;
                for j=i+1:length(y)
                    sigma=1;
                    mu=0.1;
                    p=parameter;
                    rn=norm([x(j)-x(i) y(j)-y(i)],2);
                    pn=(1/(sigma*sqrt(2*pi)))*exp(-(1/2)*(rn-mu)^2/sigma^2);
                    
                    if (pn>p)
                        A(i,j)=1;
                        A(j,i)=1;
                    end
                end
            end
            
    end
    
    %% Algebric Connectivity
    
    for i=1:N
        L(i,i)=max(0,(sum(A(:,i))-1));
        for j=i+1:N
            L(i,j)=-1*A(i,j);
            L(j,i)=-1*A(j,i);
        end
    end
    
    lambda=eig(L);
    Ac=lambda(2);
    
end


end


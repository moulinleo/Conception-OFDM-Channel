function [A_LOS,A_NLOS] = Beamforming(X,Y,Z,theta,phi,lambda,h_LOS,h_NLOS)
% Use the beamforming function to compute the direction of arrival of each
% wave. By deducind the set of incident plane waves, a new channel model
% can be build.


    %% Compute all positions of the antenna

    r=zeros(length(X)*length(Y)*length(Z),3);
    i=1;
    for x=1:length(X)
        for y=1:length(Y)
            for z=1:length(Z)  
                r(i,:) = 0.02*[X(x) Y(y) Z(z)];
                i=i+1;
            end
        end
    end

    %% Compute the 3D beamformer function

    B = zeros(length(r),length(theta),length(phi));
    for pos=1:length(r)
        for t=1:length(theta)
            for p=1:length(phi) 
                beta =(2*pi/lambda) * [sin(theta(t))*cos(phi(p)) sin(theta(t))*sin(phi(p)) cos(theta(t))];
                B(pos,t,p)=exp(1j*(beta*r(pos,:).'));
            end
        end     
    end

    %% Compute the amplitude of each incident wave in the direction (theta,phi)

    [~,NbTaps]=size(h_LOS);
    A_LOS =zeros(NbTaps,length(theta),length(phi));
    A_NLOS =zeros(NbTaps,length(theta),length(phi));
    for n=1:NbTaps
        for t=1:length(theta)
            for p=1:length(phi)
                Den=sum(abs(B(:,t,p)).^2,1); 
                A_LOS(n,t,p) = sum((h_LOS(:,n).*conj((B(:,t,p)))),1)./Den;
                A_NLOS(n,t,p) = sum((h_NLOS(:,n).*conj((B(:,t,p)))),1)./Den;
            end
        end
    end
end


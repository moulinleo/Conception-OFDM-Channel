function [H_new_LOS, H_new_NLOS] = Compute_Channel(A_LOS,A_NLOS,theta,phi,X,Y,Z,lambda)
% Build new channel model from the set of incident waves.

    H_new_LOS = zeros(10,10,10,10);
    H_new_NLOS = zeros(10,10,10,10);
    phase = 2*pi*rand(1000,1);
    save_LOS = 0;
    save_NLOS = 0;
    i =1;
    for n = 1:10
        for x = 1:10
        for y = 1:10
        for z = 1:10
            for t=1:length(theta)
                for p=1:length(phi)
                   beta = (2*pi/lambda)*(0.02*X(x).*sin(theta(t))*cos(phi(p))+0.02*Y(y).*sin(theta(t))*sin(phi(p))+0.02*Z(z).*cos(theta(t)));
                   beta1 = exp(1j*(phase(i)-beta));
                   save_LOS = save_LOS + A_LOS(n,t,p).*beta1;
                   save_NLOS = save_NLOS + A_NLOS(n,t,p).*beta1;
                end
            end
            H_new_LOS(x,y,z,n)=save_LOS;    
            H_new_NLOS(x,y,z,n)=save_NLOS;    
            i = i+1;
        end
        end
        end
        i=1;
    end

end


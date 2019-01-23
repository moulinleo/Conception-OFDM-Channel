function index = findLocalMaxima(A, treshold)
%threshold: -50dB
% Detects local maxima of 2-dimensional matrix. The returned 2-line
% matrix 'index' containts the indexes of the local maxima on each column.
% To be detected, a maxima should be higher than the minimum value defined
% in 'treshold'.


[dimLine,dimCol] = size(A);
cptLocalMax = 0;
line = 1; col = 1;

for k = 2:dimLine-1
    for m = 2:dimCol-1
        if (A(k,m)>treshold)         
            % condition pour les maxima locaux sur une matrice 3D
            if( (A(k,m)>A(k-1,m-1)) &&  (A(k,m)>A(k-1,m)) && (A(k,m)>A(k-1,m+1)) && (A(k,m)>A(k,m-1)) && (A(k,m)>A(k,m+1)) && (A(k,m)>A(k+1,m-1)) && (A(k,m)>A(k+1,m)) && (A(k,m)>A(k+1,m+1)) )
                % disp('je suis la');
                % local maxima detected!
                cptLocalMax = cptLocalMax+1;
                line(cptLocalMax) = k;
                col(cptLocalMax) = m;     
            end  
        end
    end
end

index(1,:) = line;
index(2,:) = col;




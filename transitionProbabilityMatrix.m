function P = transitionProbabilityMatrix( G,  A,  domA, size_Dom_a, w, numNeighborVertices, numNeighborAttributeVertices )
    N = size(G, 1);
    m = size(A, 2);
    
    P = zeros(N+sum(size_Dom_a));
    
    for i = 1:N
        p_i = 1/(sum(w) + numNeighborVertices(i));
        
        % prijelaz od vi do vj kroz strukturalni vrh
        % G = matrica koja sadrzi info o povezanosti = sadrzi jedinice na
        % mjestu gdje su povezani pa se zato mnozi sa G(??)
        P(i, 1:N) = G(i, 1:N) * p_i;
        
        % prijelaz kroz atributni vrh...v_i do v_jk i v_ik do v_j
        d = 0;
        for j = 1:m
            for k = 1:size_Dom_a(j)
                % ako je vrijednost na (i,j)-tom mjestu u matrici A jednaka
                % vrijednosti Dom_a_d+k (??)
                if A(i,j) == domA(d+k)
                    % N+d+k = N jer suatributni vrhovi nakon N vrhova
                    % grafa, k jer je od k-tog atributa, a d je od velicine
                    % domene j-tog atributa
                    P(i, N+d+k) = w(j) * p_i;
                    P(N+d+k,i) = 1/numNeighborAttributeVertices(j,k);
                end
            end
            d = d + size_Dom_a(j);
        end
    end
end
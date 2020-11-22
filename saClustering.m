function [ cluster, centroids ] = saCluster( G, A, size_dom_a, domA, omega, iteration, z, l, c, k )
% 1. inicijalizirati w1,...wm na 1 i fiksirati w0=1
% 2. izracunati matricu ujedinjene slučajne šetnje R
% 3. izabrati k inicijalnih centroida s najvecom gustocom
% 4. ponavljati dok objektna funkcija konvergira
%   5. dodijeliti svaki vrh vi klasteru C* sa centroidom c* gdje je
%       c* = argmax d(vi,cj)
%   6. update centroide sa najcentralinijm tockama u klasteru pomocu
%       R(vi_potez, vj) = 1/abs(vi) suma_R(vk,vj)_po_svim_vk_iz_Vi
%       ci^(t+1) = argmin ||R(vj)-R(vi_crtice)||, 
%       R(vi_crtica) = prosjecni vektor slucajne setnje za klaster Vi
%   7. update tezine w1,...,wm prema fromuli (***)
%   8. ponovno izracunati R
% vratiti k klastera V1,...,Vk

% G = matrica koja sadzi podatke o povezanosti vrhova u atributnom grafu
% A = matrica koja je sastavljena od vektora vrijednsoti na atributima
% a1,...,am za pojedini vrh
% size_dom_a = [n1,...,nm] = skup abs(Dom_ai)
% domA = [Dom_a1,...,Dom_am]  = skup domena Dom_ai = {a_i1,...,a_ini} =
% sve vrijednosti atributa

% m = broj atributa
% N = broj vrhova
% numNeighborVertices = broj susjednih koji su strukturalni vrhovi
% numNeighborAttributeVertices = broj susjeda koj su atributni vrhovi
m = size(A, 2);
N = size(G, 1);
%  ne razumijem sto napravi sum(G') !!!!!
numNeighborVertices =  sum(G');
numNeighborAttributeVertices = zeros(m, max(size_dom_a)); 

f_density = zeros(1, N);
cluster = zeros(1, k);

% 1. inicijalizacija w1,...,wm na 1
w = ones(m,1);
 
% 2. racunanje matrice R ujednjene slučajne setnje
% potrebni su matrica vjerojatnosti prelaska za atributni graf P, startna
% vjerojatnost c i maksimalna duljina slucajne setnje l
% skuziti sto je z !!!!!!
if z == 1 
    P = transitionProbabiliyMatrix(G, A, domA, size_dom_a, w, numNeighborVertices, numNeighborAttributeVertices );
else
    P = transitionProbabiliyMatrix(G, A, domA, size_dom_a, zeros(1,m), numNeighborVertices, numNeighborAttributeVertices );
end  
R = unifiedRandomWalkMatrix(l, c, P);

% 3. odabir k inicijalnih centroida s najvecom gustocom
for i = 1:N
    f_density(i) = density(N, R, omega, i);
end

[ f_densitySort, index ] = sort(f_density, 'descend');
% koristimo da bi imali indekse da kasnije lakse pretrazujemo centroide
centroids = index(1:k);

% vote = mjera koja odreduje dijelel li 2 vrha istu vrijednost atributa
% V = matrica koja sadrzi sve vrijednosti mjere vote
for i = 1:k
    V(:, (i-1)*N+1:i*N) = vote(N, A, i);
end

% 4. iteriracije
for it = 1:iteration
    alone = 0;
    
    % 8. ponovno racunanje R, jer smo izvan petlje izracunali vec R pa da
    % radimo s dobrim R treba prije svega to izracunati, sto je zapravo na
    % kraju prethodne iteracije
    if it~=1
        if z == 1 
            P = transitionProbabiliyMatrix(G, A, domA, size_dom_a, w, numNeighborVertices, numNeighborAttributeVertices );
        else
            P = transitionProbabiliyMatrix(G, A, domA, size_dom_a, zeros(1,m), numNeighborVertices, numNeighborAttributeVertices );
        end  
        R = unifiedRandomWalkMatrix(l, c, P);
    end
     
    % 5. dodijeliti  vrhove klasterima
    for i = 1:N
        % ako nije centroid
        % nisam sigurna zasto imaju size(find, 2)
        if size(find(centroids == i)) == 0
            % trazmimo maksimalnu setnju jer sto je veca setnja to je veca
            % mogucnost da su u istom klasteru
            maxWalk = -1;
            % idemo po centroidima, tj indeksima
            for j = 1:k
                if R(i, centroids(j)) > maxWalk
                    maxWalk = R(i, centroids(j));
                    indexMaxWalk = centroids(j);
                end
            end
            cluster(i) = indexMaxWalk;
        % inace je centoid pa je sam sebi centroid i sam sa sobom je u
        % klasteru
        else
            cluster(i) = i;
        end
    end
    
    % 6. update cenrtoida
    % racunamo prosjecne setnje u klasterima
    averageWalk = zeros(k, N);
    howManyInCluster = zeros(k, 1);
    % po svim centroidima
    for i = 1:k
        % po svim vrhovima
        for v = 1:N
            if cluster(v) == centroids(i)
                howManyInCluster(i) = howManyInCluster(i) + 1;
                averageWalk(i,:) = averageWalk(i,:) + R(v,1:N);
            end
        end
        
        if howManyInCluster(i) == 1
            % sam je sebi centroid i sam je u klasteru, zapamtimo ga 
            alone = i;
        else
            averageWalk(i,:) = averageWalk(i,:) / howManyInCluster(i);
        end
    end
    
    % centroid koji je sam u klasteru stavimo u najblizi klaster 
    % taj slucaj nisam sigurna sto su radili i trebamo li to 
    % !!!!
    
    % novi centroid = prosjecna setnja mu je najbliza vektoru prosjecne setnje (averageWalk) klastera
    for i = 1:k
        min = 2;
        for v = 1:N
            if cluster(v) == centroids(i)
                if norm(R(v, 1:N) - averageWalk(i,:)) < min
                    min = norm(R(v, 1:N) - averageWalk(i,:));
                    minIndex = v;
                end
            end
        end
        centroids(i) = minIndex;
    end
    
    % 7. update tezine w1,...,wm
    sumVote = 0;
    vote_i = zeros(m,1);
    for p = 1:m
        for j = 1:k
            for v = 1:N
                if claster(v) == centroids(j)
                    % brojnik
                    sumVote = sumVote + V(centroids(j), (p-1)*N+i);
                    % nazivnik
                    vote_i(p) = vote_i(p) + V(centroids(j), (p-1)*N+i);
                end
            end
        end
    end
    delta_w = (m*vote_i) / sumVote;
    w = 1/2 * (w + delta_w);

    if z == 1
        f_objective(it) = objectiveFunction (N, R, cluster, centroids);
    end
end

end

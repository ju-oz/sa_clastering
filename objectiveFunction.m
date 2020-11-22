function f_objective = objectiveFunction( N, R, cluster, centroids )
% vraca vrijednost objektivne funkcije u iteraciji u kojoj je pozvana
% N = broj vrhova
% R = vektor slucajnih setnji
% clusters = klasteri u ovoj iteraciji
% centroids = centroidi u ovoj iteraciji

k = max(size(centroids));
f_objective = 0;
clasterCardinalitets = zeros(k,1);

% racunamo vroj vrhova u svakom klasteru
for v = 1:N
    for j = 1:k
        if cluster(v) == centroids(j)
            clasterCardinalitets(j) = clasterCardinalitets(j) + 1;
        end
    end
end

% racunamo objektivnu fju
for i = 1:N
    for j = 1:N
        if cluster(i) == cluster(j)
            clust= find(centroids == cluster(j));
            f_objective = f_objective + R(i,j)/(clasterCardinalitets(clust) * clasterCardinalitets(clust));
        end
    end
end
end
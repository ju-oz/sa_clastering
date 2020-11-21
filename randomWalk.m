function R = randomWalk( l, c, P )
% l = zadana duljina slucajne setnje
% c = iz (0,1), vjerojatnost restarta
% P = matrica vjerojatnosti prijelaza grafa

    R = zeros(size(P,1), size(P,2));
    for gama=1:l
        R = R + c*(1-c)*P^gama;
    end
end
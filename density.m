function f = density( N, R, omega, i )
% N = vroj vrhova
% f = funkcija gustoce
% i = vrh ciju fju gustoce racunamo
% omega = zadani parametar
% R = matrica slucajne setnje

% iz nekog razloga baca error!!!!!!!!

    f = 0;
    
    for j = 1:N
      f = f + (1-exp(-R(i,j)^2/(2*omega^2)));  
    end
end

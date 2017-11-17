function Sbus = ComputeSbus(Lineas, Gen, Cargas)

%   MATPOWER
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2011 by Power System Engineering Research Center (PSERC)

on = find(Gen(:,1) > 0);      %% which generators are on?
gbus = Gen(:,1);                %% what buses are they at?

%% form net complex bus power injection vector
nb = max(max(Lineas(:,1)),max(Lineas(:,2)));
ngon = size(on, 1);
nload = Cargas(:,1);
Cg = sparse(gbus, (1:ngon)', ones(ngon, 1), nb, ngon);  %% connection matrix
                                                        % element i, j is 1 if
                                                        %% gen on(j) at bus i is ON
Sbus =  ( Cg * (Gen(gbus, 2) + 1j * Gen(gbus, 4) ) ... %% power injected by generators
           - (Cargas(:, 3) + 1j * Cargas(:, 4)) );
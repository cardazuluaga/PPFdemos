function Ybus = MakeYbus(Lineas)

%   MATPOWER
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2011 by Power System Engineering Research Center (PSERC)

nb = max(max(Lineas(:,1)),max(Lineas(:,2)));          %% number of buses
nl = size(Lineas, 1);                                 %% number of lines

stat = ones(nl,1);
Ys = stat ./ (Lineas(:,3) + 1j * Lineas(:,4));  %% series admittance
Bc = stat .* Lineas(:,5);                       %% line charging susceptance
tap = ones(nl, 1);                              %% default tap ratio = 1
i = find(Lineas(:, 6));                         %% indices of non-zero tap ratios
tap(i) = Lineas(i, 6);                          %% assign non-zero tap ratios
tap = tap .* exp(1j*pi/180 * Lineas(:, 7)); %% add phase shifters
Ytt = Ys + 1j*Bc/2;
Yff = Ytt ./ (tap .* conj(tap));
Yft = - Ys ./ conj(tap);
Ytf = - Ys ./ tap;

%% i.e. Ysh = Psh + j Qsh, so ...
Ysh = zeros(nb,1); %% vector of shunt admittances

%% build connection matrices
f = Lineas(:, 1);                           %% list of "from" buses
t = Lineas(:, 2);                           %% list of "to" buses
Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses

%% build Yf and Yt such that Yf * V is the vector of complex branch currents injected
%% at each branch's "from" bus, and Yt is the same for the "to" bus end
i = [1:nl; 1:nl]';                              %% double set of row indices
Yf = sparse(i, [f; t], [Ytt; Ytf], nl, nb);
Yt = sparse(i, [f; t], [Yft; Yff], nl, nb);

% Yf = sparse(i, [f; t], [Yff; Yft], nl, nb);
% Yt = sparse(i, [f; t], [Ytf; Ytt], nl, nb);

% Yf = spdiags(Yff, 0, nl, nl) * Cf + spdiags(Yft, 0, nl, nl) * Ct;
% Yt = spdiags(Ytf, 0, nl, nl) * Cf + spdiags(Ytt, 0, nl, nl) * Ct;

%% build Ybus
Ybus = Cf' * Yf + Ct' * Yt + ...                %% branch admittances
               sparse(1:nb, 1:nb, Ysh, nb, nb);        %% shunt admittance

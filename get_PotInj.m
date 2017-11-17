function [PQ] = get_PotInj(x,Ybus,ng)

N = size(Ybus,1);
Vnt = x(N+1:end)';
TheT = x(1:N)';
Vbus = Vnt.*(exp(sqrt(-1).*TheT));
Ibus = Ybus * Vbus;
Sbus = Vbus .* conj(Ibus);
Pbus = real(Sbus);
Qbus = imag(Sbus);
PQ = [Pbus(2:end);Qbus(ng+2:end)];

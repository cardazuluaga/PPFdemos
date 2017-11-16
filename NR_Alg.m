function Res = NR_Alg(Datos)

nb = max(max(Datos.Lineas(:,1)),max(Datos.Lineas(:,2)));
V0  = ones(nb,1) .* exp(sqrt(-1) * pi/180 * zeros(nb,1));
gbus = Datos.Gen(:,1);
ref = find(Datos.Cargas(:,2) == 3);
pv = find(Datos.Cargas(:,2) == 2);
pq = find(Datos.Cargas(:,2) == 1);
mpopt = mpoption;
V0(gbus) = Datos.Gen(gbus, 3);
Ybus = HacerYbus(Datos.Lineas);
Sbus = CalcularSbus(Datos.Lineas, Datos.Gen, Datos.Cargas);
[V, success, iterations] = newtonpf_Alg(Ybus, Sbus, V0, ref, pv, pq, mpopt);
Res.V = V;
Res.success = success;
Res.iter = iterations;
Res.YBUS = Ybus;
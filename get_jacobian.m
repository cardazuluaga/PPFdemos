function [dx,J] = get_jacobian(V,PQ,B,Ybus,pv,pq)

F = B - PQ;
    
[dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);
    
j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
j12 = real(dSbus_dVm([pv; pq], pq));
j21 = imag(dSbus_dVa(pq, [pv; pq]));
j22 = imag(dSbus_dVm(pq, pq));

J = [   j11 j12;
        j21 j22;    ];

dx = (J \ F);   
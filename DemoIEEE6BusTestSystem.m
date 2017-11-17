clc
clear
close all

s = RandStream('mt19937ar', 'Seed', 6.0e6); 
RandStream.setGlobalStream(s);

Sbase = 100;

Datos.Lineas = [1	2	0.1     0.2    0.0200   0     0    
                1	4	0.05	0.2    0.0200   0     0
                1	5	0.08	0.3    0.0300   0     0
                2   3   0.05    0.25   0.0300   0     0
                2   4   0.05    0.1    0.0100   0     0
                2   5   0.1     0.3    0.0200   0     0
                2   6   0.07    0.2    0.0250   0     0
                3   5   0.12    0.26   0.0250   0     0
                3   6   0.02    0.1    0.0100   0     0
                4   5   0.2     0.4    0.0400   0     0
                5   6   0.1     0.3    0.0300   0     0];   
%              bus     P        V
Datos.Gen = [   1      1        1.05     0
                2      0.5      1.05     0
                3      0.6      1.07     0];          
            
%               F   T    P        Q
Datos.Cargas = [1   3   0       0
                2   2   0       0
                3   2   0       0
                4   1   0.7     0.7
                5   1   0.7     0.7
                6   1   0.5     0.7];   %

Res = NR_Alg(Datos);
Vnt = abs(Res.V(:,:));
TheT = angle(Res.V(:,:));
PdT = Datos.Cargas(:,3);
QdT = Datos.Cargas(:,4);
VgT = Datos.Gen(2:end,3);
PgT = Datos.Gen(2:end,2);
N = size(Res.YBUS,1);
Ng = Datos.Gen(2:end,1);
Nd = find(Datos.Cargas(:,2) == 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pdc = zeros(6,1);
Pdc(Ng) = PgT;
Pdc = Pdc - PdT;
Pdc = Pdc(2:end);
DatosDC = Datos;
DatosDC.Lineas(:,3) = zeros(length(DatosDC.Lineas(:,1)),1);
DatosDC.Lineas(:,5) = zeros(length(DatosDC.Lineas(:,1)),1);
ResDC = NR_Alg(DatosDC);
Ybusdc = ResDC.YBUS(2:end,2:end);
Bbusdc = -imag(Ybusdc);
ThetaDCtem = Bbusdc\Pdc;
ThetaDC = ThetaDCtem; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
MostarResul = 1;
Ns = 1000;
MCS = 1;
ABC = 1;
JABC = 1;
ABCSMC = 1;
JABCSMC = 1;

sig2volprior = 0.0015;
DatosPrior.sig2Vol = sig2volprior; 
DatosPrior.aUni = -0.3;
DatosPrior.bUni = 0.3;

resamplingScheme = 1;       % The possible choices are
                            % systematic sampling (2),
                            % residual (1)
                            % and multinomial (3). 
                            % They're all O(N) algorithms.    

tsample = linspace(0.4,1.6,Ns);
tsample2 = linspace(min(TheT) - 0.5,max(TheT) + 0.5,Ns);
% Power generated
mhuPg = PgT;
sigmaPg = diag([0.025 0.03]);
Pg = mvnrnd(mhuPg',sigmaPg,Ns)';
mhuVg = VgT;
sigmaVg = 0.0001*eye(length(VgT));
Vg = mvnrnd(mhuVg',sigmaVg,Ns)';
% Loads
mhuPd = PdT(Nd);
sigmaPd = 0.035*eye(length(mhuPd));
Pd = mvnrnd(mhuPd',sigmaPd,Ns)';
mhuQd = QdT(Nd);
sigmaQd = 0.035*eye(length(mhuQd));
Qd = mvnrnd(mhuQd',sigmaQd,Ns)';

% These random variable were defined according to "Probabilistic Load-Flow 
% Computation using point estimate method" by Chun-Lien Su.
%% MCS
% We used MATPOWER to implement MCS.
% See http://www.pserc.cornell.edu/matpower/
if MCS ==1
x1 = zeros(Ns,N-1);
x2 = zeros(Ns,N-1);
DatosMCS = Datos;
Tppfi = cputime;
for n = 1:Ns
    DatosMCS.Gen(Ng,[1 2 3]) = [Ng  Pg(:,n)  Vg(:,n) ];
    DatosMCS.Cargas(Nd,[1 3 4]) = [Nd Pd(:,n) Qd(:,n)];
    Restemp = NR_Alg(DatosMCS);
    x1(n,:) = abs(Restemp.V(2:end,:))';
    x2(n,:) = angle(Restemp.V(2:end,:))';
    fprintf('\n Iteration #%d \n',n)
end
Tppff = cputime - Tppfi;

VolCPPF = x1;
ThetaCPPF = x2;

nvariv = 5;
nvarith = 5;
for k = 1:nvariv
   eval(['CPPFVol.v',int2str(k+1),'=ksdensity(x1(:,k),tsample);']); 
end

for k = 1:nvarith
   eval(['CPPFThe.th',int2str(k+1),'=ksdensity(x2(:,k),tsample2);']); 
end
end
%% ABCRej
if ABC == 1
    
epsilon = 0.7;
N = size(Res.YBUS,1);
x = zeros(Ns,2*N);
x(:,N+1:end) = ones(Ns,N);
x(:,2*N/2+1) = Vnt(1)*ones(Ns,1);
x(:,(2*N/2 + Ng)) = Vg';
ng = length(Ng);
NnodoPQtemp = Datos.Cargas(find(Datos.Cargas(:,2)== 1),1);
NnodoPQ = NnodoPQtemp;
nd = length(Nd);
linfT = ThetaDC - 0.07;
lsupT = ThetaDC + 0.07;
b = zeros((2*N - ng - 2),Ns);
b(1:ng,:) = b(1:ng,:) + Pg;
b(Nd-1,:) = b(Nd-1,:) - Pd;
b((N-ng-2+NnodoPQ),:) = b((N-ng-2+NnodoPQ),:) - Qd;
Rho = zeros(Ns,1);
xacep = zeros(1,2*N);
j = 1;
n = 1;
TppfABCi = cputime;
while j <= Ns
    x(n,2:N) =  unifrnd(linfT,lsupT,(N-1),1)';
    x(n,(N + Ng(end) +1):end) = ones(1,(N-length(Ng)-1)) + sqrt(sig2volprior)*randn(1,(N-length(Ng)-1));
    [PQ] = get_PotInj(x(n,:),Res.YBUS,ng);
    Rho(n) = sqrt(mean((b(:,n) - PQ).^2));
    if Rho(n)<=epsilon
        xacep(j,:) = x(n,:);
        fprintf('\n Number of parameters accepted #%d \n',j)
        j = j + 1;
    end
    if n == Ns
       n = 1;
    end
    fprintf('\n Iteration #%d \n',n)
    n = n + 1;
end
TppfABCf = cputime - TppfABCi;

ThetaABC = xacep(:,2:N);
VolABC = xacep(:,(N + Ng(end) +1):end);

nvariv = 3;
nvarith = 5;

for k = 1:nvariv
   eval(['ABCVol.v',int2str(k+3),'=ksdensity(VolABC(:,k),tsample);']); 
end

for k = 1:nvarith
   eval(['ABCThe.th',int2str(k+1),'=ksdensity(ThetaABC(:,k),tsample2);']); 
end
end

%% Jacobian ABC
if JABC == 1

eptrue = 0.7;

N = size(Res.YBUS,1);
x = zeros(Ns,2*N);
x(:,N+1:end) = ones(Ns,N);
x(:,2*N/2+1) = Vnt(1)*ones(Ns,1);
x(:,(2*N/2 + Ng)) = Vg';
ng = length(Ng);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xini = [ThetaDC' ones(1,N-length(Ng)-1)];
xbefore = zeros(Ns,2*N - ng - 2);
xbefore(1,:) = Xini;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NnodoPQtemp = Datos.Cargas(find(Datos.Cargas(:,2)== 1),1);
NnodoPQ = NnodoPQtemp;
nd = length(Nd);
b = zeros((2*N - ng - 2),Ns);
b(1:ng,:) = b(1:ng,:) + Pg;
b(Nd-1,:) = b(Nd-1,:) - Pd;
b((N-ng-2+NnodoPQ),:) = b((N-ng-2+NnodoPQ),:) - Qd;

sigma2vol = 0.000001;
sigma2theta = 0.0000001;
Sigma = zeros(2*N - ng - 2);
Sigma(1:N-1,1:N-1) = sigma2theta*eye(N-1);
Sigma(N:end,N:end) = sigma2vol*eye(N-length(Ng)-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pv = find(Datos.Cargas(:,2) == 2);
pq = find(Datos.Cargas(:,2) == 1);
npv = length(pv);
npq = length(pq);
j1 = 1;         j2 = npv;           %% j1:j2 - V angle of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      %% j3:j4 - V angle of pq buses
j5 = j4 + 1;    j6 = j4 + npq;      %% j5:j6 - V mag of pq buses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RhoC = zeros(Ns,1);
xacep = zeros(1,2*N);
j = 2;
n = 1;
m = 1;
TppfCABCi = cputime;
while j <= Ns
    xinitemp = mvnrnd(xbefore(j-1,:)',Sigma);
    x(n,2:N) = xinitemp(1:N-1);
    x(n,(N + Ng(end) +1):end) = xinitemp(1,N:end);
    [PQ] = get_PotInj(x(n,:),Res.YBUS,ng);
    RhoC(n) = sqrt((1/size(b(:,n),1))*(b(:,n) - PQ)'*(b(:,n) - PQ));%
    if RhoC(n)<=eptrue
        xacep(j,:) = x(n,:);
        V = x(n,N+1:end)'.*(exp(sqrt(-1).*x(n,1:N)'));
        [dx,~] = get_jacobian(V,PQ,b(:,n),Res.YBUS,pv,pq);
        xbefore(j,:) = xinitemp + dx';
        fprintf('\n Number of parameters accepted #%d \n',j)
        j = j + 1;
    end
    if n == Ns
       n = 1;
    end
    fprintf('\n Iteration #%d \n',m)
    n = n + 1;
    m = m + 1;
end
TppfCABCf = cputime - TppfCABCi;

ThetaJABC = xbefore(:,1:N-1);
VolJABC = xbefore(:,N:end);

nvariv = 3;
nvarith = 5;

for k = 1:nvariv
   eval(['JABCVol.v',int2str(k+3),'=ksdensity(VolJABC(:,k),tsample);']); 
end

for k = 1:nvarith
   eval(['JABCThe.th',int2str(k+1),'=ksdensity(ThetaJABC(:,k),tsample2);']); 
end

end
%% ABCSMC
if ABCSMC == 1
epT = [3.0 2.0 1.0 0.9 0.7];
NepT = length(epT);

N = size(Res.YBUS,1);
nd = length(Nd);
ng = length(Ng);
nbus = 2*N - ng - 2;

NnodoPQtemp = Datos.Cargas(find(Datos.Cargas(:,2)== 1),1);
NnodoPQ = NnodoPQtemp;
linfV = min(Vnt) - 0.01;
lsupV =  max(Vnt) + 0.01;
linfT = ThetaDC - 0.1;
lsupT = ThetaDC + 0.1;
b = zeros((2*N - ng - 2),Ns);
b(1:ng,:) = b(1:ng,:) + Pg;
b(Nd-1,:) = b(Nd-1,:) - Pd;
b((N-ng-2+NnodoPQ),:) = b((N-ng-2+NnodoPQ),:) - Qd;
xsmc1 = zeros(Ns,nbus);
xsmc2 = xsmc1;
xsmcacpe = zeros(1,nbus);
xsmcacpe2 = zeros(1,nbus);
Rho31 = zeros(Ns,1);
Rho32 = zeros(Ns,1);

xsmc = zeros(Ns,2*N);
xsmc(:,N+1) = Vnt(1)*ones(Ns,1);
xsmc(:,(N + Ng)) = Vg';
posiaceptasmc = zeros(1,1);
Wsmc = 1e-99*ones(NepT,Ns);
Ksmc = ones(NepT,Ns);
Vsmc = ones(NepT,Ns);

sigma2vol = 0.00001;
sigma2theta = 0.0001;
Sigma = zeros(nbus);
Sigma(1:N-1,1:N-1) = sigma2theta*eye(N-1);
Sigma(N:end,N:end) = sigma2vol*eye(N-length(Ng)-1);

U = unifrnd(0,1,1,Ns);
j = 1;
n = 1;
m = 1;
TppfABCSMCi = cputime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ABC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while j <= Ns
   xsmc1(n,1:N-1) = unifrnd(linfT,lsupT,N-1,1)';   
   xsmc1(n,N:end) = ones(1,(N-length(Ng)-1)) + sqrt(sig2volprior)*randn(1,(N-length(Ng)-1));
   
   xsmc(n,2:N) = xsmc1(n,1:N-1);
   xsmc(n,(N + Ng(end) +1):end) = xsmc1(n,N:end);
   
   [PQ] = get_PotInj(xsmc(n,:),Res.YBUS,ng);
   Rho31(n) = sqrt((1/size(b(:,n),1))*(b(:,n) - PQ)'*(b(:,n) - PQ));
   if Rho31(n) <= epT(1)
      xsmcacpe(j,:) = xsmc1(n,:);
      posiacepta(j) = n;
      Wsmc(1,j) = 1/Rho31(n);
      fprintf('\n Number of parameters accepted #%d \n',j)
      j = j + 1;
   end
   if n == Ns
      n = 1;
   end
   fprintf('\n Simulation #%d \n',m)
   n = n + 1;
   m = m + 1;
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ABC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wsmc(1,:) = Wsmc(1,:)./sum(Wsmc(1,:));
p = 1;
n = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SMC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
for j = 2:NepT
   if resamplingScheme == 1
     outIndex = residualR(1:length(Wsmc(j-1,:)),Wsmc(j-1,:)');        % Residual resampling.
   elseif resamplingScheme == 2
     outIndex = systematicR(1:length(Wsmc(j-1,:)),Wsmc(j-1,:)');      % Systematic resampling.
   else  
     outIndex = multinomialR(1:length(Wsmc(j-1,:)),Wsmc(j-1,:)');     % Multinomial resampling.  
   end;
   xsmcpre = xsmcacpe(outIndex,:);
while p <= Ns %
    xsmc2(n,:) = mvnrnd(xsmcpre(n,:)',Sigma);
    Ksmc(j,n) = mvnpdf(xsmc2(n,:),xsmcpre(n,:),Sigma);
    
    xsmc(n,2:N) = xsmc2(n,1:N-1);
    xsmc(n,(N + Ng(end) +1):end) = xsmc2(n,N:end);
    
    [PQ] = get_PotInj(xsmc(n,:),Res.YBUS,ng);
    sigPQ =  mad(PQ);
    Rho32(n) = sqrt((1/size(b(:,n),1))*(b(:,n) - PQ)'*(b(:,n) - PQ));
   if Rho32(n) <= epT(j)
      xsmcacpe2(p,:) = xsmc2(n,:);
      posiaceptasmc(p) = n;
      pitheta = (prod(unifpdf(xsmc2(1,1:N-1)',linfT,lsupT)))...
                *mvnpdf(ones(1,N-length(Ng)-1),xsmc2(n,N:end),sig2volprior*eye(3));
      denWsmc = sum(Wsmc(j-1,:).*Ksmc(j-1,:));
      Wsmc(j,p) = pitheta/denWsmc + 1e-99;
      fprintf('\n Number of parameters accepted #%d \n',p)
      p = p + 1;  
   end
   if n == Ns
      n = 1; 
   end
   fprintf('\n Simulation #%d \n',m)
   n = n + 1;
   m = m + 1;
end
Wsmc(j,:) = Wsmc(j,:)./sum(Wsmc(j,:));
p = 1;
xsmcacpe = xsmcacpe2;
n = 1;
end
TppfABCSMCf = cputime - TppfABCSMCi;
xacepmcmc = xsmcacpe;

ThetaABCSMC = xacepmcmc(:,1:N-1);
VolABCSMC = xacepmcmc(:,N:end);

for k = 1:nvariv
   eval(['ABCSMCVol.v',int2str(k+3),'=ksdensity(VolABCSMC(:,k),tsample);']); 
end

for k = 1:nvarith
   eval(['ABCSMCThe.th',int2str(k+1),'=ksdensity(ThetaABCSMC(:,k),tsample2);']); 
end

end
%% JABCSMC
if JABCSMC == 1
epT = [3.0 2.0 1.0 0.9 0.7];

NepT = length(epT);

N = size(Res.YBUS,1);
nd = length(Nd);
ng = length(Ng);
nbus = 2*N - ng - 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xini = [ThetaDC' ones(1,N-length(Ng)-1)];
xbefore = zeros(Ns,2*N - ng - 2);
xbefore(1,:) = Xini;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NnodoPQtemp = Datos.Cargas(find(Datos.Cargas(:,2)== 1),1);
NnodoPQ = NnodoPQtemp;
linfV = min(Vnt) - 0.01;
lsupV =  max(Vnt) + 0.01;
linfT = ThetaDC - 0.03;
lsupT = ThetaDC + 0.03;
b = zeros((2*N - ng - 2),Ns);
b(1:ng,:) = b(1:ng,:) + Pg;
b(Nd-1,:) = b(Nd-1,:) - Pd;
b((N-ng-2+NnodoPQ),:) = b((N-ng-2+NnodoPQ),:) - Qd;
xsmc1 = zeros(Ns,nbus);
xsmc2 = xsmc1;
xsmcacpe = zeros(1,nbus);
xsmcacpe2 = zeros(1,nbus);
Rho41 = zeros(Ns,1);
Rho42 = zeros(Ns,1);

xsmc = zeros(Ns,2*N);
xsmc(:,N+1) = Vnt(1)*ones(Ns,1);
xsmc(:,(N + Ng)) = Vg';
posiaceptasmc = zeros(1,1);
Wsmc = 1e-99*ones(NepT,Ns);
Ksmc = ones(NepT,Ns);
Vsmc = ones(NepT,Ns);

sigma2vol = 0.000001;
sigma2theta = 0.0000001;
Sigma = zeros(nbus);
Sigma(1:N-1,1:N-1) = sigma2theta*eye(N-1);
Sigma(N:end,N:end) = sigma2vol*eye(N-length(Ng)-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pv = find(Datos.Cargas(:,2) == 2);
pq = find(Datos.Cargas(:,2) == 1);
npv = length(pv);
npq = length(pq);
j1 = 1;         j2 = npv;           %% j1:j2 - V angle of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      %% j3:j4 - V angle of pq buses
j5 = j4 + 1;    j6 = j4 + npq;      %% j5:j6 - V mag of pq buses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = unifrnd(0,1,1,Ns);
j = 2;
n = 1;
m = 1;
TppfJABCSMCi = cputime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%JABC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while j <= Ns
   xinitemp = mvnrnd(xbefore(j-1,:)',Sigma);
   xsmc1(n,1:N-1) = xinitemp(1:N-1);    
   xsmc1(n,N:end) = xinitemp(1,N:end); 
   
   xsmc(n,2:N) = xsmc1(n,1:N-1);
   xsmc(n,(N + Ng(end) +1):end) = xsmc1(n,N:end);
   
   [PQ] = get_PotInj(xsmc(n,:),Res.YBUS,ng);
   Rho41(n) = sqrt((1/size(b(:,n),1))*(b(:,n) - PQ)'*(b(:,n) - PQ));
   if Rho41(n) <= epT(1)
      xsmcacpe(j,:) = xsmc1(n,:);
      V = xsmc(n,N+1:end)'.*(exp(sqrt(-1).*xsmc(n,1:N)'));
      Va = angle(V);
      Vm = abs(V);
      [dx,~] = get_jacobian(V,PQ,b(:,n),Res.YBUS,pv,pq);
      xbefore(j,:) = xinitemp + dx';
      Wsmc(1,j) = 1/Rho41(n);
      fprintf('\n Number of parameters accepted #%d \n',j)
      j = j + 1;
   end
   if n == Ns
      n = 1;
   end
   fprintf('\n Iteration #%d \n',m)
   n = n + 1;
   m = m + 1;
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%JABC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wsmc(1,:) = Wsmc(1,:)./sum(Wsmc(1,:));
p = 1;
n = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SMC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
for j = 2:NepT
   if resamplingScheme == 1
     outIndex = residualR(1:length(Wsmc(j-1,:)),Wsmc(j-1,:)');        % Residual resampling.
   elseif resamplingScheme == 2
     outIndex = systematicR(1:length(Wsmc(j-1,:)),Wsmc(j-1,:)');      % Systematic resampling.
   else  
     outIndex = multinomialR(1:length(Wsmc(j-1,:)),Wsmc(j-1,:)');     % Multinomial resampling.  
   end;
   xsmcpre = xbefore(outIndex,:);
while p <= Ns
    xsmc2(n,:) = mvnrnd(xsmcpre(n,:)',Sigma);
    Ksmc(j,n) = mvnpdf(xsmc2(n,:),xsmcpre(n,:),Sigma);
    
    xsmc(n,2:N) = xsmc2(n,1:N-1);
    xsmc(n,(N + Ng(end) +1):end) = xsmc2(n,N:end);
    
    [PQ] = get_PotInj(xsmc(n,:),Res.YBUS,ng);
    Rho42(n) = sqrt((1/size(b(:,n),1))*(b(:,n) - PQ)'*(b(:,n) - PQ));
   if Rho42(n) <= epT(j)
      xsmcacpe2(p,:) = xsmc2(n,:);
      posiaceptasmc(p) = n;
      pitheta = (mvnpdf(ThetaDC',xsmc2(n,1:N-1),sigma2theta*eye(5)))...
                *mvnpdf(ones(1,N-length(Ng)-1),xsmc2(n,N:end),sigma2vol*eye(3));
      denWsmc = sum(Wsmc(j-1,:).*Ksmc(j-1,:));
      Wsmc(j,p) = pitheta/denWsmc + 1e-99;
      fprintf('\n Number of parameters accepted #%d \n',p)
      p = p + 1;  
   end
    if n == Ns
      n = 1;
   end
   fprintf('\n Iteration #%d \n',m)
   n = n + 1;
   m = m + 1;
end
Wsmc(j,:) = Wsmc(j,:)./sum(Wsmc(j,:));
p = 1;
xbefore = xsmcacpe2;
n = 1;
end
TppfJABCSMCf = cputime - TppfJABCSMCi;
xacepmcmc = xbefore;

ThetaJABCSMC = xacepmcmc(:,1:N-1);
VolJABCSMC = xacepmcmc(:,N:end);

for k = 1:nvariv
   eval(['JABCSMCVol.v',int2str(k+3),'=ksdensity(VolJABCSMC(:,k),tsample);']); 
end

for k = 1:nvarith
   eval(['JABCSMCThe.th',int2str(k+1),'=ksdensity(ThetaJABCSMC(:,k),tsample2);']); 
end

end
%% Results

figure
hax=axes;
hold on
if MCS == 1
plot(tsample,CPPFVol.v6,'--r','linewidth',2)
end
if ABC == 1
plot(tsample,ABCVol.v6,'--b','linewidth',2)
end
if JABC == 1
plot(tsample,JABCVol.v6,'b','linewidth',2)
end
if ABCSMC == 1
plot(tsample,ABCSMCVol.v6,'--m','linewidth',2)
end
if JABCSMC == 1
plot(tsample,JABCSMCVol.v6,'m','linewidth',2)
end
line([Vnt(6,1) Vnt(6,1)],get(hax,'YLim'),'color',[0 0 1],'linewidth',2)
xlim([0.9 1.1])
legend('MCS','ABC','JABC','ABCSMC','JABCSMC','DetSol')
xlabel('Voltages [V]')
ylabel('Density')

figure
hax=axes;
hold on
if MCS == 1
plot(tsample2,CPPFThe.th6,'--r','linewidth',2)
end
if ABC == 1
plot(tsample2,ABCThe.th6,'--b','linewidth',2)
end
if JABC == 1
plot(tsample2,JABCThe.th6,'b','linewidth',2)
end
if ABCSMC == 1
plot(tsample2,ABCSMCThe.th6,'--m','linewidth',2)
end
if JABCSMC == 1
plot(tsample2,JABCSMCThe.th6,'m','linewidth',2)
end
line([TheT(6,1) TheT(6,1)],get(hax,'YLim'),'color',[0 0 1],'linewidth',2)
xlim([-0.3 0.2])
legend('MCS','ABC','JABC','ABCSMC','JABCSMC','DetSol')
xlabel('Angles [rad]')
ylabel('Density')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('--------------------------------------------------------------\n')
fprintf('Comp. Time\n')
if MCS == 1
fprintf('CTime for MCS = %f\n',Tppff)
end
if ABC == 1
fprintf('CTime for ABC = %f\n',TppfABCf)
end
if JABC == 1
fprintf('CTime for JABC = %f\n',TppfCABCf)
end
if ABCSMC == 1
fprintf('CTime for ABCSMC = %f\n',TppfABCSMCf)
end
if JABCSMC == 1
fprintf('CTime for JABCSMC: %f\n',TppfJABCSMCf)
end
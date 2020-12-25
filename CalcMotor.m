clear
clc
close

%% Declaração das constantes gerais

t1 = 0:0.01:150; %tempo da simulação após o termino da queima de propelente
Pe = 101325; %pressão externa (ambiente) (Pa)
g = 9.81; %gravidade (m/s^2)
pho_ar = 1.2754; %massa especifica do ar (kg/m^3)
RF = 0.060/2; %raio do foguete (m)

%% Declaração das constantes do motor/propelente

h = 0.01; %altura da garganta (m)
en = 0.004; %espessura da parede do nozzel(m)
em = 0.001325; %espessura da parede do motor(m)
ei = 0.000; %espessura do isolamento(m)
ee = 0.006; %espessura dos separadores somados(m)
alpha = pi/3; %angulo seção convergente (rad)
beta = pi/6; %angulo seção convergente (rad)
pho_n = 7860; %dencidade do metal do nozzel (kg/m^3)
%pho_m = 2698.79071; %dencidade do metal(6061-T6) do motor (kg/m^3)
pho_m = 7861.09294; %dencidade do metal(A36) do motor (kg/m^3)

n = 8; %numero de grãos
alt = 0.07; %altura do grão (m)
R0 = 0.015/2; %raio interno do grão (m)
RE = 0.059/2; %raio externo do grão (m)
pho = 1841; %densidade do propelente (kg/m^3)
C1 = 930; %velocidade cracteristica do propelente (C*) (m/s)
r = 0.0113; %taxa de queima (m/s)
P0i = 6890100; %pressão da camera esperada (Pa)
k = 1.042; %constante k do propelente

%% Declaração das constantes do foguete

cd = 0.5; %coeficiente de arrasto
M0 = 5; %massa seca (kg)

%% Calculo das informações basicas do foguete

L = n*alt; %altura total do propelente
S = pi*(RF)^2; %area molhada (para calculo de arrasto) (m^2)
K = cd*S*pho_ar/2; %constante de arrasto (kg/m)
Mp = pi*(RE^2-R0^2)*L*pho; %massa de propelente (kg)

%% Funções do foguete

Aq = @(t) (2*pi*(R0 + r*t)*(L-2*n*r*t)+pi*(RE^2-(R0+r*t)^2)); %área de queima com as pontas dos grão não inibidas
%Aq = @(t) (2*pi*(R0 + r*t)*(L)); %área de queima com as pontas dos grão inibidas
Ae = r * pho * C1 * Aq(0) / P0i; %área da garganta
R = sqrt(Ae/pi); %raio da garganta
As = Ae/(((k+1)/2)^(1/(k-1))*(Pe/P0i)^(1/k)*sqrt((k+1)/(k-1)*(1-(Pe/P0i)^((k-1)/k)))); %área de saida do nozzel
RS = sqrt(As/pi); %raio da saida do nozzel

%%Calculo da massa do nozzel
Vncc = pi*tan(alpha)*(RE-R)*(RE^2 + RE*R + R^2); %volume da seção convergente "cheia" do nozzel
Vnc = Vncc - pi*tan(alpha)*((RE-en)-(R-en))*((RE-en)^2 + (RE-en)*(R-en) + (R-en)^2); %volume da seção convergente
Vng = pi*(RE^2 - R^2)*h; %volume da seção da garganta do nozzel
Vndc = pi*tan(beta)*(R-RS)*(R^2 + R*RS + RS^2); %volume da seção divergente "cheia" do nozzel
Vnd = Vndc - pi*tan(beta)*((R-en)-(RS-en))*((R-en)^2 + (R-en)*(RS-en) + (RS-en)^2); %volume da seção divergente
Vn = Vnc + Vng + Vnd; %volume total do nozzel
Mn = Vn*pho_n; %massa total do nozzel
%%Fim do calculo da massa do nozzel

%%Calculo da massa do motor
Mm = pi*((RE+ei+em)^2-(RE+ei)^2)*(L+(n-1)*ee)*pho_m;
%%Fim do calculo da massa motor

P0 = @(t) (r*pho*C1*Aq(t)/Ae); %pressão da camera no decorrer do tempo
Cf = @(t) (sqrt((2*k^2)/(k-1)*(2/(k+1))^((k+1)/(k-1))*(1-(Pe/P0(t))^((k-1)/k)))); %coeficiente de força no decorrer do tempo
m = @(t) (-r*pho*Aq(t)); %fluxo de massa no decorrer do tempo
M = @(t) (M0+Mn+Mm+Mp+m(t)*t); %massa total no decorrer do tempo


%% Forças atuantes no foguete

E = @(t) (Cf(t)*Ae*P0(t)); %empuxo
P = @(t) (-M(t)*g); %peso
D = @(v) (-K*v^2); %arrasto

F = @(t,v) (E(t)+P(t)+D(v)); %força resultante

%% Equações de movimento

%baseado na equação da força F = a*M + v*m, substituindo a força pela
%força resultante e isolando a aceleração obtem-se:

a = @(t,v) ((F(t,v)-v*m(t))/M(t)); %aceleração dependente da velocidade e do tempo

%%Resolução da edo por ode45

if(((RE-R0)/r) < (L/(2*n*r)))
    tmax = (RE-R0)/r;
else
    tmax = L/(2*n*r);
end
yp=@(t,y)[y(2);a(t,y(2))]; 
y0=[0;0];
tspan=[0 tmax];
options=odeset('RelTol',1e-12);
[t,y]=ode45(yp,tspan,y0,options);

%%Fim da resolução da edo, y(:,1) é o deslocamento e y(:,2) é a velocidade

M2 = (M0+Mn+Mm+Mp); %massa total depois da queima
P2 = (-M2*g); %peso
F2 = @(v) (P2+D(v)); %força depois da queima
a2 = @(v) ((F2(v))/M2); %aceleração depois da queima

%%Resolução da edo por ode45
yp=@(t,y)[y(2);a2(y(2))]; 
y0=[y(length(y),1);y(length(y),2)];
tspan=[tmax 40];
options=odeset('RelTol',1e-12);
[t2,y2]=ode45(yp,tspan,y0,options);
%%Fim da resolução da edo, y2(:,1) é o deslocamento e y2(:,2) é a velocidade

%X = @(t) (y(length(y),1) + y(length(y),2)*t - g*t.^2/2); %equação horaria do movimento após a queima
%V = @(t) (y(length(y),2) - g*t); %equação horaria da velocidade após a queima

%% Plots

plot(t,y(:,1),'r') %plot do deslocamento no tempo durante a queima
hold on
plot(t2,y2(:,1),'b') %plot do deslocamento no tempo após a queima
title('Deslocamento x Tempo')
xlabel('t')
ylabel('x(t)')
 
figure

plot(t,y(:,2),'r') %plot da velocidade no tempo durante a queima
hold on
plot(t2,y2(:,2),'b') %plot da velocidade no tempo após a queima
title('Velocidade x Tempo')
xlabel('t')
ylabel('v(t)')

figure

Em = [];

for i = 1:length(t)
    Em = [Em;E(t(i))];
end

plot(t,Em) %plot do empuxo no tempo durante a queima
title('Empuxo x Tempo')
xlabel('t')
ylabel('E(t)')

figure

P0m = [];

for i = 1:length(t)
    P0m = [P0m;P0(t(i))];
end

plot(t,P0m) %plot da pressão da camera no tempo durante a queima
title('Pressão na camera x Tempo')
xlabel('t')
ylabel('P0(t)')

%% Resultados
disp('Informações de voo')
disp(['Apogeu máximo(m): ',num2str(max(y2(:,1)))])
disp(['Velocidade máxima(m/s): ',num2str(max(y(:,2)))])
disp(['Velocidade máxima(mach): ',num2str(max(y(:,2))/340)])
disp(' ')
disp('Informações do Nozzel')
disp(['Área da garganta(m^2): ',num2str(Ae)])
disp(['Diametro da garganta(m): ',num2str(2*R)])
disp(['Área da saida(m^2): ',num2str(As)])
disp(['Diametro da saida(m): ',num2str(2*RS)])
disp(' ')
disp('Informações do Motor')
disp(['Tempo de queima(s): ',num2str(max(t))])
disp(['Empuxo máximo(N): ',num2str(max(Em))])
disp(['Empuxo mínimo(N): ',num2str(min(Em))])
disp(['Pressão máxima da camera(Pa): ',num2str(max(P0m))])
disp(['Pressão mínimo da camera(Pa): ',num2str(min(P0m))])
disp(['Massa de seca (sem estrutura)(kg): ',num2str(Mn+Mm)])
disp(['Altura do motor(m): ',num2str(L+(n-1)*ee)])

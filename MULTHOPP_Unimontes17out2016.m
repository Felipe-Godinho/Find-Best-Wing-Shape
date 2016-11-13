% Methodo Multhopp Generico para N aeronaves.
% Utilizado para encontrar a distribuicao de sustentacao ao longo de um asa
% qualquer, inclinação da curva CL x alfa, alem do calculo do arrasto
% induzido respectivo.
 
 
%% DADOS DE ENTRADA
M = 1.0901515;                 % Massa especifica do ar a 1200m de altitude( Kg/m^3)
g = 9.8067;                    % aceleracao da gravidade (m/s^2)
Mi = 1.7516*10^-5;             % Viscosidade dinamica do ar a 1200 metros de altitude (N.s/m^2) 
alfa = 12.7;                    % angulo geometrico da asa em graus
alfaL0= 0;                     % angulo de sustentação nula do perfil Naca0012
a0 = 5.6;                      % inclinação da curva cl x alfa do perfil Naca0012 (rad^-1)
m=1001;                         % Número de estações (sempre ímpar)

%% DADOS AERONAVE 
b=20;                            % envergadura da asa (m)
S=40;                           % area da asa (m^2)
cordas=[3 1];                   % vetor de cordas de raiz das semi-asas (m)
afilamento=cordas(2)/cordas(1);

%% DECLARACAO DAS VARIAVEIS
alphai=zeros(1,m);                              % criando a matriz de angulos induzidos
%alpha = ones(1,m)*(alfaL0 - alfa)/ 57.3;        % Angulo de ataque efetivo para cada estacao
theta = ([m:-1:1].*pi/(m+1));                   % Angulo encontrado devido a posicao das 'm' cordas para uma distribuicao eliptica
eta=cos(theta);                                 % Posicao das 'm' cordas em relacao ao centro da asa
y = 0.5 * b * eta;                              % Posicao absoluta das 'm' Cordas
delta=(y(1,2:m)-y(1,1:m-1));                    % intervalos entre as 'm' cordas para o Multhopp - variavel 'y'
deta=(eta(2:m)-eta(1:m-1));                     % intervalos entre as 'm' cordas - variavel 'eta'

%% ENCONTRANDO AS CORDAS EM CADA ESTACAO "m" - Exemplo Pullin - página: II-8
cv=3-(2*abs(eta));                              % denominando as cordas da asa para as "m" estações
% cv=cordas(1)*ones(1,m);

% Outros Parâmetros Multhopp
lambda=b./cv;                                   % lambda na estação "nu"
bv=2.*lambda./a0;                               % envergadura na estação "nu"
Sasa=((cv(:,1:m-1)+cv(:,2:m))/2).*delta;        % area do somatório das "nu" secoes (m^2)
%% INICIO DO PROGRAMA DE CIRCULAÇÃO - MULTHOPP

clAlfa = zeros(31,1);
ciAlfa = zeros(31,1);
j = 1;
for i=-15:15
    alfa = ones(1,m)*deg2rad(alfaL0 - 10);
    for t=1:m
        for k = 1:m
            if k == t
                bvn(k,t) = -((m+1)/(4*sin(theta(k))));
                B(k,t)=-((m+1)/(4*sin(theta(k)))+bv(1,k));
            else
                B1(k,t)=( sin(theta(t)) / (cos(theta(t)) - cos(theta(k))).^2 );
                B2(k,t) = ( (1-(-1)^(t-k)) / (2*(m+1)) );
                B(k,t) = B1(k,t) * B2(k,t);
                bvn(k,t) = B1(k,t) * B2(k,t);
            end
        end
    end   
    % Calculo da Circulacao
    
    circ = inv(B) * (alfa');

    % Calculo de Coeficiente Sustentacao
    
    CL = zeros(m,1);
    CLmed = zeros(m,1);
    for t=1:m

        CL(t) = 2 * circ(t) * lambda(1,t);

        if t>1

            CLmed(t-1)= (CL(t-1) + CL(t))/2;

        end

    end

    %% RESULTADOS DA ASA GERADA
        
        alongamento = (b^2)./S;                                                             % alongamento da asa
        alphai(1,:) = -(bvn*circ)';                                                         % angulo induzido em cada estacao
        clAlfa(j) = (sum(CLmed(1,:).*Sasa(1,:)))./(sum(Sasa));                                % coeficiente de sustentacao da asa
        gama=circ';                                                                         % efeito de circulacao em cada estacao, provocando os vortices
        Fx_obj = gama.*alphai(1,:);                                                         % determinando a funcao a ser integrada
        Cdilocal=alongamento*((Fx_obj(2:(m+1)/2)+Fx_obj(1:(m+1)/2-1)).*deta(1:(m-1)/2)/2);  % induzido em cada estacao
        ciAlfa(j)=2*sum(Cdilocal(1,:));                                                      % coeficiente de arrasto induzido da asa
        E = clAlfa(j)/ciAlfa(j)                                                               % eficiencia da asa devido o induzido.
        j = j+1;
end
        
%%PLOTANDO RESULTADOS
figure(2);subplot(4,1,1); plot(y(1,:),cv(1,:)/4,y(1,:),-cv(1,:)*3/4);
figure(2);subplot(4,1,2); plot((-15:15),clAlfa);
figure(2);subplot(4,1,3); plot((-15:15),ciAlfa);
figure(2);subplot(4,1,4); plot(y(1,:),CL,'o-r');

        
   
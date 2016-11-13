function[E,L] = fitAsa(envergaduras,cordas)
%% DADOS DE ENTRADA

    load sefx;
    perfil = sefx;
    [~,pos] = min(abs(perfil(:,2)));  % angulo de sustentação nula do perfil
    alfaL0 = perfil(pos,1);
    [~,a] = max(perfil(:,2)./perfil(:,3));
    alfa = perfil(a,1);              %angulo de maxima eficiencia para decolagem
    m = 51;                          % Número de estações (sempre ímpar)
    alpha = ones(1,m)*deg2rad(alfaL0 - alfa);
    M = 1.0901515;                    % Massa especifica do ar a 1200m de altitude( Kg/m^3)
    g = 9.8067;                       % aceleracao da gravidade (m/s^2)
    Mi = 1.7516*10^-5;                % Viscosidade dinamica do ar a 1200 metros de altitude (N.s/m^2)
    v = 15 ;                          % velocidade maxima em m/s
    aux = find(sefx(:,1)==9);       
    a0 = rad2deg((sefx(aux,2)-sefx(pos,2))/(sefx(aux,1)-sefx(pos,1)));                           

    %% DADOS AERONAVE  
    n=size(envergaduras,2);            %número de seções 
    b = 1;                             %definir a envergadura da asa m
    b1 = envergaduras(:,1:n);  
    envergaduras = 2*b.*envergaduras(:,1:n);
    S = envergaduras(1,1)*cordas(1);       % area da região retangular
    for i=1:n-1
        S = S + (cordas(i)+cordas(i+1))*envergaduras(1,i+1)*0.5; %area da região trapesoidal
    end
    %afilamento = 0.618;
    %% DECLARACAO DAS VARIAVEIS
    alphai=zeros(1,m);                              % criando a matriz de angulos induzidos
    theta = ((m:-1:1).*pi/(m+1));                   % Angulo encontrado devido a posicao das 'm' cordas para uma distribuicao eliptica
    eta=cos(theta);                                 % Posicao das 'm' cordas em relacao ao centro da asa
    y = 0.5 * b * eta;                              % Posicao absoluta das 'm' Cordas
    delta=(y(1,2:m)-y(1,1:m-1));                    % intervalos entre as 'm' cordas para o Multhopp - variavel 'y'
    deta=(eta(2:m)-eta(1:m-1));                     % intervalos entre as 'm' cordas - variavel 'eta'

    %% ENCONTRANDO AS CORDAS EM CADA ESTACAO "m" 
    cv=ones(1,m);
    for i=1:m
        if abs(eta(i))<=2*b1(1,1)
            cv(1,i)=cordas(1);                              % denominando as cordas da asa para as "m" estações    
        elseif 2*b1(1,1)<abs(eta(i))&&abs(eta(i))<=2*sum(b1(1,1:2))
            cv(1,i)=.5*abs(eta(i))*(cordas(2)-cordas(1))/(sum(b1(1,1:2))-b1(1,1)) + ((-b1(1,1)*(cordas(2)-cordas(1))/(sum(b1(1,1:2))- b1(1,1))) + cordas(1));
        elseif 2*sum(b1(1,1:2))<abs(eta(i))
            cv(1,i)=.5*abs(eta(i))*(cordas(3)-cordas(2))/(sum(b1(1,1:3))- sum(b1(1,1:2))) + ((-sum(b1(1,1:2))*(cordas(3)-cordas(2))/(sum(b1(1,1:3))- sum(b1(1,1:2)))) + cordas(2));
        end    
    end    
    % Outros Parâmetros Multhopp
    lambda=b./cv;                                   % lambda na estação "nu"
    bv=2.*lambda./a0;                               % envergadura na estação "nu"
    Sasa=((cv(:,1:m-1)+cv(:,2:m))/2).*delta;        % area do somatório das "nu" secoes (m^2)

    %% INICIO DO PROGRAMA DE CIRCULAÇÃO - MULTHOPP

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
        
        circ = inv(B) * (alpha');
        
        % Calculo de Coeficiente Sustentacao
        
        for t=1:m
            
            CL(1,t) = 2 * circ(t) * lambda(1,t);
            
            if t>1
                
                CLmed(1,t-1)= (CL(1,t-1) + CL(1,t))/2;
                
            end
            
        end
        
%% RESULTADOS DA ASA GERADA
        
        alongamento = (b^2)./S;                                                             % alongamento da asa
        alphai(1,:) = -(bvn*circ)';                                                         % angulo induzido em cada estacao
        CLtotal = (sum(CLmed(1,:).*Sasa(1,:)))./(sum(Sasa));                                % coeficiente de sustentacao da asa
        gama=circ';                                                                         % efeito de circulacao em cada estacao, provocando os vortices
        Fx_obj = gama.*alphai(1,:);                                                         % determinando a funcao a ser integrada
        Cdilocal=alongamento*((Fx_obj(2:(m+1)/2)+Fx_obj(1:(m+1)/2-1)).*deta(1:(m-1)/2)/2);  % induzido em cada estacao
        Cditotal=2*sum(Cdilocal(1,:));                                                      % coeficiente de arrasto induzido da asa
        L = CLtotal*S*M*.5*v^2;
        %CD = CDp*Se/Sw
        Cdpasa = perfil(a,4)*sum(Sasa)/b*0.1*0.3;
        E = CLtotal/(Cditotal+Cdpasa);
        if L<30
            E = -Inf;
        end    
end       
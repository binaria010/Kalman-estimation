final_k = 1000;
n = 3;
rank(obsv(dyn_op(1),obs_op(1)))
eig(dyn_op(1))




x_est = zeros(n,1) %Vector columna
%n = size(x_est);
P = 10*eye(n); %Covarianza de la estimacion (generalmente es de la forma M*eye(n))
R0 = cov_obs_noise(0); %Covarianza del ruido en la observacion
m = 1; %Dimension del espacio de mediciones
y = zeros(final_k+1,m); 
utilde =  0; 
utt = zeros(n,1);
Rtilde = R0;
u = zeros(final_k+1,1)%rand(final_k+1,1); %la coordenada k corresponde al tiempo k-1

%S(0)

pi_op = eye(n);

x = rand(n,1);% Inicializo el estado real del sistema
x0 = x;
y(1)= obs_op(0)* x + obs_noise(0); %Observacion en tiempo 0 
er_vect = zeros(final_k+1,1);
for k = 0:final_k
    %Observacion
    Htilde = obs_op(k)*pi_op; %Htilde en tiempo k
    Rtilde = cov_obs_noise(k); %Rtilde en tiempo k, esto no esta implementado ruido en la dinamica.
    K = (P*Htilde')/(Htilde*P*Htilde'+Rtilde); %Faltan los terminos correspondientes al S
    ytilde = y(k+1)-utilde; %y(k+1) es la observacion en tiempo k
    utt = dyn_op(k)*utt+cont_op(k)*u(k+1); %u(k+1) es el control en tiempo k, utt es en tiempo k
    utilde = obs_op(k+1)*utt; %utilde en tiempo k 
    x_est = x_est + K * (ytilde - Htilde*x_est); %estimacion k-esima
    P = (eye(n)-K*Htilde)*P; %Faltan los terminos del Sk cuando hay ruido en la dyn
    pi_op = dyn_op(k)*pi_op; %calculo pi_op en tiempo k
    
    er_vect(k+1) = norm(x_est-x0);
    
    %Evolucion del sistema
    x = dyn_op(k)*x+cont_op(k)*u(k+1); %Estado en tiempo k+1
    y(k+2) = obs_op(k+1)*x + obs_noise(k+1);%Observacion en tiempo k+1, se podria hacer recursivo
end

x
x0
x_est
norm(x_est-x0)
P
plot(er_vect)

function A_k = dyn_op(k)
    A_k = (1/16)*[1 2 3; 2 5 6; 7 8 9]; 
end

function h_k = obs_op(k)
    n = 3; %Dar la dimension del espacio de estados
    l = floor(n/2);
    h_k = zeros(1,n);
    h_k(l) = 1;
end

function B_k = cont_op(k)
    %Debe ser una matriz de tamaño n * el tamaño del control
    n = 3; %Dar la dimension del espacio de estados
    B_k = ones(n,1);
end

function R_k = cov_obs_noise(k)
    m = 1; %Dar la dimension del espacio de mediciones
    sigma = 0.02;
    R_k = sigma*eye(m);
end

function v_k = obs_noise(k)
    %Devuelve el ruido para la medicion en el instante k
    m = 1; %Dar la dimension del espacio de estados
    v_k = mvnrnd(zeros(m,1),cov_obs_noise(k),1); 
end

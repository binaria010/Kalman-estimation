%% Fixed Point Smoother for estimation at fixed time j=0

close all
clear all

% initialize filter:
j=1;
num_obs = 100;
m=1; n=2;
T = 0.1;
a= 0.2;
Q =a^2 * [(T^4)/4 (T^3)/2;(T^3)/2 T^2]; % noise in the dynamics
sigma = 1; P = sigma*eye(n); %Covarianza inicial del error 
Sigma_j = P; Pi_j = P; 
rng(1)
x = 100*rand(n,1);% Inicializo el estado real del sistema
x_0 = x;
xhat = zeros(n,1); % este x es la inicializacion de arriba, no se si esta bien hacer eso
xhat_j = xhat;
y = zeros(num_obs +1 ,m); 
y(1)= obs_op(0)* x + obs_noise(0); %Observacion en tiempo 0

% figure
% subplot(2,1,1),
% %plot(x_0,'r'), title("Estimacion estado inicial")
% hold on
% pause(0.05)

error1=zeros(num_obs+1,1);
error2=zeros(num_obs+1,1);

for k = j-1:num_obs
    % gain matrices:
    L = dyn_op(k)*P*(obs_op(k))'/(obs_op(k)*P*obs_op(k)' + cov_obs_noise(k));
    lambda = Sigma_j*obs_op(k)'/(obs_op(k)*P*obs_op(k)' + cov_obs_noise(k));
    
    pi_trace(k+1) = trace(Pi_j);
    p_trace(k+1) = trace(P);
    
    %update estimations
    xhat_j = xhat_j +lambda*(y(k+1)-obs_op(k)*xhat);
    xhat = dyn_op(k)*xhat + L*(y(k+1)-obs_op(k)*xhat);
%     subplot(2,1,1),
%     plot(xhat_j,'g'), title("Estimacion estado inicial")
%     subplot(2,1,2),
%     plot(xhat-x,'g'), title("Error estimacion estado actual")
    %pause(0.05)
    
    error1(k+1)=abs(xhat_j(1)-x_0(1));
    error2(k+1)=abs(xhat_j(2)-x_0(2));
    
    % update covariances
    P = dyn_op(k)*P*(dyn_op(k) - L*obs_op(k))' + Q;
    Pi_j = Pi_j - Sigma_j*obs_op(k)'*lambda';
    Sigma_j = Sigma_j*(dyn_op(k) - L*obs_op(k))';
    
    % propagate
    x = dyn_op(k)*x + dyn_noise(k)';
    y(k+2) = obs_op(k+1)*x + obs_noise(k+1); %Observacion en tiempo k+1
end

figure(1)
subplot(2,1,1),
plot(pi_trace), title("trace of the smoothed covariance")
subplot(2,1,2),
plot(p_trace), title("trace of the covariance")

figure(2)
plot(0:0.1:10,error1,'g',0:0.1:10,error2,'b'), title("error en la primera y segunda coordenadas")
%%

function A_k = dyn_op(~)
    T = 0.1;
    A_k = [1 T; 0 1]; 
end

function h_k = obs_op(~)
    n = 2; %Dar la dimension del espacio de estados
    l = floor(n/2);
    h_k = zeros(1,n);
    h_k(l) = 1;
end

function B_k = cont_op(~)
    %Debe ser una matriz de tamaño n * el tamaño del control
    n = 2; %Dar la dimension del espacio de estados
    B_k = ones(n,1);
end

function R_k = cov_obs_noise(~)
    m = 1; %Dar la dimension del espacio de mediciones
    sigma = 0.01;
    R_k = sigma*eye(m);
end

function v_k = obs_noise(k)
    %Devuelve el ruido para la medicion en el instante k
    m = 1; %Dar la dimension del espacio de estados
    v_k = mvnrnd(zeros(m,1),cov_obs_noise(k),1); 
end

function w_k = dyn_noise(~)
    %noise in the dynamics at time k
    n=2;
    %Q = 0.00*eye(n);
    T = 0.1;
    a= 0.2;
    Q =a^2 * [(T^4)/4 (T^3)/2;(T^3)/2 T^2];
    w_k = mvnrnd(zeros(n,1),Q,1); 
end
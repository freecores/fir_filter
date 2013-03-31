clear all;
format long g;

%%%%%%%%%%%%%%%%%%%% h[n] %%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%               OBTENCION DE COEFICIENTES PARA
%                       filtro_FIR.vhd
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

fs = 32e6
orden = 31                  % = coeficientes -1
fc = 10e3                  % frec de corte

bits_entrada = 8

% NOTA: las ilustraciones que se muestran tienen los valores 
%       redondeados por lo que para verlos correctamente (en la 
%       linea de comandos) descomentar el penultimo parrafo al 
%       final de este archivo (y_u, y_delta, y_rampas)

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


wn = 2*fc/fs

h=fir1(orden,wn) 			% coeficientes del filtro FIR

stem(h);
title('h(n) = coeficientes filtro FIR');grid on;

%%%%%%%%%%%%%%%%%%%% x[n] %%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);

num_muestras=64;
magnitud=2^bits_entrada-1;            % amplitud de las entradas 0:1:255 (A)
% entrada impulsional delta[n]

delta=[];
delta(1:num_muestras)=zeros(1);
delta(num_muestras/2)=ones(1)*magnitud;
subplot(2,3,1); axis([0 64 0 max(delta)]); stem(delta); title('x[n]=A delta[n]'); grid on;

% ventana rectangular u[n]

u=[];
u(1:num_muestras)=zeros(1);
u(num_muestras/4 :num_muestras/4+7)=ones(1)*magnitud;
subplot(2,3,2); axis([0 64 0 max(u)]); stem(u); title('x[n]=A u[n]'); grid on;

% rampa rampas[n]

rampas=[];
rampas(1:num_muestras)=zeros(1);
for R = 1:(num_muestras/4)
  rampas(R)=magnitud/(num_muestras/4)*R;
end
for R = (num_muestras/4+1):num_muestras/2
  rampas(R)=magnitud/(num_muestras/4)*(R-(num_muestras/4));
end
subplot(2,3,3); axis([0 64 0 max(rampas)]); stem(rampas); title('x[n]=A rampas[n]'); grid on;

%%%%%%%%%%%%%%%%%%%%% y[n] %%%%%%%%%%%%%%%%%%%%%%%%

y_u=filter(h,1,u);			% respuesta a: A u[n]
y_delta=filter(h,1,delta); 	% respuesta a: A delta[n]
y_rampas=filter(h,1,rampas);% respuesta a: A rampas[n]

subplot(2,3,4); axis([0 64 0 max(y_delta)]); stem(y_delta); title('y[n]'); grid on;
subplot(2,3,5); axis([0 64 0 max(y_u)]); stem(y_u); title('y[n]'); grid on;
subplot(2,3,6); axis([0 64 0 max(y_rampas)]); stem(y_rampas); title('y[n]'); grid on;

format;



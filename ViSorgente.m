function [U]=ViSorgente(Centro, Estremo_1, Estremo_2, L2G_TransfMatrix, G2L_TransfMatrix)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Questa funzione permette di calcolare il vettore di velocità indotta da
% % una linea di sorgenti. 
% % Input:    Centro: Vettore colonna di coordinate del punto in cui va
% %                   calcolata la velocità
% %           Estremo_1/2: vettori colonna di coordinate degli estremi del
% %                        pannello di sorgenti
% %           L2G_TransfMatrix: Matrice 2x2 di trasformazione da coordinate
% %                             locali a coordinate globali. Dovete
% %                             calcolarla voi esternamente.
% %           G2L_TransfMatrix: Matrice 2x2 di trasformazione da coordinate
% %                             globali a coordinate locali. Dovete
% %                             calcolarla voi esternamente.
% % Output:   U : Vettore colonna contenente le componenti x ed y della
% %               velocità indotta.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if nargin == 5
%     is2Print = false;
% end

% Trasformo da coordinate globali a coordinate locali
Centro=G2L_TransfMatrix{1,1}.*Centro{1,1};

Estremo_1=G2L_TransfMatrix{1,1}.*Estremo_1{1,1};

Estremo_2=G2L_TransfMatrix{1,1}.*Estremo_2{1,1};

%% Calcolo u e v in coordinate locali
% Calcolo r1 (congiungente punto indotto - estremo 1)
r1=Centro-Estremo_1;

% Calcolo theta_1 (angolo che la congiungente r1 forma con l'asse x locale)
theta_1=atan2(r1(2), r1(1));

% Calcolo r2 (congiungente punto indotto - estremo 2)
r2=Centro-Estremo_2;

% Calcolo theta_2 (angolo che la congiungente r2 forma con l'asse x locale)
theta_2=atan2(r2(2), r2(1));

% Fix in caso di auto-induzione
if (abs(theta_1)<10^(-12) && abs(theta_2)>3); theta_1=0; theta_2=pi; end
if (abs(theta_2)<10^(-12) && abs(theta_1)>3); theta_2=0; theta_1=-pi; end

% Calcolo le componenti della velocità
u=-(0.5/pi)*log(norm(r2)/norm(r1));
v=theta_2-theta_1;
v=v/(2*pi);

%% Converto da coordinate locali a coordinate globali
U=L2G_TransfMatrix*[u;v];

if abs(U(1))<10^(-12)
    U(1)=0;
end
if abs(U(2))<10^(-12)
    U(2)=0;
end

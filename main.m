%% Test Hess Smith
clc
close all
clear

%% Input

U_inf=1;  % Velocità all'infinito [m/s]
alpha=0;   % Angolo di incidenza [°]
U_inf_x=U_inf*cos(deg2rad(alpha));
U_inf_y=U_inf*sin(deg2rad(alpha));

U_inf=[U_inf_x; U_inf_y];
U_inf_normal=[-U_inf(2); U_inf(1)];
U_inf_normal=U_inf_normal./norm(U_inf_normal);

TestCase=0;
NCorpi=1;  % Numero di corpi da analizzare

CodiceProfilo=cell(NCorpi,1);
CodiceProfilo{1}='6412';
Chord=[1];
NPannelli=[102];

LE_X_Position=[0];
LE_Y_Position=[0];

%% Creazione profilo
% Numero profilo:
i=1;

[x,y]=CreateProfile(CodiceProfilo{i},NPannelli(i),Chord(i));

Corpi{i}.x=x;
Corpi{i}.y=y;

%% Creazione di una struttura di pannelli
Centro=cell(NCorpi,1);
Normale=cell(NCorpi,1);
Tangente=cell(NCorpi,1);
Estremo_1=cell(NCorpi,1);
Estremo_2=cell(NCorpi,1);
alpha=cell(NCorpi,1);
lunghezza=cell(NCorpi,1);
L2G_TransfMatrix=cell(NCorpi,1);
G2L_TransfMatrix=cell(NCorpi,1);
for i=1:NCorpi
    [Centro{i},Normale{i},Tangente{i},Estremo_1{i},Estremo_2{i},alpha{i},lunghezza{i},L2G_TransfMatrix{i},G2L_TransfMatrix{i}]=CreaStrutturaPannelli(Corpi{i});
end
        
figure (1)    %Corpo_i
plot(Centro{i,i}(:,1),Centro{1,1}(:,2),'-o','LineWidth',1)
title(['PANNELLIZZAZIONE NACA',CodiceProfilo{1}]);
grid on
axis equal

%% Inizializzazione matrici e vettori
% Ora che ho i pannelli, posso inizializzare la matrice ed i vettori

NCols=sum(NPannelli)+NCorpi;
NRows=NCols;
matrixA=zeros(NRows,NCols);
TermineNoto=zeros(NRows,1);

%%
us=ViSorgente(Centro, Estremo_1, Estremo_2, L2G_TransfMatrix, G2L_TransfMatrix);
uv=ViVortice(Centro, Estremo_1, Estremo_2, L2G_TransfMatrix, G2L_TransfMatrix);

matrixAs=dot(us,Normale);
vectorcs=dot(us,Tangente(1))+dot(us,Tangente(end));
%vectorav=sum;
cv=sum(vectorcs);

bv=-U_inf_normal*(Tangente(1)+Tangente(end));
bs=dot(-U_inf,normale);

matrixA=[];
TermineNoto=[];

%%




%% Risoluzione sistema lineare
Soluzione=linsolve(matrixA,TermineNoto);


%% Calcolo del cp e della velocità sui pannelli
U_Pannelli=cell(NCorpi, 1);
Ut_Pannelli=cell(NCorpi, 1);
Un_Pannelli=cell(NCorpi, 1);
Cp=cell(NCorpi, 1);

for Corpo_i=1:NCorpi
    U_Pannelli{Corpo_i}=zeros(NPannelli(Corpo_i),2);
    Ut_Pannelli{Corpo_i}=zeros(NPannelli(Corpo_i),1);
    Un_Pannelli{Corpo_i}=zeros(NPannelli(Corpo_i),1);
end

for Corpo_i=1:NCorpi
    for i=1:NPannelli(Corpo_i)

        U_Pannelli{Corpo_i}(i,:)=U_inf'; 
        Centro_qui = Centro{Corpo_i}(i,:)';
        Tangente_qui=Tangente{Corpo_i}(i,:)'; 
        Normale_qui=Normale{Corpo_i}(i,:)'; 
    
        for Corpo_j=1:NCorpi
            for j=1:NPannelli(Corpo_j)
                Estremo_1_qui=Estremo_1{Corpo_j}(j,:)';             
                Estremo_2_qui=Estremo_2{Corpo_j}(j,:)';
                L2G_TransfMatrix_qui=squeeze(L2G_TransfMatrix{Corpo_j}(j,:,:));
                G2L_TransfMatrix_qui=squeeze(G2L_TransfMatrix{Corpo_j}(j,:,:));

                U_Sorgente=ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
                U_Vortice=ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);

                U_Pannelli{Corpo_i}(i,:)=U_Pannelli{Corpo_i}(i,:)+sigma_mia{Corpo_j}(j).* U_Sorgente'+gamma_mia(Corpo_j).*U_Vortice';
            end
        end
        Ut_Pannelli{Corpo_i}(i) = dot(U_Pannelli{Corpo_i}(i, :)', Tangente_qui);
        Un_Pannelli{Corpo_i}(i) = dot(U_Pannelli{Corpo_i}(i, :)', Normale_qui);
    end
    Cp{Corpo_i} = 1-Ut_Pannelli{Corpo_i}.^2/norm(U_inf)^2;
end


Cl = cell(NCorpi, 1);
for Corpo_i = 1:NCorpi
    Cl_qui = 0;
    for i =1:NPannelli(Corpo_i)

        Normale_qui = Normale{Corpo_i}(i, :)';
        lunghezza_qui = lunghezza{Corpo_i}(i);

        Cl_qui = Cl_qui + (-Cp{Corpo_i}(i)*( lunghezza_qui.*dot(Normale_qui, U_inf_normal )));
    end
    
    Cl{Corpo_i} = Cl_qui / Chord(Corpo_i);
end

%%


figure
legend_string = cell(NCorpi, 1);
for Corpo_i = 1:NCorpi
%     plot(Centro{Corpo_i}(:, 1) - min(Centro{Corpo_i}(:, 1)), -Cp{Corpo_i}, '-')
    plot(Centro{Corpo_i}(:, 1), -Cp{Corpo_i}, '-')
    hold on
    legend_string{Corpo_i} = strcat("Corpo ", num2str(Corpo_i));
end
title('$C_P$', 'interpreter', 'latex')

% Profilo_1 = importdata("Profilo_1.dat");
% Profilo_2 = importdata("Profilo_2.dat");
% 
% plot(Profilo_1(:, 1)-min(Profilo_1(:, 1)), -Profilo_1(:, 2), '*')
% plot(Profilo_2(:, 1)-min(Profilo_2(:, 1)), -Profilo_2(:, 2), '*')
% legend_string{3} = strcat("Ale, Corpo ", num2str(1));
% legend_string{4} = strcat("Ale, Corpo ", num2str(2));
legend(legend_string)



figure
legend_string = cell(NCorpi, 1);
for Corpo_i = 1:NCorpi
%     plot(Centro{Corpo_i}(:, 1) - min(Centro{Corpo_i}(:, 1)), -Cp{Corpo_i}, '-')
    plot(Corpi{Corpo_i}.x, Corpi{Corpo_i}.y, '-')
    hold on
    legend_string{Corpo_i} = strcat("Corpo ", num2str(Corpo_i));
end
axis equal
legend(legend_string, 'interpreter', 'latex')


%%


ifSaveFigures=false;

if ifSaveFigures

    xMin = 2000;
    xMax = -2000;
    yMin = 2000;
    yMax = -2000;
    
    for Corpo_i = 1:NCorpi
        
        xMin = min(xMin, min(Centro{Corpo_i}(:, 1)));
        xMax = max(xMax, max(Centro{Corpo_i}(:, 1)));
        
        yMin = min(yMin, min(Centro{Corpo_i}(:, 2)));
        yMax = max(yMax, max(Centro{Corpo_i}(:, 2)));
    end
    
    xMin = xMin - 1;
    xMax = xMax + 1;
    yMin = yMin - 1;
    yMax = yMax + 1;
    



    
    Nx = 400;
    Ny = 400;
    
    
    x = linspace(xMin,xMax, Nx);
    y = linspace(yMin, yMax, Ny);

    [X, Y] = meshgrid(x, y);

    isIn = zeros(Nx, Ny);
    
    t = cputime;
    parfor(i = 1:Nx)
%         i
        for(j = 1:Ny)
            for( Corpo_i = 1:NCorpi)
                Boundary = [ Corpi{Corpo_i}.x Corpi{Corpo_i}.y];
                if (inpolygon(X(i, j), Y(i, j), Boundary(:, 1), Boundary(:, 2)))
                    isIn(i, j) = 1;
                end
            end
        end
    end
    cputime - t


    U_Mesh = zeros(Nx, Ny);
    V_Mesh = zeros(Nx, Ny);
    U_Mesh_Mag = zeros(Nx, Ny);
    Cp_Mesh = zeros(Nx, Ny);

    t = cputime;

        parfor PointIndex_i = 1:Nx
    %         PointIndex_i
            for PointIndex_j = 1:Ny

                if(~isIn(PointIndex_i, PointIndex_j))

                    U = U_inf'; 
                    Centro_qui = [X(PointIndex_i, PointIndex_j); Y(PointIndex_i, PointIndex_j)];

                    for Corpo_j = 1:NCorpi
                        for j = 1:NPannelli(Corpo_j)

                            Estremo_1_qui = Estremo_1{Corpo_j}(j, :)';             
                            Estremo_2_qui = Estremo_2{Corpo_j}(j, :)';
                            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix{Corpo_j}(j, :, :));
                            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix{Corpo_j}(j, :, :));

                            U_Sorgente = ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
                            U_Vortice = ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);

                            U = U + sigma_mia{Corpo_j}(j) .* U_Sorgente' + gamma_mia(Corpo_j) .* U_Vortice';

                        end
                    end

                    U_Mesh(PointIndex_i, PointIndex_j) = U(1);
                    V_Mesh(PointIndex_i, PointIndex_j) = U(2);
                    U_Mesh_Mag(PointIndex_i, PointIndex_j) = norm(U);
                    Cp_Mesh(PointIndex_i, PointIndex_j) = 1 - norm(U) / norm(U_inf);

                else
                    U_Mesh(PointIndex_i, PointIndex_j) = NaN;
                    V_Mesh(PointIndex_i, PointIndex_j) = NaN;
                    U_Mesh_Mag(PointIndex_i, PointIndex_j) = NaN;
                    Cp_Mesh(PointIndex_i, PointIndex_j) = NaN;
                end
            end
        end


    SavingNameStart = "./Figures/test_";


    UMag_fig = figure;
    contourf(X, Y, U_Mesh_Mag, 100,'LineStyle','None');
    colorbar
    hold on
    for Corpo_j = 1:NCorpi
        plot(Corpi{Corpo_j}.x, Corpi{Corpo_j}.y, '-k')
    end
    axis equal
    xlabel('x', 'interpreter', 'latex')
    ylabel('y', 'interpreter', 'latex')
%     legend("$U_{MAG}$", 'interpreter', 'latex')
    title('Contorno di $U_{MAG}$', 'interpreter', 'latex')
    SavingName = strcat(SavingNameStart, '_UMag.eps');
    saveas(UMag_fig, SavingName);
    
    
    U_fig = figure;
    contourf(X, Y, U_Mesh, 100,'LineStyle','None');
    colorbar
    hold on
    for Corpo_j = 1:NCorpi
        plot(Corpi{Corpo_j}.x, Corpi{Corpo_j}.y, '-k')
    end
    axis equal
    xlabel('x', 'interpreter', 'latex')
    ylabel('y', 'interpreter', 'latex')
    title('Contorno di $U$', 'interpreter', 'latex')
    SavingName = strcat(SavingNameStart, '_U.eps');
    saveas(U_fig, SavingName);
    
    
    V_fig = figure;
    contourf(X, Y, V_Mesh, 1000,'LineStyle','None');
    colorbar
    hold on
    for Corpo_j = 1:NCorpi
        plot(Corpi{Corpo_j}.x, Corpi{Corpo_j}.y, '-k')
    end
    axis equal
    xlabel('x', 'interpreter', 'latex')
    ylabel('y', 'interpreter', 'latex')
    title('Contorno di $V$', 'interpreter', 'latex')
    SavingName = strcat(SavingNameStart, '_V.eps');
    saveas(V_fig, SavingName);


    Cp_fig = figure;
    contourf(X, Y, Cp_Mesh, 100,'LineStyle','None');
    colormap(flipud(hot));
    colorbar
    hold on
    for Corpo_j = 1:NCorpi
        plot(Corpi{Corpo_j}.x, Corpi{Corpo_j}.y, '-k')
    end
    axis equal
    xlabel('x', 'interpreter', 'latex')
    ylabel('y', 'interpreter', 'latex')
    title('Contorno di $U_{MAG}$', 'interpreter', 'latex')
    SavingName = strcat(SavingNameStart, '_Cp.eps');
    saveas(Cp_fig, SavingName);
    
    
    Streamlines_fig = figure;
    contourf(X, Y, U_Mesh_Mag, 100,'LineStyle','None');
    colormap(flipud(hot));
    colorbar
    hold on
    streamslice(X, Y, U_Mesh, V_Mesh, 10);
    for Corpo_j = 1:NCorpi
        plot(Corpi{Corpo_j}.x, Corpi{Corpo_j}.y, '-k')
    end
    axis equal
    xlabel('x', 'interpreter', 'latex')
    ylabel('y', 'interpreter', 'latex')
    title('Contorno di $U_{MAG}$ e linee di corrente', 'interpreter', 'latex')
    SavingName = strcat(SavingNameStart, '_Streamlines.eps');
    saveas(Streamlines_fig, SavingName);
    
end
    
    

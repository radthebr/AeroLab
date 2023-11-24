%% 
clc, close all
format long g

%%
path_run_script='data/com1foil.txt';
naca='NACA 0012';

AoA=[0, 5, 10];

%% Wtite exe Script
nangles=Length(AoA);
xfoil(nangles)=struct();

for iangle=1:nangles
    polar_file=['data/polar_',num2str(iangle),'.txt'];
    cp_file=['data/cp_',num2str(iangle),'.txt'];

    fid=fopen(path_run_script,'w');
    printf(fid, '%s\n',' ');
    printf(fid, '%s\n', naca);
    printf(fid, '%s\n','OPER');
    printf(fid, '%s\n','PACC');
    printf(fid, '%s\n', polar_file);
    printf(fid, '%s\n',' ');
    printf(fid, '%s\n', ['ALFA ', num2str(AoA(iangle))]);
    printf(fid, '%s\n', ['CPWR ', cp_file]);
    printf(fid, '%s\n',' ');
    printf(fid, '%s\n',' ');
    printf(fid, '%s\n','QUIT');
    fclose(fid);

    system('xfoil <data/com1foil.txt>');
    %pause(0.1)
    data=importdata(cp_file,' ', 1);
    xfoil(iangle).cp=data.data;
    data=importdata(polar_file, ' ', 12);
    xfoil(iangle).polar=data.data;

    %system('rm data/*');
end

%% Execute Script
figure
plot(xfoil(1).cp(:,1), -xfoil(1).cp(:,2));
hold on;
plot(xfoil(2).cp(:,1), -xfoil(2).cp(:,2));
plot(xfoil(3).cp(:,1), -xfoil(3).cp(:,2));
hold off

%% Read Data

%% Clean Data





















function [x,y]=CreateProfile(Profilo,NPannelli,Chord)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

        % Sfruttiamo XFoil per la creazione del profilo, sperando che sia
        % all'interno del database
        %Profilo=char(Profilo);
        fileID1=fopen('launch.txt','w');
        fprintf(fileID1,'open -a XQuartz\n');
        fprintf(fileID1,'export DISPLAY=:0.0\n');
        fprintf(fileID1,'cd /Applications/Xfoil-for-Mac-main/runs\n');
        fprintf(fileID1,'../bin/xfoil < XFoilInput.txt\n');
        fclose(fileID1);

        fileID2=fopen('XFoilInput.txt','w');
        fprintf(fileID2,'NACA %s\n',Profilo);
        fprintf(fileID2,'pane\n');
        fprintf(fileID2,'gdes\n');
        fprintf(fileID2,'tgap 0 0\n');
        fprintf(fileID2,'exec\n\n');
        fprintf(fileID2,'ppar\n');
        fprintf(fileID2, ['n ',num2str(NPannelli+1),'\n\n\n']);
        filename = strcat('NACA_',Profilo,'.dat');
        fprintf(fileID2,'save %s\n',filename);
        fprintf(fileID2,'Y\n');
        fprintf(fileID2,'quit\n');
        fclose(fileID2);
        
        system('mv /Users/lafrazz/Locale/XfoilEXAMPLES/Laboratorio/XFoilInput.txt /Applications/Xfoil-for-Mac-main/runs/XFoilInput.txt');
        system('sh /Users/lafrazz/Locale/XfoilEXAMPLES/Laboratorio/launch.txt');
        system(['mv /Applications/Xfoil-for-Mac-main/runs/NACA_',Profilo,'.dat /Users/lafrazz/Locale/XfoilEXAMPLES/Laboratorio/NACA_',Profilo,'.dat']);
        
        Corpo=ImportXFoilProfile(filename);
        
        % Prima flippa i vettori
        x=flipud(Corpo.x);
        y=flipud(Corpo.y);
        
        x=x.*Chord;
        y=y.*Chord;
        
        system('rm /Users/lafrazz/Locale/XfoilEXAMPLES/Laboratorio/launch.txt');
        system(['rm /Users/lafrazz/Locale/XfoilEXAMPLES/Laboratorio/NACA_',Profilo,'.dat']);
end
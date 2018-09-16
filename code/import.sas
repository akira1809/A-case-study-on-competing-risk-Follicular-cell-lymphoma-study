PROC IMPORT OUT= WORK.follicular 
            DATAFILE= "C:\akira\school\KUMC\2018 Spring\BIOS845\Project\
datasets\follic.txt" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

function [Nr,Nrx,Nry]=RIS_UPAMapping(Kr)
    if Kr==1
        Nr=240;Nrx=15;Nry=16;
    elseif Kr==2
        Nr=120;Nrx=10;Nry=12;
    elseif Kr==3
        Nr=80;Nrx=8;Nry=10;
    elseif Kr==4
        Nr=60;Nrx=6;Nry=10;
    elseif Kr==5
        Nr=48;Nrx=6;Nry=8;
    elseif Kr==6
        Nr=40;Nrx=5;Nry=8;
    else
        fprintf("Error: Undefined RIS UPA configuration!");
    end
end
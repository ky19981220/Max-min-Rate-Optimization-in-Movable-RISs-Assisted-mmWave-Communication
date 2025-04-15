function [POS_UE_mc]=getUEPOS(Niter,maxKu,radius,center,UEheight)
    POS_UE_mc=cell(Niter,1);
    for ii=1:Niter
        num=0;
        while(1)
            pos_rand=2*rand(1,2)-1;
            if norm(pos_rand)<=1
                num=num+1;
                POS_UE_mc{ii}=[POS_UE_mc{ii};center+[radius*pos_rand,UEheight]];
            end
            if num==maxKu
                break;
            end
        end
    end
end
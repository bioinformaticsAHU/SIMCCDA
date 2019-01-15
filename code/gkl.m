function [result_circ,result_dis]=gkl(nc,nd,inter)
%gkl compute Gaussian interaction profile kernel
%	Inputs:
%		nc: the number of circRNAs
%		nd:	the number of diseases
%		inter: an nc*nd association matrix between circRNAs and diseases
%
%	Outputs:
%		result_circ: Gaussian interaction profile kernel of circRNAs
%		result_dis: Gaussian interaction profile kernel of diseases

    for i=1:nc
        sl(i)=norm(inter(i,:))^2;   %Maximum singular value
    end
    gamal=nc/sum(sl')*1;
    [pkl]=zeros(nc,nc);
    for i=1:nc
        for j=1:nc
            pkl(i,j)=exp(-gamal*(norm(inter(i,:)-inter(j,:)))^2);
        end
    end        
    for i=1:nd
        sd(i)=norm(inter(:,i))^2;
    end
    gamad=nd/sum(sd')*1; 
    for i=1:nd
        for j=1:nd
            pkd(i,j)=exp(-gamad*(norm(inter(:,i)-inter(:,j)))^2);
        end
    end 
    result_circ=pkl;
    result_dis=pkd;
    
end
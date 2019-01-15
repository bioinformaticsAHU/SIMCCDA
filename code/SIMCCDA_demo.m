function [M_recover]=SIMCCDA_demo(cD,circsim,dissim,a,b)
%% configuration
addpath('SIMC');

%% compute Gaussian interaction profile kernel of circRNAs and diseases
[nc,nd]=size(cD);
[CC,DD]=gkl(nc,nd,cD);

%% extract primary feature vectors of circRNAs and diseases
circ_feature=PCA((circsim+CC)/2,a);
dis_feature=PCA((dissim+DD)/2,b);

%% inductive matrix completion
Omega=find(cD==1);
M_recover=SIMC(cD,Omega,circ_feature,dis_feature);   
end

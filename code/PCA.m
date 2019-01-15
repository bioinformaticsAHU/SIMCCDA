function [feature_vecs]=PCA(simMat,para)
%PCA extracting primary feature vectors based on energy strategy 
%	Inputs:
%		simMat: kernel of circRNAs or diseases
%		para:	percent of energy for extracting primary feature vectors
%	Outputs:
%		feature_vecs: primary feature vectors of circRNAs or diseases

    pca_rank=0;
    singular_mat=svd(simMat);
    for i=1:rank(simMat)
        if sum(singular_mat(1:i))>para*sum(svd(simMat))
            pca_rank=i;
            break;
        end
    end   
    [~,~,V]=svds(simMat,pca_rank);
    feature_vecs=V;     
end
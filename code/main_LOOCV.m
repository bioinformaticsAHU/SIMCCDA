clc;
clear;
%%
% SIMCCDA: Prediction of circRNA-disease associations based on inductive matrix completion
%  
% chr_diseasematrix.csv: an n*m association matrix between circRNAs and diseases, n is the number of circRNAs, and m is the number of diseases
% seqsimilarity.csv: an n*n sequence similarity matrix of circRNAs
% dissimilarity.csv: an m*m semantic similarity matrix of diseases

%% load data
cD=importdata('../dataset/Dataset-1/chr_diseasematrix.csv'); 
circsim=importdata('../dataset/Dataset-1/seqsimilarity.csv');  
dissim=csvread('../dataset/Dataset-1/dissimilarity.csv',1);
%% LOOCV
A_ori = cD;
[row,column]=size(A_ori);
[score_ori] = zeros(row,column);
for i=1:column 
    i
    [test]=zeros(row,column);
    for j=1:row
        if cD(j,i)==1           
            cD(j,i)=0;       
            [result]= SIMCCDA_demo(cD,circsim,dissim,0.6,0.9);
            score_ori(j,i) = result(j,i);
            cD=A_ori; 
            for m=1:row
                if cD(m,i)==0
                    test(m,i)=max(result(m,i),score_ori(m,i)); 
                end
            end           
        end 
    end
    for kk=1:row
        if cD(kk,i)==0
            score_ori(kk,i)=test(kk,i);
        end
    end
end

%% compute AUC
pre_label_score = score_ori(:);
label_y = A_ori(:);
[auc,topi]=roc(pre_label_score,label_y);




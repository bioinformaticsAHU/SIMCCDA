function [auc,topi] = roc(deci,label_y)

[~,ind] = sort(deci,'descend');
roc_y = label_y(ind);
topi=[sum(roc_y(1:10) == 1),sum(roc_y(1:30) == 1),sum(roc_y(1:50) == 1),sum(roc_y(1:100) == 1)];

stack_x = cumsum(roc_y == 0)/sum(roc_y == 0);
stack_y = cumsum(roc_y == 1)/sum(roc_y == 1);

auc=sum((stack_x(2:length(roc_y))-stack_x(1:length(roc_y)-1)).*stack_y(2:length(roc_y)));
plot(stack_x,stack_y,'--')
xlabel('1-Specificity');
ylabel('Sensitivity');
title(['ROC curve of (AUC = ' num2str(auc) ' )'])
end
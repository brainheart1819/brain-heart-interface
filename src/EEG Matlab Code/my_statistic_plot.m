%% function used for statistic data plot
function my_statistic_plot(data, number, Title, xLabels, legendLabels)
b = bar(data);
color = ["red","green","blue","yellow"];
bSet = [];
for i = 1:number/2
    b(2*i-1).FaceColor = char(color(i));
    b(2*i).FaceColor = char(color(i));
    bSet = [bSet b(2*i-1)];
end
set(gca,'XTickLabel',xLabels);
title(Title);
legend(bSet, legendLabels);

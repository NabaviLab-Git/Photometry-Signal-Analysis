function [xmat,ymat] = scatterjittered(x,xpos,width);
%scatterjittered: Add ordered jitter to scatterplot
[counts,edges,bin] = histcounts(x,'BinMethod','fd');
offset =   width/(max(counts));

if mod(max(counts),2) == 0
    ymat = NaN(length(counts),max(counts)+1);
    xmat = linspace(xpos-0.5*width,xpos+0.5*width,max(counts)+1);
else
    ymat = NaN(length(counts),max(counts));
    xmat = linspace(xpos-0.5*width,xpos+0.5*width,max(counts));
end

for i = 1:length(counts)
    ymat(i,1:counts(i)) = x(bin == i);
    ymat(i,:) = circshift(ymat(i,:),round(0.5*max(counts)-0.5*sum(~isnan(ymat(i,:)))));
    
end

end
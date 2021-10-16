function [b,stats]=plot_fit(x,y,robust,dots,line)
% usage: [b,stats]=plot_fit(x,y,0,'ko','r')

if robust==1
    [b,stats]=robustfit(x,y);
else
    [b,~,stats]=glmfit(x,y);
end

yFitted=b(1)+b(2)*x;

plot(x(:,end),y,dots,x(:,end),yFitted,line,'MarkerSize',24);



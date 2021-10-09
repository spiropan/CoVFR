function [b,stats]=plot_fit(x,y,robust,dots,line)
% usage: [b,stats]=plot_fit(x,y,0,'ko','r')
% pass empty b if want to estimate fit from scratch

if robust==1
    [b,stats]=robustfit(x,y);
else
    [b,~,stats]=glmfit(x,y);
end

yFitted=b(1)+b(2)*x;

%yFitted=glmval(b,x,'identity')
%x(:,end)

%plot(x,y,colstr,x,i2,col,'MarkerSize',7);

%plot(xFitting,yFitted,col,x,y,'g.','MarkerSize',2);
% THE OLD ONE BELOW
%plot(x,y,dots, xFitting,yFitted,line,'MarkerSize',8);
%plot(x,y,dots,x,yFitted,line,'MarkerSize',8);
%plot(x(:,end),y,dots,x(:,end),yFitted,line,'MarkerSize',8);
plot(x(:,end),y,dots,x(:,end),yFitted,line,'MarkerSize',8);



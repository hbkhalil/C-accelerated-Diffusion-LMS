function [] = print_network( xy,A,options )
%This function prints the network 
%   inputs:
%       xy:         Nodes coordinates (two column matrix)
%       A:          Topology matrix
%       options:    Display option:
%                                  options.EdgeColor        (default='k')
%                                  options.MarkerColor      (default=[10 222 208]/255)
%                                  options.Markersize       (default=24)
%                                  options.LineWidth        (default=1.5)
%                                  options.Fontsize         (default=10)


%%
x=xy(:,1);
y=xy(:,2);

options=initoptions('print_network',options);


options.MarkerEdgeColor='k';
options.MarkerFaceColor=[10 222 208]/255;

Fontsize=options.Fontsize;


N=size(A,1);

g=graph(A-eye(N));
p=plot(g,'NodeLabel',{});
text(xy(:,1),xy(:,2),cellstr(num2str([1:N]','%d'))...
    ,'HorizontalAlignment','center','Fontsize',Fontsize,'Interpreter', 'latex');

%axis([min(xy(:,1))-0.1 max(xy(:,1))+0.1 min(xy(:,2))-0.1 max(xy(:,2))+0.1]);
%axis([0 1 0 1])
%set(gca,'position',[0 0 1 1],'units','normalized')

p.XData = xy(:,1);
p.YData = xy(:,2);
p.Marker = 'o';
p.NodeColor = options.MarkerColor;
p.MarkerSize= options.Markersize;
p.EdgeColor= options.EdgeColor;
p.LineWidth= options.LineWidth;
axis off


end


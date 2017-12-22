function adjustFigures(h)

a = 270;
img_width = a / 300*325;
img_height = a;

h.Children(1).Position(1) = 0.87;
h.Children(1).Position(3) = 0.03;
h.Children(2).Position(3) = 0.7;


set(h,'Position',[0,0,img_width,img_height]);



end
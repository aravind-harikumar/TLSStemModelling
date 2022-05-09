function mat_points = pointsincircle(n, center, radius, z)
    %require in input: n = number of points that are needed
    %                  center = center coordinates of the circle 
    %                  radius = radius of the circle
    %                  z = height of the slice
    %the output is a matrix containig the coordinates X, Y, Z of
    %the points that has been created
    angle = 2*pi*rand(n,1);
    r = radius*sqrt(rand(n,1));
    X = r.*cos(angle)+center(1);
    Y = r.*sin(angle)+center(2);
    Z = z*ones(n,1);
    mat_points = [X Y Z];
    %figure;
    %plot(X,Y,'.','MarkerSize',5,'Color',[0 0.5 0]);
end
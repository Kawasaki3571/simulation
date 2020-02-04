xrange = 0.4;
yrange = 0.4;
dx = xrange / 10;
dx = 0.017;
ix = round(xrange / dx);
jx = round(yrange / dx);
p1 = ones(ix, jx);

center_x = xrange/2;
center_y = yrange/2;
kei = xrange / 2;
for x1_s = 0 : dx : xrange
	y1_s = center_y + sqrt(kei^2 - (x1_s - center_x)^2);
    i1 = round(x1_s/dx) + 1;
    j1 = round(y1_s/dx) + 1;
    u1(i1, j1) = 0;
    v1(i1, j1) = 0;
    for jstair = j1 : jx + 1
    	p1(i1, jstair) = 0;
        %v1(i1, jstair) = 0;
    end
end
for x1_s = 0 : dx : xrange
	y1_s = center_y - sqrt(kei^2 - (x1_s - center_x)^2);
    i1 = round(x1_s/dx) + 1;
    j1 = round(y1_s/dx) + 1;
    %p1(i1, j1) = 0;
    %v1(i1, j1) = 0;
	for jstair = 1 : j1
        p1(i1, jstair) = 0;
        %v1(i1, jstair) = 0;
    end
end

x = 1 : ix;
y = 1 : jx;
imagesc(x, y, p1(x, y))
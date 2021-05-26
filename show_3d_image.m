function show_3d_image(x_grid, y_grid, z_grid, image)

cla
colorbar('off');
set(gca, 'YDir','normal');

half_max_value = max(image(:))*0.5;
% tick_space = 1;
% figure,
p = patch(isosurface(x_grid, y_grid, z_grid, image, half_max_value)); 
isonormals(x_grid, y_grid, z_grid, image, p);
p.FaceColor = 'blue'; p.EdgeColor = 'none';
camlight;lighting gouraud; view(37.5,30);
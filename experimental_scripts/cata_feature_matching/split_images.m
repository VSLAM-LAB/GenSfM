clear

r_outer = 1250
r_inner = 350
pp = [2190 1463]

im_path = 'images_test/';
out_path = 'images_split/';

files = dir([im_path '*.JPG'])

for file_k = 1:length(files)
filename = files(file_k).name;

im = imread([im_path filename]);

%%
im_warp = {};
for cube_idx = 1:8
phase = pi/4 * (cube_idx-1);
n = 100;
[pphi,rr] = meshgrid(linspace(0+phase,pi/2+phase,n),linspace(r_inner,r_outer,n));

x = pp(1) + rr .* cos(pphi);
y = pp(2) + rr .* sin(pphi);

k = round(r_outer * 2 / pi);
target_y = rr - r_inner;
target_x = (pphi-phase) * k;


xx = [x(:)'; y(:)'];
yy = [target_x(:)'; target_y(:)'];

n_pts = 10;
tform_method = @(a,b) fitgeotrans(a,b,'polynomial',4);

tform1 = tform_method(xx',yy');
[im_warp{cube_idx},RB] = imwarp(im,tform1,'cubic','OutputView',imref2d([r_outer-r_inner ceil(pi/2*k)]));

filename_out = [filename sprintf('_%d.png',cube_idx-1)];

fprintf('Writing %s...\n',filename_out);
imwrite(im_warp{cube_idx},[out_path filename_out],'png')

end

fprintf('params=%d,%d,%d,%d,%d,%d,%d\n',size(im,2),size(im,1),pp(1),pp(2),r_inner,r_outer,k)

end


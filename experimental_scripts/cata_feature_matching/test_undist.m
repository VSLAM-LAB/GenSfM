clear

im = imread('DSC_1276.JPG');

r_outer = 1250
r_inner = 250
r_inner = 350
imagesc(im)

sz = size(im);
pp = [sz(2)/2;sz(1)/2];
pp = [2190 1463]
imagesc(im)
hold on
plot(pp(1),pp(2),'r.')
plotc(pp,r_outer,100)
plotc(pp,r_inner,100)
hold off
%%
im_warp = {};
colors = 'rgbcymkw'
close all
for cube_idx = 1:8
phase = pi/4 * (cube_idx-1);
n = 100
[pphi,rr] = meshgrid(linspace(0+phase,pi/2+phase,n),linspace(r_inner,r_outer,n))

x = pp(1) + rr .* cos(pphi);
y = pp(2) + rr .* sin(pphi);

subaxis(3,3,1,'Margin',0,'Spacing',0,'Padding',0)
if cube_idx == 1
    imagesc(im)
    hold on
    plotc(pp,r_outer,100)
    plotc(pp,r_inner,100)
end
plot(x(1:7:end),y(1:7:end),[colors(cube_idx) '.'],'MarkerSize',1)
axis equal off

k = round(r_outer * 2 / pi);
% k = r_outer / kkk
target_y = rr - r_inner;
target_x = (pphi-phase) * k;



xx = [x(:)'; y(:)'];
yy = [target_x(:)'; target_y(:)'];
% yy(2,:) = max(yy(:))-yy(2,:);

n_pts = 10;
% tform_method = @(a,b) fitgeotrans(a,b,'lwm',n_pts);
% tform_method = @(a,b) fitgeotrans(a,b,'pwl');
tform_method = @(a,b) fitgeotrans(a,b,'polynomial',4);


tform1 = tform_method(xx',yy');

% im_warp = imwarp(im,tform1,'OutputView',imref2d(size(im)));
[im_warp{cube_idx},RB] = imwarp(im,tform1,'cubic','OutputView',imref2d([r_outer-r_inner ceil(pi/2*k)]));

%[max(target_y(:)) max(target_x(:))]
% figure(1+cube_idx)
subaxis(3,3,1+cube_idx,'Margin',0,'Spacing',0,'Padding',0)
imagesc(im_warp{cube_idx})
axis equal off
end

%% Transfer coordinate back

fprintf('params=%d,%d,%d,%d,%d\n',pp(1),pp(2),r_inner,r_outer,k)

im_k = 6
%x_undist = [1041.7; 844.09];
imagesc(im_warp{im_k});

x_undist = ginput(1)

r = x_undist(2) + r_inner;
phi = x_undist(1)/k + pi/4 * (im_k-1);

x_dist = pp + r * [cos(phi) sin(phi)];


figure(1)
subplot(1,2,1)
imagesc(im_warp{im_k});
hold on
plot(x_undist(1),x_undist(2),'r.')
hold off
subplot(1,2,2)
imagesc(im);
hold on
plot(x_dist(1),x_dist(2),'r.')
hold off




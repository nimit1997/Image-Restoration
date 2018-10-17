%Nimit Kapadia

kernel = imread('Kernel2G_crop.png');
%kernel_gray = rgb2gray(kernel)
kernel_crop = kernel(23:40,23:40)
imwrite(kernel_crop,'Kernel2G_center.png')
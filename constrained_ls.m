%Nimit Kapadia

psf_kernel = im2double(imread('kernel4G_center.png')); %reading the kernel image and bringing the array values between 0 and 1
%using im2double
norm_factor = sum(sum(psf_kernel)); %normalizing factor
psf_kernel = psf_kernel./norm_factor; %normalization
blurd = imread('VBlur.jpg'); %reading the blurred image
%R = psf_kernel(:,:,1) 

p = [0,-1,0;-1,4,-1;0,-1,0]; %laplacian operator
P = fftshift(fft2(p,PQ(1),PQ(2))); %dft of the laplacian operator
%P = dft2(p,PQ(1),PQ(2));

PQ = paddedsize(size(blurd)); %computing padded sizes useful for FFT-based filtering
H = fftshift(fft2(double(psf_kernel), PQ(1), PQ(2))); %dft of the kernel
%H = dft2(psf_kernel, PQ(1), PQ(2));
gammafactor = input('gamma = ');
gammafactor = typecast(gammafactor,'double');
while(gammafactor>-1)
    GR = fftshift(fft2(double(blurd(:,:,1)), PQ(1), PQ(2))); %dft of the blurred image (Red)
    %GR = dft2(double(blurd(:,:,1)), PQ(1), PQ(2));
    GG = fftshift(fft2(double(blurd(:,:,2)), PQ(1), PQ(2))); %dft of the blurred image (Green)
    %GG = dft2(double(blurd(:,:,2)), PQ(1), PQ(2));
    GB = fftshift(fft2(double(blurd(:,:,3)), PQ(1), PQ(2))); %dft of the blurred image (Blue)
    %GB = dft2(double(blurd(:,:,3)), PQ(1), PQ(2));
    F = zeros(PQ(1),PQ(2),3); %zero matrix
    F(:,:,1) = GR; %substituting dft values of red channel to the F matrix
    F(:,:,2) = GG; %substituting dft values of green channel to the F matrix
    F(:,:,3) = GB; %substituting dft values of blue channel to the F matrix
    
    %Constrained Least Square filtering formula
    HMODSQR = abs(H).*abs(H);
    PMODSQR = abs(P).*abs(P);
    DEN = bsxfun(@plus,HMODSQR,gammafactor*PMODSQR);
    G = conj(H)./DEN;
  
    %Inverse filtering using constrained least square filtering formula by putting gamma = 0
%     dummy = ones(PQ(1),PQ(2));
%     G = dummy./H;
    
    %constrained least square filtering on all three channels
    F(:,:,1) = G.*GR;
    F(:,:,2) = G.*GG;
    F(:,:,3) = G.*GB;
    
    %inverse dft on all three channels
    R_filtered = ifft2(fftshift(F(:,:,1)));
    %R_filtered = idft2(F(:,:,1));
    R_filtered = R_filtered(2:size(blurd,1)+1, 2:size(blurd,2)+1);
    G_filtered = ifft2(fftshift(F(:,:,2)));
    %G_filtered = idft2(F(:,:,2));
    G_filtered = G_filtered(2:size(blurd,1)+1, 2:size(blurd,2)+1);
    B_filtered = ifft2(fftshift(F(:,:,3)));
    %B_filtered = idft2(F(:,:,3));
    B_filtered = B_filtered(2:size(blurd,1)+1, 2:size(blurd,2)+1);

    %Display results (show all values)

    RGB_filtered = zeros(size(R_filtered,1),size(R_filtered,2),3,'uint8');

    RGB_filtered(:,:,1) = real(R_filtered);
    RGB_filtered(:,:,2) = real(G_filtered);
    RGB_filtered(:,:,3) = real(B_filtered);

    figure, imshow(RGB_filtered,[])
    gammafactor = input('gamma = ');
    
    %calculates peak signal to noise ratio
    %peaksnr = psnr(RGB_filtered,imread('GroundTruth2_1_1.jpg'))
    
    %calculates ssim value
    %ssimval = ssim(RGB_filtered,imread('GroundTruth2_1_1.jpg'))
end
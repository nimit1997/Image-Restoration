psf_kernel = im2double(imread('kernel1G_center.png'));
norm_factor = sum(sum(psf_kernel));
psf_kernel = psf_kernel./norm_factor;
blurd = imread('Blurry2_1.jpg');
%R = psf_kernel(:,:,1) 

p = [0,-1,0;-1,4,-1;0,-1,0];
P = fftshift(fft2(p,PQ(1),PQ(2)));
%P = dft2(p,PQ(1),PQ(2));

PQ = paddedsize(size(blurd));
H = fftshift(fft2(double(psf_kernel), PQ(1), PQ(2)));
%H = dft2(psf_kernel, PQ(1), PQ(2));
gammafactor = input('gamma = ');
gammafactor = typecast(gammafactor,'double');
while(gammafactor>-1)
    GR = fftshift(fft2(double(blurd(:,:,1)), PQ(1), PQ(2)));
    %GR = dft2(double(blurd(:,:,1)), PQ(1), PQ(2));
    GG = fftshift(fft2(double(blurd(:,:,2)), PQ(1), PQ(2)));
    %GG = dft2(double(blurd(:,:,2)), PQ(1), PQ(2));
    GB = fftshift(fft2(double(blurd(:,:,3)), PQ(1), PQ(2)));
    %GB = dft2(double(blurd(:,:,3)), PQ(1), PQ(2));
    F = zeros(PQ(1),PQ(2),3);
    F(:,:,1) = GR;
    F(:,:,2) = GG;
    F(:,:,3) = GB;
    
    HMODSQR = abs(H).*abs(H);
    PMODSQR = abs(P).*abs(P);
    DEN = bsxfun(@plus,HMODSQR,gammafactor*PMODSQR);
    G = conj(H)./DEN;
  
%     dummy = ones(PQ(1),PQ(2));
%     G = dummy./H;
    
    F(:,:,1) = G.*GR;
    F(:,:,2) = G.*GG;
    F(:,:,3) = G.*GB;
    
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
    peaksnr = psnr(RGB_filtered,imread('GroundTruth2_1_1.jpg'))
    ssimval = ssim(RGB_filtered,imread('GroundTruth2_1_1.jpg'))
end
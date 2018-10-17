%Nimit Kapadia

psf_kernel = im2double(imread('Kernel4G_center.png')); %reading the kernel image and bringing the array values between 0 and 1
%using im2double
%disp(max(transpose(card_orig)))
norm_factor = sum(sum(psf_kernel)); %normalizing factor
psf_kernel = psf_kernel./norm_factor; %normalization
blurd = imread('VBlur.jpg'); %reading the blurred image
%R = psf_kernel(:,:,1) 

PQ = paddedsize(size(blurd)); %computing padded sizes useful for FFT-based filtering
H = fftshift(fft2(psf_kernel, PQ(1), PQ(2))); %dft of the kernel
%H = dft2(psf_kernel, PQ(1), PQ(2));
RO = input('Ro = ');
RO = typecast(RO,'double');
while(RO>0)
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
    for R = 1:PQ(1)
        for C=1:PQ(2)
            Rdash = R - PQ(1)/2;
            Cdash = C - PQ(2)/2;
            if (Rdash^2 + Cdash^2)<(RO^2) %applying the constraint for truncation
            %if PQ(1)<0
                %inverse filtering on all three channels
                F(R,C,1) = GR(R,C)/H(R,C); 
                F(R,C,2) = GG(R,C)/H(R,C);
                F(R,C,3) = GB(R,C)/H(R,C);
%             else
                  %everything 0 beyond the truncation radius
%                 F(R,C,1) = 0;
%                 F(R,C,2) = 0;
%                 F(R,C,3) = 0;
            end
        end
    end
    
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
    RO = input('Ro = ');
    %calculates peak signal to noise ratio
    %peaksnr = psnr(RGB_filtered,imread('GroundTruth1_1_1.jpg'))
    
    %calculates ssim value
    %ssimval = ssim(RGB_filtered,imread('GroundTruth1_1_1.jpg'))
end
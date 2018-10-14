psf_kernel = im2double(imread('Kernel1G_center.png'));
%disp(max(transpose(card_orig)))
norm_factor = sum(sum(psf_kernel));
psf_kernel = psf_kernel./norm_factor;
blurd = imread('Blurry1_1.jpg');
%R = psf_kernel(:,:,1) 

PQ = paddedsize(size(blurd));
H = fftshift(fft2(psf_kernel, PQ(1), PQ(2)));
%RO = input('Ro = ');
%RO = typecast(RO,'double');
%while(RO>0)
    GR = fftshift(fft2(double(blurd(:,:,1)), PQ(1), PQ(2)));
    GG = fftshift(fft2(double(blurd(:,:,2)), PQ(1), PQ(2)));
    GB = fftshift(fft2(double(blurd(:,:,3)), PQ(1), PQ(2)));
    F = zeros(PQ(1),PQ(2),3);
    F(:,:,1) = GR;
    F(:,:,2) = GG;
    F(:,:,3) = GB;
    for R = 1:PQ(1)
        for C=1:PQ(2)
%            Rdash = R - PQ(1)/2;
%            Cdash = C - PQ(2)/2;
%            if (Rdash^2 + Cdash^2)<(RO^2)
            %if PQ(1)<0
                F(R,C,1) = GR(R,C)/H(R,C);
                F(R,C,2) = GG(R,C)/H(R,C);
                F(R,C,3) = GB(R,C)/H(R,C);
%            end
        end
    end

    filtered_image = butterworthbpf(F,0,10,10);
    R_filtered = ifft2(fftshift(F(:,:,1)));
    R_filtered = R_filtered(2:size(blurd,1)+1, 2:size(blurd,2)+1);
    G_filtered = ifft2(fftshift(F(:,:,2)));
    G_filtered = G_filtered(2:size(blurd,1)+1, 2:size(blurd,2)+1);
    B_filtered = ifft2(fftshift(F(:,:,3)));
    B_filtered = B_filtered(2:size(blurd,1)+1, 2:size(blurd,2)+1);
%Display results (show all values)

    RGB_filtered = zeros(size(R_filtered,1),size(R_filtered,2),3,'uint8');

    RGB_filtered(:,:,1) = real(R_filtered);
    RGB_filtered(:,:,2) = real(G_filtered);
    RGB_filtered(:,:,3) = real(B_filtered);

    figure, imshow(RGB_filtered,[])
%    RO = input('Ro = ');
%    peaksnr = psnr(RGB_filtered,'GroundTruth1_1_1.jpg');
%    disp(peaksnr)
%end
function F=dft2(f,M1,N1)
M=size(f,1);
N=size(f,2);

%zero padding
if M< M1
    pad_size= (M1-M);
    pad= zeros(round(pad_size/2),N);
    f=[pad; f;pad];
    M= M1;
end


if N< N1
    pad_size= (N1-N);
    pad= zeros(M,round(pad_size/2));
    f=[pad f pad];
    N=N1;
end

H=zeros(M,N);
% shift DFT to centre
for k=1:1:M
    for m=1:1:N
      H(k,m)=((-1)^(k+m));
    end
end

f=f.*H;
Fx= dft(f);
F= transpose(dft(transpose(Fx)));

% DFT corresponds to image with origin at centre of image  
F=F.*H;
      
end
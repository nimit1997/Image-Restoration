function z=idft(x)

N=size(x,2);
ee=exp(2*pi*1i/N);
exp_matrix=zeros(N,N);
for i=1:1:N
    for j=1:1:N
        exp_matrix(i,j)=ee^((j-1)*(i-1));
    end
end
z= 1/N*x*exp_matrix;
end
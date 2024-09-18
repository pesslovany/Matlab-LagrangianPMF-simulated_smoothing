function [pdfIn, kernel] = convFftN(pdfIn, kernel, Npa, nx, kernelFFT)
% Calculates convolution by FFT
%INPUTS:
% pdfIn - pdf for convolution
% kernel - kernel for convolution
% Npa - number of points per axis
% nx - dimension
% kernelFFT - is kernel already in frequency space?
%OUTPUTS:
% pdfIn - convolution result

dims = 1:1:nx;
% Will be used to truncate the padding need the do the convolution
ifun = @(m,n) ceil((n-1)/2)+(1:m);
subs(1:ndims(pdfIn)) = {':'};

for dim=dims % FFT over all dimensions
    % compute the FFT length with padding
    l = Npa(dim)+Npa(dim)-1;
    pdfIn = fft(pdfIn,l,dim); % FFT of the PDF
    if ~kernelFFT
        kernel = fft(kernel,l,dim); % FFT of the kernel
    end
    subs{dim} = ifun(Npa(dim),Npa(dim)); % Padding indices
end

% Perform convolution
pdfIn = pdfIn.*kernel; % without mex file
%inplaceprod(pdfIn, kernel); % with mex file

% Back to state space
for dim=dims
    pdfIn = ifft(pdfIn,[],dim);
end

% Make sure the result is real
pdfIn = real(pdfIn(subs{:}));

end


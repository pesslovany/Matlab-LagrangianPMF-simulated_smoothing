function [predPdf,predGrid,gridStep,convKerTens] = gbfTimeUpdateFFT(F,measPdf,measGrid,gridStep,k,Npa,invQ,predDenDenomW,nx,u)
%PMFUPDATEFFT time update in the Lagrangian formulation
%INPUTS:
% F - model dynamics matrix
% measPdf - measurement pdf
% measGrid - measurement grid
% gridStep - measurement grid step per dimension (and all old grid steps)
% k - time step
% Npa - number of points per axis
% invQ - inversion of dynamics noise
% predDenDenomW - constant of the dynamics gaussian noise function 
% nx - dimension of the state
% u - input to the model dynamics
%OUTPUTS:
% predPdf - predictive pdf
% predGrid - predictive grid
% gridStep - all grids steps with added predictive grid step
% convKerTens - convolution kernel in frequency space

% Pred Grid
predGrid = F*measGrid + u; % Predictive grid
gridStep(:,k+1) = F*gridStep(:,k); % Predictive grid step size

measPdfDotDeltas = (measPdf*prod(gridStep(:,k))); % measurement PDF * measurement PDF step size
mesPdfDotDeltasTens = reshape(measPdfDotDeltas,Npa); % Into tensor space

halfGridInd = ceil(length(predGrid)/2); % Index of middle point of predictive grid
distPoints = (predGrid(:,halfGridInd)'-(predGrid)'); % Distance of transformed grid points to the new grid points
convKer = ((exp(sum(-0.5*distPoints*invQ.*distPoints,2)))/predDenDenomW)';% Convolution kernel values
convKerTens = reshape(convKer,Npa); % Convolution kernel values in tensor format

[predDensityProb2cub, convKerTens] = convFftN(mesPdfDotDeltasTens, convKerTens, Npa, nx, 0); % FFT convolution

predPdf = reshape(predDensityProb2cub,length(predGrid),1); % back to vector for easier manipulation
predPdf = predPdf./(sum(predPdf)*prod(gridStep(:,k+1)))'; % Normalizaton (theoretically not needed)
predPdf(predPdf<0) = 0; % FFT approach sometimes can produce very small negative values

end


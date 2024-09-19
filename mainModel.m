%%
%   Lagrangian based grid-based filter with smoothing for simulated models.
%   authors: Jakub Matousek, pesslovany@gmail.com
%            Jindrich Dunik, dunikj@kky.zcu.cz, University of West Bohemia
% See paper: TBD doi.

%% Parameters and system simulation
clc
clear variables
close all
format shortG

modelChoose = 1; % Which model to initialize
smoothing = 1; % Do you want to perform smoothing

model = initModel(modelChoose); % Initialize model

% Unpack model structure variables
fields = fieldnames(model);
for i = 1:numel(fields)
    eval([fields{i} ' = model.' fields{i} ';']);
end

clear model

% PMF initialization and parameters

% Initial grid
[predGrid, predGridDelta, gridDimOld, gridCenter, gridRotation] = gridCreation(meanX0,varX0,sFactor,nx,Npa);

% Initial PMD
fixTerm = (predGrid-meanX0);
denominator = sqrt((2*pi)^nx*det(varX0));
predPdf = ((exp(sum(-0.5*fixTerm'/(varX0).*fixTerm',2)))/denominator); % Initial Gaussian Point mass density (PMD)

% Auxiliary variables
predDenDenomW = sqrt((2*pi)^nx*det(Q)); % Denominator for convolution in predictive step
halfGrid = ceil(N/2); % Middle row of the TPM matrix index

% For cycle over time stes
for k = 1:1:kf

    % Grid-based Filter
    tic
    
    % Measurement update
    [filtPdf] = gbfMeasUpdate(predGrid,nz,k,z(:,k),V,predPdf,predGridDelta(:,k),hfunct); % Measurement update
    
    % Filtering mean and var
    filtMeanPMF(:,k) = predGrid*filtPdf*prod(predGridDelta(:,k)); %#ok<*SAGROW> % Measurement update mean
    chip_ = (predGrid-filtMeanPMF(:,k));
    chip_w = chip_.*repmat(filtPdf',nx,1);
    filtVarPMF(:,:,k) = chip_w*chip_' * prod(predGridDelta(:,k)); % Measurement update variance
    filtGrid = predGrid;
    
    % Meas PDF interp
    [filtPdf,gridDimOld, GridDelta(:,k), measGridNew, eigVect] = measPdfInterp(filtPdf,...
        gridDimOld,gridCenter,gridRotation,Q,sFactor,nx,Npa,filtMeanPMF(:,k),filtVarPMF(:,:,k),F);
    gridCenter = F*filtMeanPMF(:,k) + u(:,k);
    gridRotation = F*(eigVect);
    
    % Save variables for smoothing
    if smoothing
        measPdfOut(:,k) = filtPdf;
        filtGridOut{k} = measGridNew;
        filtGridDelta{k} = GridDelta(:,k);
        filtGridDim{k} = gridDimOld;
        gridRotOut{k} = eigVect;
        gridCentOut{k} = filtMeanPMF(:,k) + u(:,k);
    end
    
    % Time Update
    [predPdf,predGrid,predGridDelta,TPMrowCubPom{k+1}] = gbfTimeUpdateFFT(F,filtPdf,measGridNew,GridDelta,k,Npa,invQ,predDenDenomW,nx,u(:,k));

    % Save variables for smoothing
    if smoothing
        predDensityProbOut(:,k+1) = predPdf;
    end
    
    tocPMF(k) = toc; % Time evaluation

end

%% Rauch–Tung–Striebel smoothing for Lagrangian Grid based Filters

if smoothing

    % First smoothing estimate is filtering 
    smoothMean(:,kf) = filtMeanPMF(:,kf);
    smoothVar(:,:,kf) = filtVarPMF(:,:,kf);

    % First smoothed PDF is filtering
    smoothed{kf} = reshape(measPdfOut(:,kf),Npa);

    % Go backward through time
    for k = kf-1:-1:1

        % Last smoothed PDF (in k+1) has to be interpolated on the prediction k+1
        % grid so that they can be divided as is in the smoothing equation

        % For interpolation the whole space is rotated and moved so that
        % the grid from which we are interpolating has boundaries aligned
        % with state space axes and its center is in the origin of the
        % coordinate system. Then efficient griddedInterolation can be
        % used. More in the accompanying paper.
        
        % Grid to interp on \xi_{k+1}
        gridInterp = F*filtGridOut{k};
        
        % Interpolant - \xi_{k+1}^interp (rotated and moved to origin)
        Fint = griddedInterpolant(filtGridDim{k+1},reshape(smoothed{k+1},Npa),"linear","none");

        % Grid to be interpolated rotated and moved to the interpolation
        % space
        inerpOn = inv(gridRotOut{k+1})*(gridInterp - gridCentOut{k+1});

        % Interpolation
        interpSmooth = Fint(inerpOn')'; 
        interpSmooth(isnan(interpSmooth)) = realmin; % "Zeros" for extrapolation, otherwise artifacts would appear

        % Last smooth divided by predictive
        smoothDivPred = reshape(interpSmooth,Npa)./reshape(predDensityProbOut(:,k+1),Npa);

        % Put smoothDivPred on filtering grid (i.e. calculate the integral/convolution in the smoothing equation by FFT)
        [fixTerm] = convFftN(smoothDivPred, TPMrowCubPom{k+1}, Npa, nx, 1);

        % Correct the filtering estimate to get the smoothing estimate
        smoothed{k} = abs(fixTerm.*reshape(measPdfOut(:,k),Npa));
        smoothed{k} = smoothed{k}/(sum(smoothed{k},"all")*prod(filtGridDelta{k}));

        % Calculate mean and covariance of smoothed estimate
        smoothMean(:,k) = (filtGridOut{k})*reshape(smoothed{k},N,1)*prod(filtGridDelta{k}); % Measurement update mean
        chip_ = (filtGridOut{k}-smoothMean(:,k));
        chip_w = chip_.*repmat(reshape(smoothed{k},1,N),nx,1);
        smoothVar(:,:,k) = chip_w*chip_' * prod(filtGridDelta{k}); % Measurement update variance

    end

end

rmsePMF = sqrt(mean((x(:,1:kf-1)-filtMeanPMF(:,1:kf-1)).^2,2)); %#ok<*SAGROW>
rmsePMFsmooth = sqrt(mean((x(:,1:kf-1)-smoothMean(:,1:kf-1)).^2,2)); %#ok<*SAGROW>

% Smooth
plot(smoothMean(1,:),smoothMean(2,:),"LineWidth",2,"Color","green")
hold on
plot(filtMeanPMF(1,:),filtMeanPMF(2,:),"LineWidth",2,"Color","red")
plot(x(1,:),x(2,:),Color="red",Marker=".")
legend('smooth','filt','true state')

 

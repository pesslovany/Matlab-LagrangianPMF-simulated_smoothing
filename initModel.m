function [model] = initModel(modelChoose)
%initModel initialize variables for estimation that are model dependent

switch modelChoose
    case 1 % Simple example 2D model
        % System parameters
        model.kf = 30;
        model.nx = 2; % state dimension
        model.nz = 2; % measurement dimension
        model.dt = 1; % time step
        model.Q = 10*eye(2); % system noise
        model.u = zeros(model.nx,model.kf); % Input
        model.invQ = inv(model.Q);
        model.R = eye(2); % measurement noise covariance for both modes
        model.invR = inv(model.R);
        % PMF parameters
        model.Npa = [201 201]; % number of points per axis
        model.N = prod(model.Npa); % number of points - total
        model.sFactor = 4; % scaling factor (number of sigmas covered by the grid)
        model.meanV = [0;0]; % Mean values of components of meas noise
        model.wV = 1; % weights
        % Initial condition - Gaussian
        model.meanX0 = [100; 100];% initial cond
        model.varX0 = [100 0;
            0 100]; % initial cond variance
        theta = deg2rad(5);
        model.F = [cos(theta), -sin(theta); 
         sin(theta),  cos(theta)]; % model dynamics
        model.x = mvnrnd(model.meanX0,model.varX0,1)';
        w = mvnrnd(zeros(model.nx,1),model.Q, model.kf)';
        model.hfunct = @(x,v,k) [0.1 0.9; 0.8 0.3]*x + v; % measurement equation
        for k = 1:model.kf-1
            model.x(:,k+1) = model.F*model.x(:,k) + w(:,k);
        end
        % measurement generation
        model.z = model.hfunct(model.x,0,0)+sqrt(model.R)*randn(model.nz,length(model.x));
        model.V.pdf = gmdistribution(model.meanV',repmat(model.R,1,1,1),model.wV); % Meas noise pdf
end

end
function [speedFitted,tStartFitted] = findSpeed(Bridge,t0,Doz,posAcc,tImpact_guess,varargin)
% [speedFitted,tStartFitted] =
% findSpeed(Bridge,t0,Doz,posAcc,tImpact,varargin) estimate the arrival
% time and the speed of a vehicle crossing a bridge.
%
% Input
%
%   - Bridge: structure Bridge including its parameter (cf. example1.mlx)
%   - t0: [1xN] double: time vector (in seconds)
%   - Doz: [1xN] double: Vertical response of the bridge at one location, ideally at midspan.
%   - posAcc: [1x1]: relative position of the accelerometer along the span
%       (between 0 and 1)
%   - tImpact_guess: Intiial guess for the arrival time of the vehicle
%   (obtained from the clustering analysis)
%   - varargin:
%       - orderFilt: order of the high-pass and low-pass filters
%       - fcUp: [1x1] double: cutoff frequency of the low-pass filter (in Hz)
%       - fcLow: [1x1] double: cutoff frequency of the high-pass filter (in Hz)
%       - newN: [1x1] integer: Number of datapoint for the simulated bridge response (lower means faster but also less accurate)
%       - meanU: [1x1] double:Mean wind speed, required if the aerodynamic damping is accounted for.
%       - plotData: Display in a graph the fitted vs measured response. It is useful for validation purpose
%       - gradS: [1x1] double: termination tolerance for the speed estimate, i.e. the speed is identified at +/- gradS (in m/s)
%       - gradT: [1x1] double: termination tolerance for the arrival time estimate, i.e. the speed is identified at +/- gradS (in m/s)
%       - RMSE_threshold: [1x1] double: RMSE value above which the vehicle
%       direction is suspected to be wrongly accounted for, leading to a
%       new fitting but with an opposite direction.
%       -  deltaT: Data located at "tImpact_guess +/-deltaT" are including
%       for the fitting -> This value depends on the bridge
%
% Output
%   - speedFitted: [1xM]: estimated speed of the M vehicles identified
%   - tStartFitted: [1xM]: estimated arrial time of the M vehicles identified
%
%
% Author info: E. Cheynet  - UiS/UiB - last modified: 02-08-2020
% see also dynaResp_vehicle_TD findSpeed getVehicleID

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('orderFilt',4); % order of the high-pass and low-pass filter
p.addOptional('fcUp',0.15); % cut-off frequency (low pass filter) -> must be lower than the first eigen frequency
p.addOptional('fcLow',0.04); % cut-off frequency (high pass filter) -> must be above the lower frequency correctly measured by the accelerometer
p.addOptional('plotData',1);
p.addOptional('RMSE_threshold',0.3);
p.addOptional('gradS',0.5); % absolute change in speed (m/s)
p.addOptional('gradT',0.5); % absolute change in time (s)
p.addOptional('newN',150);
p.addOptional('deltaT',30);
p.addOptional('meanU',0);
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
newN = p.Results.newN ;
orderFilt = p.Results.orderFilt ;
fcUp = p.Results.fcUp;
fcLow = p.Results.fcLow;
gradS = p.Results.gradS;
gradT = p.Results.gradT;
meanU = p.Results.meanU;
deltaT = p.Results.deltaT;
plotData = p.Results.plotData;
RMSE_threshold = p.Results.RMSE_threshold;

%% Initalization of the data and parameters definition
Nvehicle = numel(tImpact_guess);
speedFitted = zeros(1,Nvehicle);
tStartFitted = zeros(1,Nvehicle);

% Guess values of the vehicle
vehicle.mass = 2000;
vehicle.direction = 1;
vehicle.wn = 0;

% For the 2D fitting, the surface response is calculated on a Nt x Ns
% matrix
Nt = 8; % number of step for tImpact
Ns = 8; % number of step for speed

% Maximal and minimal speed accepted for the fitting
minSpeed = 25; % km/h
maxSpeed = 80; % km/h


%% Check for inconsitent parameters
[Nsensors,N]= size(Doz);

if Nsensors>1 && N==1
    Doz = Doz';
    [Nsensors,N]= size(Doz);
else
    error('This version of the code does not include the possibility to include multiple accelerometers');
end

if numel(posAcc) ~=Nsensors, error('numel(posAcc) ~= size(DOz,1)'); end
if max(abs(posAcc))>1, error('"posAcc" should be the relative position of the accelerometer, not its absolte position');end

[~,indY] = min(abs(posAcc-Bridge.x));


%% Least-square fitting algorithm

for ii=1:Nvehicle
    
    t1 = tImpact_guess(ii)-deltaT; % get lower boundary for arrival time
    t2 = tImpact_guess(ii)+deltaT; % get upper boundary for arrival time
    
    newT = linspace(t1,t2,newN); % new time vector with reduced number of datapoint
    newFs = 1/median(diff(newT)); % Get the new sampling frequency
    Wind.u = meanU + zeros(Bridge.Nyy,newN); % Include mean wind speed if needed
    Wind.w = zeros(Bridge.Nyy,newN);
    Wind.t = newT;
    
    % Interpolate the bridge response over the new time vector
    rz = interp1(t0,Doz,newT);
    
    % Isolate the background response
    rz = filterMyData(rz,newFs,orderFilt,fcUp,fcLow);
    
    % Normalzie the response
    rz = rz./max(abs(rz));
    
    % Check if there exist NaNs
    if any(isnan(rz)), warning('There exist nans values in the rz time series. Consider replacing them by interpolated values'); end
    
    
    tic
    % We assume one direction first (the one sepcified by the vehicle)
    [vehicle,rSim,newRMSE_0,X1_0,X2_0,speedFitted(ii),tStartFitted(ii)] = fitData(gradS,gradT,vehicle,rz,newFs,indY,t1,tImpact_guess(ii));
    
    myRMSE = RMSE(rSim(:),rz(:));
    
    % If the RMSE is larger than the acceptable threshold value, this may be
    % due to the wrong assessement of the vehicle direction.
    if myRMSE>RMSE_threshold
        fprintf(['RMSE value is ',num2str(myRMSE,2),'. New attempt is made by changing the vehicle direction \n']);
        vehicle.direction = -1*vehicle.direction; % Change the direction of the vehicle
        [vehicle,rSim,newRMSE_0,X1_0,X2_0,speedFitted(ii),tStartFitted(ii)] = fitData(gradS,gradT,vehicle,rz,newFs,indY,t1,tImpact_guess(ii));
        myRMSE = RMSE(rSim(:),rz(:));
        fprintf(['new RMSE value is ',num2str(myRMSE,2),'. \n']);
    end
    toc
    
    if plotData==1
        if Nvehicle<=4
            subplot(2,2,ii)
        elseif Nvehicle<10 &&  Nvehicle>4
            subplot(3,3,ii)
        else
            subplot(Nvehicle,1,ii)
        end
        imagesc(X1_0(1,:),X2_0(:,1),newRMSE_0);
        axis xy
        colorbar; shading flat; hold on;
        h = plot(speedFitted(ii),tStartFitted(ii),'r*');
        colormap pink
        hh = colorbar;
        ylabel(hh, 'RMSE (a.u.)')
        caxis([min(newRMSE_0(:)) max(newRMSE_0(:))])
        xlabel(' Vehicle speed (m/s)')
        ylabel(' Starting time (s)')
        set(gcf,'color','w')
        legend(h,'Minimum','location','best');
        axis tight
    end
end

%% Nested functions
    function [vehicle,rSim,rmse0,X1,X2,speedFitted,TimpactFitted] = fitData(gradS,gradT,vehicle,rz,newFs,indY,t1,tImpact)
% [vehicle,rSim,newRMSE_0,X1_0,X2_0,speedFitted,tStartFitted] =
% fitData(gradS,gradT,vehicle,rz,newFs,indY,t1,tImpact) estimate the
% vehicle speed and arrival time by fitting the bridge response to a moving mass model to the
% recorded filtered displacement records
% 
% Inputs:
% 
%   - gradS: [1x1] double: termination tolerance for the speed estimate, i.e. the speed is identified at +/- gradS (in m/s)
%   - gradT: [1x1] double: termination tolerance for the arrival time estimate, i.e. the speed is identified at +/- gradS (in m/s)
%   - vehicle: structure including the vehicle parameters
%   - rz: [1xnewN] Bridge displacement response for a given cluster
%   - newFs: sampling frequency
%   - indY: relative position along the span of the accelerometer
%   - tImpact: guess of the arrival time.
% 
% Outputs:
% 
%   - vehicle: structure including the updated vehicle parameters
%   - rSim: COmputed bridge displacement response for a given cluster
%   - newRMSE_0: RMSE value between the fitted and recorded filtered
%   displacement response
%   - X1_0: range of  vehicle's speed for plotting purpose only
%   - X2_0: range of vehicles' arrival time for plotting purpose only
%   - speedFitted: estimated vehicle speed
%   - TimpactFitted: estimated arrival time
% 
% Author: E. Cheynet - UiB - last modified: 03-08-2020s
% 
        dS = 1; % Initial value of the gradient for the speed
        dT = 1; % Initial value of the gradient for the impact time
        count = 1;
        
        ArrivalTarget = linspace(t1,tImpact,Nt); % vector of possible arrival time
        SpeedTarget = linspace(minSpeed,maxSpeed,Ns)./3.6;
        speed0 = nanmean(SpeedTarget); % Initial value for target speed (guess)
        tStart0 = nanmean(ArrivalTarget);  % Initial value for the impact time  (guess)
        
        while dS>=gradS &&  dT>=gradT
            % while the termination tolerance is not reached, keep fitting
            if count==1 % First iteration
                [rmse0,speed1,tStart1,X1,X2]= getRMSE_data(SpeedTarget,ArrivalTarget,newFs,rz,indY,vehicle.direction);
            else
                [~,speed1,tStart1,~,~]= getRMSE_data(SpeedTarget,ArrivalTarget,newFs,rz,indY,vehicle.direction); 
            end
            dS  = abs(speed1-speed0); % absolute difference
            dT = abs(tStart1-tStart0); % absolute difference
            
            speed0 = speed1; % old speed is updated
            tStart0 = tStart1; % old time is updated
            minT1 = max([ArrivalTarget(1), tStart1-(mean(ArrivalTarget)-ArrivalTarget(1))/1.5]); % recalculate the lower boundary for the arrival time
            maxT1 = min([ArrivalTarget(end), tStart1+(mean(ArrivalTarget)-ArrivalTarget(1))/1.5]); % recalculate the upper boundary for the arrival time
            minS1 =  max([SpeedTarget(1), speed1-(mean(SpeedTarget)-SpeedTarget(1))/1.5]); % recalculate the lower boundary for the vehicle speed
            maxS1 = min([SpeedTarget(end), speed1 + (mean(SpeedTarget)-SpeedTarget(1))/1.5]); % recalculate the upper boundary for the vehicle speed
            ArrivalTarget = linspace(minT1,maxT1,Nt);
            SpeedTarget = linspace(minS1,maxS1,Ns);
            count = count+1;
            if count>=15
                warning('Fitting failed');
                speedFitted = zeros(1,Nvehicle);
                TimpactFitted = zeros(1,Nvehicle);
                return
            end
        end
        speedFitted = speed1;
        TimpactFitted = tStart1;
        vehicle.speed = speedFitted;
        vehicle.tStart = TimpactFitted;
        
        % Compute the bridge response with the fitted vehicle parameters
        [Do1,~] = dynaResp_vehicle_TD(Bridge,Wind,'vehicle',vehicle,'comments','off'); % both wind and traffic
        [~,indY0] = min(abs(posAcc-Bridge.x));
        rSim = squeeze(Do1(2,indY0,:)); % indice 2 for the vertical motion
        rSim = filterMyData(rSim(:)',newFs,orderFilt,fcUp,fcLow); % isolate the background response
        rSim = rSim./max(abs(rSim)); % normalize the response
    end
    function [valRMSE,speedFitted,TimpactFitted,X11,X22]= getRMSE_data(Speed,Timpact,fs,rz,indY,direction)
        % [newRMSE,speedFitted,TimpactFitted,X11,X22]=
        % getRMSE_data(SpeedTarget,Tsync,fs,rz,indY,direction) computes the RMSE
        % between the fitted and target data.
        %
        % Inputs:
        %   - Speed: [1xNs] target vehicles' speed used for the surface response
        %   computation
        %   - Timpact: [1xNt]: target vehicles' arrival time used for the surface
        %   response computation
        %   - fs: [1x1] sampling frequency
        %   - rz: [1xN]: Bridge response at one location
        %   - indY: relative position of the accelerometer
        %   - direction: Direction of the vehicle (1 or -1)
        %
        % Outputs:
        %
        % - valRMSE: [1x1]  RMSE value
        % - speedFitted: estimated vehicles' speed associated with the lwoest RMSE value [1xM]
        % - TimpactFitted: estimated vehicles' arrival time  associated with the lwoest RMSE value [1xM]
        % - X11: range of  vehicle's speed for plotting purpose only
        % - X22: range of vehicles' arrival time for plotting purpose only
        %
        % Author: E. Cheynet - UiB - last modified: 03-08-2020
        
        surfaceResponse = nan(Ns,Nt);
        for i1=1:Ns
            for i2=1:Nt
                vehicle.speed = Speed(i1);
                vehicle.tStart = Timpact(i2);
                vehicle.direction = direction;
                
                [dummyRz,~] = dynaResp_vehicle_TD(Bridge,Wind,'vehicle',vehicle,'comments','off'); % displacement data
                dummyRz = squeeze(dummyRz(2,indY,:))'; % Isolate the vertical response (the lateral response is the first dimension and the torsional response is the third one)
                dummyRz = filterMyData(dummyRz(:)',fs,orderFilt,fcUp,fcLow); % Isolate the background response
                dummyRz = dummyRz./max(abs(dummyRz)); % Normalize the response
                [surfaceResponse(i1,i2)] = RMSE(dummyRz(:),rz(:)); % Compute the surface response in terms of RMSE value
            end
        end
        
        % Oversample the surface response for visualization purpose only
        [X11,X22] = meshgrid(linspace(Speed(1),Speed(end),200),linspace(Timpact(1),Timpact(end),200));
        [x1,x2]= meshgrid(Speed,Timpact);
        valRMSE = interp2(x1,x2,surfaceResponse',X11,X22,'linear');
        
        % Get the lowest RMSE value and associate it to the fitted vehicle
        % speed and arrival time.
        speedFitted = X11(valRMSE==min(min(valRMSE)));
        TimpactFitted =X22(valRMSE==min(min(valRMSE)));
    end
end

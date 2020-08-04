function [rz0,t0,tImpact,tCluster,rzCluster] = findVehicleID(rz,t,varargin)
% [rz0,t0,tImpact,tCluster,rzCluster,flag] = findVehicleID(rz,t,varargin)
% identify the transient bridge response events associated with the passage
% of vehicles.  It relies on an outlier detection algorithm and a
% clustering analysis.
% 
% Input
%   - rz: [1xN] double: Vertical response of the bridge at one location, ideally at
%   midspan.
%   - t: [1xN] double: time vector (in seconds)
%   - varargin:
%       - orderFilt: order of the high-pass and low-pass filters
%       - fcUp: [1x1] double: cutoff frequency of the low-pass filter (in Hz)
%       - fcLow: [1x1] double: cutoff frequency of the high-pass filter (in Hz)
%       - coefDecimate:  [1x1] integer: used for decimation of the time series to speed-up
%       the identification of the vehicles
%       - CuttOff:  cutoff value for the clusters, constructed using a hierarchical cluster tree (cf Matlab function "cluster")
%       - 'tMin':[1x1] double:  clusters located at t<tMin are ignored
%       - 'tMax':[1x1] double:  clusters located at t>tmax are ignored
%       - 'minNp'[1x1] integer: clusters with less than 'minNp' points are
%       ignored
%       -'deltaTmax': [1x1] double: clusters located less than 'deltaTmax'
%       seconds are ignored (they can lead to mixed response, not easily
%       separated)
%       - 'Outlier': string: method for the outlier algorithm ('gesd' by
%       default)
% 
% Output
%   - rz0: background response isolated from the total response (used for
%   plotting purpose only)
%   - t0: time vector associated with t0 (numel(t0)<numel(t))
%   - tImpact: [1xM] double: first guess of the arrival time, calculated for
%   each cluster, where M is the number of identified clusters.
%   - tCluster: cell [1xM] where each element is the time vector of one
%       cluster.
%   - rzCluster: cell [1xM] where each element is the bridge response
%     associated with one cluster.
% 
% Author info: E. Cheynet  - UiB - last modified: 02-08-2020
% see also dynaResp_vehicle_TD findSpeed findMass



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('orderFilt',4); % order of the high-pass and low-pass filter
p.addOptional('fcUp',0.15); % cut-off frequency (low pass filter) -> must be lower than the first eigen frequency
p.addOptional('fcLow',0.04); % cut-off frequency (high pass filter) -> must be above the lower frequency correctly measured by the accelerometer
p.addOptional('coefDecimate',round(1/(5*median(diff(t)))));
p.addOptional('CuttOff',300); % distance for the cluster analysis (number of points)
p.addOptional('tMin',50); % starting time
p.addOptional('tMax',20); % ending time
p.addOptional('minNp',30); % number of cluster points
p.addOptional('deltaTmax',80); % maximal time between 2 vehicles (sec)
p.addOptional('Outlier','gesd'); % algorithm for outlier detection
p.parse(varargin{:});
orderFilt = p.Results.orderFilt ;
fcUp = p.Results.fcUp;
fcLow = p.Results.fcLow;
coefDecimate = p.Results.coefDecimate;
CuttOff = p.Results.CuttOff;
tMin = p.Results.tMin;
tMax = p.Results.tMax;
minNp = p.Results.minNp;
deltaTmax = p.Results.deltaTmax;
Outlier = p.Results.Outlier;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Filtering and decimation of the data

fs = 1/median(diff(t)); % get sampling frequency
rz0 = filterMyData(rz(:)',fs,orderFilt,fcUp,fcLow); % Extract the background response
rz0 = decimate(rz0,coefDecimate,coefDecimate,'fir'); % decimation


%% Outlier analysis

t0 = decimate(t,coefDecimate,coefDecimate,'fir');
if strcmpi(Outlier,'gesd')
    [tf] = isoutlier(rz0,'gesd','MaxNumOutliers', round(0.40*numel(rz0)));
elseif strcmpi(Outlier,'median')
    [tf] = isoutlier(rz0,'median','ThresholdFactor',3);
else
    error('Outlier detection algorithm is unknown')
end

tOut = t0(tf==1); % Time vector associated with outliers
rzOut = rz0(tf==1); % bridge response associated with outliers

%% Clustering analysis
if numel(tOut)<=1 % if no outlier detected
    rz0 = [];
    t0 = [];
    tImpact = [];
    tCluster = [];
    rzCluster = [];
    return
else
    
    % Create hierarchical cluster tree.
    Z = linkage(squareform(pdist(tOut')),'single','euclidean');
    myClus = cluster(Z,'Cutoff',CuttOff,'Criterion','distance');
    uniqueCluster = unique(myClus);
    Ncluster = numel(uniqueCluster);
    
    clear tImpact rzCluster tCluster
    ll=1;
    
    for jj=1:Ncluster
        dummyT = median(tOut(myClus==uniqueCluster(jj))); % get median time for each cluster
        Ndata = numel(myClus(myClus==uniqueCluster(jj))); % get number of point for each cluster
        % if traffic a the beginning or end of the time series, dismiss it
        % If the cluster has not enough point, dismiss it
        if dummyT<t(1)+tMin || dummyT>t(end)-tMax || Ndata<minNp || isnan(uniqueCluster(jj))
            myClus(myClus==uniqueCluster(jj)) = nan;
        else
            tImpact(ll) = nanmedian(tOut(myClus==uniqueCluster(jj)));
            tCluster{ll} =  tOut(myClus==uniqueCluster(jj));
            rzCluster{ll} = rzOut(myClus==uniqueCluster(jj));
            ll=ll+1;
        end
    end
    
    
    % Check if two adjacent clusters are not two close
    if exist('tCluster','var') && exist('rzCluster','var')
        [tCluster,rzCluster,tImpact] = cleanCluster(tCluster,rzCluster,tImpact,deltaTmax);
    else
        tCluster = [];
        rzCluster = [];
        tImpact = [];
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [tCluster,rzCluster,tImpact] = cleanCluster(t0,rz0,tImpact0,deltaTmax)
        
        [tImpact,indSort]  = sort(tImpact0); % sort the clusters
        tCluster = t0(indSort);
        rzCluster = rz0(indSort);
        
        % Check if two cluster are not too close
        indExclu = [];
        qq= 1;
        for ii=1:numel(tCluster)
            deltaT = max(tCluster{ii})-min(tCluster{ii});
            if deltaT > deltaTmax,
                warning(['deltaT = ',num2str(deltaT,3),' sec. It looks like that two vehicle are too close to each other to be identified properly'])
                indExclu(qq)=ii;
                qq=qq+1;
            end
        end
        
        tCluster(indExclu) = [];
        rzCluster(indExclu) = [];
        tImpact(indExclu) = [];
        
        tCluster = tCluster(~cellfun('isempty',tCluster)) ;
        rzCluster = rzCluster(~cellfun('isempty',rzCluster)) ;
    end
end
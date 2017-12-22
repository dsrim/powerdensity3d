function threshold_plot(sigmatype,N,k,cutoff,varargin)

if ~isempty(varargin)
if strcmp(varargin{1}, 'QT')
    display(['(' mfilename ') loading QreconR ...'])
    figure; clf;
    load(['QreconR_' sigmatype '_' num2str(N) '.mat'])
    bderror4R = error4R.*(error4R < cutoff) + cutoff.*(error4R >= cutoff); 
    ezplot3D(bderror4R,3)
    title('error')
elseif strcmp(varargin{1},'sigma')
    display(['(' mfilename ') loading reconSigma ...'])
    load(['QreconSigma_' sigmatype '_' num2str(N) '.mat'])
    figure; clf;
    bdsigma = sigma.*(sigma < cutoff) + cutoff.*(sigma >= cutoff);
    ezplot3D(bdsigma,3)
    title({'Sigma'; ['threshold = ' num2str(cutoff)]})
elseif strcmp(varargin{1},'abcd')
    display(['(' mfilename ') loading a,b,c,d ...'])
    load(['abcd_' sigmatype '_' num2str(N) '.mat'])
    f = sqrt(a.^2 + b.^2 + c.^2 + d.^2);
    figure; clf;
    bdf = f.*(f < cutoff) + cutoff.*(f >= cutoff);
    ezplot3D(bdf,k)
    title({'sqrt(a^2 + b^2 + c^2 + d^2)'; ['threshold = ' num2str(cutoff)]})
end
else
    display(['(' mfilename ') loading reconR ...'])
    figure; clf;
    load(['reconR_' sigmatype '_' num2str(N) '.mat'])
    bderror4R = error4R.*(error4R < cutoff) + cutoff.*(error4R >= cutoff); 
    ezplot3D(bderror4R,3)
    title('error')
end

end
% 
% 
% load(['reconRQ_' sigmatype '_' num2str(N) '.mat'])
% bderror4RQ = error4RQ.*(error4RQ < cutoff) + cutoff.*(error4RQ >= cutoff); 
% ezplot3D(bderror4RQ,3)
% title('error from Q')

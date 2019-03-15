function [wave,power,frequencies] = modwtDecomp( data,Fs,wtype,wlevel )
% [wave,power,frequencies] = modwtDecomp( data,Fs,wtype,wlevel )
%
% Computes the maximal-overlap discrete wavelet transform (modwt) over each
% column in the n x m data matrix "data" using the wavelet type specified
% by "wtype", and the wavelet level specified by "wlevel".
%
%                   >>> INPUTS >>>
% data: 
%   n x m data matrix in column-major format
% Fs:
%   the sampling rate
% wtype:
%   string specifying the wavelet type 
% wlevel:
%   the final wavelet level to decompose 
%
%                   <<< OUTPUTS <<<
% wave:
%   n x wlevel x m tensor of wavelet coefficients from the modwt
% power:
%   |wave|^2
% frequencies:
%   wlevel+1 x 2 matrix specifying the frequency bandwidth in each level
% 
% by JMS, 5/20/2016

% maximal-overlap discrete wavelet transform (MODWT)
[n,m] = size( data );
wave = zeros( wlevel+1,n,m );
for j = 1:m
    wave(:,:,j) = modwtmra( modwt( data(:,j),wtype,wlevel ) );
end
power = abs( wave ).^2; % compute wavelet power spectrum
frequencies = wave2freq( wlevel,Fs );

end
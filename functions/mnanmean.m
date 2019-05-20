function y = mnanmean(x)
% MNANMEAN NaN protected median value.
%   MNANMEAN(X) returns the mean treating NaNs as missing values.
%   For vectors, MNANMEAN(X) is the mean value of the non-NaN
%   elements in X.  For matrices, MNANMEAN(X) is a row vector
%   containing the mean value of each column, ignoring NaNs.
%
%   See also NANMEAN, NANSTD, NANMIN, NANMAX, NANSUM.

%   by John Peter Acklam, jacklam@math.uio.no
%   
  
x=x(:);  
  
if isreal(x)  
  ii = ~isnan(x);
  if any(ii)
    y = sum(x(ii))/sum(ii);
  else
    y = NaN;
  end
else
    ii = ~isnan(x);
  if any(ii)
    yr = sum(real(x(ii)))/sum(ii);
    yi = sum(imag(x(ii)))/sum(ii);
    y=yr+i*yi;
  else
    y = NaN;
  end
  
end
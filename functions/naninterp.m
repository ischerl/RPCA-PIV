function [u,v]=naninterp(u,v,varargin)
% function [u,v]=naninterp(u,v,method,mask,x,y)
%
% Interpolates NaN's in a vectorfield. Used by GLOBFILT, MEDIANFILT and
% VALIDATE. Sorts all spurious vectors based on the number of spurous
% neighbors to a point.  Uses FILLMISS.M to replace NaN's
% [u,v]=NANINTERP(u,v) Will replace all NaN's in u and v
% [u,v]=NANINTERP(u,v,[idxw idyw],x,y) Will replace NaN's but leave out
% areas that have been masked out (using MASK.M) Using the MASK option
% requires the x and y matrices input along with the u and v.  METHOD
% should be 'linear' or 'weighted' and defaults to 'linear'.
%  See also: FILLMISS, NANINTERP2, MATPIV

% Copyright 1999 - 2001, J. Kristian Sveen (jks@math.uio.no)
% Dept. of Mathematics, Mechanics Division, University of Oslo, Norway
%
% For use with MatPIV 1.5, 
% Distributed under the terms of the GNU - GPL license
% timestamp: 22:04, 20 Feb 2001

usr=1;
if ~any(strcmp(varargin,'linear')) & ~any(strcmp(varargin,'weighted'))
    met='linear';
else
    tm=cellfun('isclass',varargin,'char');
    if sum(tm)==3
        disp('Wrong input to naninterp!'); return
    end
    met=varargin(tm(1));
end

%%%%%%%%%%%%%%%%%%%%%
if strcmp(met,'weighted')==1 & ~strcmp(met,'linear')==1
    
    if length(varargin)==4
        maske=varargin{2}; if ischar(maske) & ~isempty(maske), maske=load(maske); maske=maske.maske; end  
        if isempty(maske)
            [u,v]=naninterp2(u,v); usr=any(isnan(u(:)));
        else      
            xx=varargin{3}; yy=varargin{4};
            while usr~=0
                in2=zeros(size(xx));
                for i=1:length(maske)
                    in=inpolygon(xx,yy,maske(i).idxw,maske(i).idyw);
                    in2=in2+double(in);                    
                end
                in2=logical(in2);
                % interpolate NaN's using FILLMISS.M
                u(in2)=0; v(in2)=0;
                u=fillmiss(u); v=fillmiss(v);
                usr=any(isnan(u(:)));
                u(in2)=nan; v(in2)=nan;
            end
            u(in2)=NaN; v(in2)=NaN;
            numm=size(in2,1)*size(in2,2);
        end
    else
        while usr~=0
            numm=sum(isnan(u(:)));
            % interpolate NaN's using FILLMISS.M
            u=fillmiss(u); v=fillmiss(v);
            usr=any(isnan(u(:)));
        end
    end
    %fprintf([' -> ',num2str(numm),' NaN''s interpolated\n'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif strcmp(met,'linear')==1 & ~strcmp(met,'weighted')==1
    
    if length(varargin)==4
        maske=varargin{2}; if ischar(maske) & ~isempty(maske),maske=load(maske);maske=maske.maske; end  
        if isempty(maske)
            [u,v]=naninterp2(u,v); usr=any(isnan(u(:)));
        else
            xx=varargin{3}; yy=varargin{4};
            maske=rmfield(maske,'msk'); % this is done to avoid the large matrix 
            %being copied into the next function.
            [u,v]=naninterp2(u,v,maske,xx,yy);
        end
    else
        while usr~=0
            % interpolate NaN's using NANINTERP2.M
            [u,v]=naninterp2(u,v); usr=any(isnan(u(:)));
        end
    end
else
    disp('Something is VERY wrong with your input'); return  
end


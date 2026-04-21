

function[amp_int,distance_int,A_int, TA_int,T_int]=getamp(range_int,R,DIS,A,TA,T ,c0, c_target)
if all(abs(diff(A))<5*eps)
    mid=round((length(A)+1)/2);
    A  = A(mid);
    TA = TA(mid);
end;
if length(A)>=2
    method='linear';
    last1= nan;
    T_int=interp1(R,T,range_int, method,last1);
    A_int=interp1(R,A,range_int, method,last1);
    TA_int=interp1(R,TA,range_int, method,last1);
    distance_int=interp1(R,DIS,range_int, method,last1);
    [amp TL]=ray_loss(R, A*pi/180, TA*pi/180,c0, c_target);
    amp_int=interp1(R(1:length(TL)), amp,range_int, method,0);
elseif length(A)==1
    if abs(A)<10*eps
        T_int    = range_int/c_target;
        A_int    = A*ones(size(range_int));
        TA_int   = TA*ones(size(range_int));
        amp_int  = 1./range_int;
        distance_int=range_int;
    else
        T_int    = nan;
        A_int    = nan;
        TA_int   = nan;
        amp_int  = 0;
        distance_int=nan;
    end;
end;


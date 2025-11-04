% Code modified from Hargreaves

function [k,g,s,time,r,theta,f,pns] = SafeSpiralOut(sys,Tadc,Tgrad,N,Fcoeff,rmax,slewMargin,pnsDesignLimit)

gamma = 4258;       % [Hz/G]

oversamp = 100;		% Keep this even.
To = Tgrad/oversamp;	% To is the period with oversampling.

% Initialise variables
if isfield(sys,'resonFreq')
    resonFreq = sys.resonFreq;
else
    resonFreq = [0 0];
end

fx = zeros(3,1); fy = zeros(3,1);
if isfield(sys,'safeModel')
    if isfield(sys.safeModel,'RIV') && sys.safeModel.RIV == true
        ax = To./(sys.safeModel.tauW + To);
        ay = ax;
        Ax = sys.safeModel.AW;
        Ay = Ax;
        pnsScaling = sys.safeModel.pnsScaling*ones(1,2);
    else
        ax = To./(sys.safeModel.tauX + To);
        ay = To./(sys.safeModel.tauY + To);
        Ax = sys.safeModel.AX;
        Ay = sys.safeModel.AY;
        pnsScaling = sys.safeModel.pnsScaling;
    end
else
    ax = [1, 1, 1];  
    ay = [1, 1, 1];
    Ax = [1, 1, 1];
    Ay = [1, 1, 1];
    pnsScaling = [0, 0];
end


q0 = 0;	q1 = 0;
r0 = 0; r1 = 0;
t = 0;
count = 1;
theta = zeros(1,10000); 
r     = zeros(1,10000); 
time  = zeros(1,10000);
pns   = zeros(2,10000);
f     = zeros(1,10000);
g     = zeros(1,10000);
s     = zeros(1,10000);

gmaxHW     = 0.99*sys.gmax;
smaxHW     = 0.99*sys.smax;
gmaxCur    = gmaxHW;
smaxCur    = slewMargin*smaxHW;

while r0 < rmax
	[q2,r2,slewMin] = findq2r2(smaxCur,gmaxCur,r0,r1,To,Tadc,N,Fcoeff,rmax,smaxHW);

	% Integrate for r, r', theta and theta' 	
	q1 = q1 + q2*To;
	q0 = q0 + q1*To;
	r1 = r1 + r2*To;
	r0 = r0 + r1*To;
    t = t + To;

    % Calculate PNS and control current slew rate. 
    sx = 1/gamma*(r2*cos(q0) - q1*r1*sin(q0) - q2*r0*sin(q0) - q1*r1*sin(q0) - q1*q1*r0*cos(q0) )/1e2; %1e2 from G/cm/s to T/m/s 
    sy = 1/gamma*(r2*sin(q0) + q1*r1*cos(q0) + q2*r0*cos(q0) + q1*r1*cos(q0) - q1*q1*r0*sin(q0) )/1e2; %1e2 from G/cm/s to T/m/s 
    fx(1) = ax(1)*sx     +(1-ax(1))*fx(1);
    fx(2) = ax(2)*abs(sx)+(1-ax(2))*fx(2);
    fx(3) = ax(3)*sx     +(1-ax(3))*fx(3);
    pnsVal(1) = pnsScaling(1)*(Ax(1)*abs(fx(1)) + Ax(2)*fx(2) + Ax(3)*abs(fx(3)));
    fy(1) = ay(1)*sy     +(1-ay(1))*fy(1);
    fy(2) = ay(2)*abs(sy)+(1-ay(2))*fy(2);
    fy(3) = ay(3)*sy     +(1-ay(3))*fy(3);
    pnsVal(2) = pnsScaling(2)*(Ay(1)*abs(fy(1)) + Ay(2)*fy(2) + Ay(3)*abs(fy(3)));
    pnsRms = sqrt(pnsVal(1)^2+pnsVal(2)^2);
    if pnsRms>pnsDesignLimit
        smaxCur = slewMin;
    else 
        smaxCur = slewMargin*smaxHW;
    end
    
    % Calculate current frequency and control current gradient amplitude.
    freq = sqrt(r1^2 + q1^2*r0^2)/2/pi/r0; %[Hz]
    activeReson = false;
    for n=1:size(resonFreq,1)
        if freq<=resonFreq(n,2) && freq>=resonFreq(n,1)
            gmaxReson = resonFreq(n,1)/gamma*2*pi*r0;
            dGmax = findDeltaG(slewMargin*smaxHW,r0,r1,To,N,Fcoeff,rmax);
            gmaxCur = max([dGmax(2), gmaxReson]); % entry #2 corresponds to decreasing amplitude
            activeReson = true;
        end
    end
    if activeReson==false
        gmaxCur = gmaxHW;
    end

	% Store
    if (rem(count,oversamp)==0)
        idx = count/oversamp + 1;
        theta(idx) = q0;
	    r(idx)     = r0;
	    time(idx)  = t;
        g(idx)     = 1/gamma*(r1*exp(1i*q0) + 1i*q1*r0*exp(1i*q0)); 
        s(idx)     = 1/gamma*(r2*exp(1i*q0) + 1i*q1*r1*exp(1i*q0) + 1i*q2*r0*exp(1i*q0) + 1i*q1*r1*exp(1i*q0) - q1^2*r0*exp(1i*q0) );
        f(idx)   = freq;
        pns(:,idx) = pnsVal;
    end
    count = count+1; 
end

theta = theta(1:idx);
r     = r(1:idx);
time  = time(1:idx);
g     = g(1:idx);
s     = s(1:idx);
f     = f(1:idx);
pns   = pns(:,1:idx);

k = r.*exp(1i*theta);

return;



function [q2,r2,smin] = findq2r2(smax,gmax,r,r1,T,Ts,N,Fcoeff,rmax,smaxOrig)

gamma = 4258;			% Hz/G

F = 0;		% FOV function value for this r.
dFdr = 0;		% dFOV/dr for this value of r.
for rind = 1:length(Fcoeff)
	F = F+Fcoeff(rind)*(r/rmax)^(rind-1);
	if (rind>1)
		dFdr = dFdr + (rind-1)*Fcoeff(rind)*(r/rmax)^(rind-2)/rmax;
    end
end

GmaxFOV = 1/gamma /F/Ts;		% FOV limit on G; Ts is dwell time
Gmax = min(GmaxFOV,gmax);	%

twopiFoN = 2*pi*F/N;
twopiFoN2 = twopiFoN^2;

maxr1 = sqrt((gamma*Gmax)^2 / (1+(2*pi*F*r/N)^2));  
if (r1 > maxr1)			
	% Grad amplitude limited.  Here we
	% just run r upward as much as we can without
	% going over the max gradient.
	r2 = (maxr1-r1)/T; 
else
	%	A,B,C are coefficents of the equation which equates
	% 	the slew rate calculated from r,r1,r2 with the
	%	maximum gradient slew rate.
	%
	%	A*r2*r2 + B*r2 + C  =  0	
	%
	%	A,B,C are in terms of F,dF/dr,r,r1, N and smax.
	%
	A = 1+twopiFoN2*r*r;
	B = 2*twopiFoN2*r*r1*r1 + 2*twopiFoN2/F*dFdr*r*r*r1*r1;
	C = twopiFoN2^2*r*r*r1^4 + 4*twopiFoN2*r1^4 + (2*pi/N*dFdr)^2*r*r*r1^4 + 4*twopiFoN2/F*dFdr*r*r1^4 - (gamma)^2*smax^2;
    
	[rts] = qdf(A,B,C);	% qdf = Quadratic Formula Solution.
	r2 = rts(1);	% Use bigger root - this one corresponds to increasing 
                    % r2 which means spiralling out
end

slew = 1/gamma * (r2-twopiFoN2*r*r1^2 + 1i*twopiFoN*(2*r1^2 + r*r2 + dFdr/F*r*r1^2));
sr = abs(slew)/smaxOrig;
if (abs(slew)/smaxOrig > 1.0001)
    fprintf(1,'Slew violation, slew = %d, smax = %d, sr=%f, r=%f, r1=%f',round(abs(slew)),round(smaxOrig),sr,r,r1);
end

% Calculate minimum possible slew rate for current point
% D = B^2 - 4*A*C >= 0, otherwise slew is not feasible
A = 1+twopiFoN2*r*r;
B = 2*twopiFoN2*r*r1*r1 + 2*twopiFoN2/F*dFdr*r*r*r1*r1;
C = twopiFoN2^2*r*r*r1^4 + 4*twopiFoN2*r1^4 + (2*pi/N*dFdr)^2*r*r*r1^4 + 4*twopiFoN2/F*dFdr*r*r1^4;
smin = sqrt(C-(B^2)/(4*A))/gamma; 

%	Calculate q2 from other pararmeters.
q2 = 2*pi/N*dFdr*r1^2 + 2*pi*F/N*r2;



function dG = findDeltaG(smax,r,r1,T,N,Fcoeff,rmax)

gamma = 4258;			% Hz/G

F = 0;		% FOV function value for this r.
dFdr = 0;		% dFOV/dr for this value of r.
for rind = 1:length(Fcoeff)
	F = F+Fcoeff(rind)*(r/rmax)^(rind-1);
	if (rind>1)
		dFdr = dFdr + (rind-1)*Fcoeff(rind)*(r/rmax)^(rind-2)/rmax;
    end
end

twopiFoN = 2*pi*F/N;
twopiFoN2 = twopiFoN^2;

A = 1 + twopiFoN2*r*r;
B = 2*twopiFoN2*r*r1*r1 + 2*twopiFoN2/F*dFdr*r*r*r1*r1;
C = twopiFoN2^2*r*r*r1^4 + 4*twopiFoN2*r1^4 + (2*pi/N*dFdr)^2*r*r*r1^4 + 4*twopiFoN2/F*dFdr*r*r1^4 - (gamma*smax)^2;
[rts] = qdf(A,B,C);	% qdf = Quadratic Formula Solution.
r2 = rts;

rnew = r + (r1+r2*T)*T;
dG = (r1 + r2*T)/gamma .* sqrt(1+(twopiFoN*rnew).^2);



function [roots] = qdf(a,b,c)

d = b^2 - 4*a*c;
if(d<0)
    d = 0;
end
roots(1) = (-b + sqrt(d))/(2*a);
roots(2) = (-b - sqrt(d))/(2*a);
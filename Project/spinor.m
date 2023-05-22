function [g_out] = spinor(varargin)
%   spinor.m -- plots a circle with a line pointing in a particular
%   angular direction.  Here are the mandatory inputs:
%
%   x -- vector of x positions 
%   y -- vector of y positions
%   mass -- vector of masses.
%   angl -- vector of angles.
%
%   spinor creates a circle of radius=(mass)^1/2.  In this manner, area is
%   preserved as per mass.  The circle is plotted at x and y, and a line is
%   plotted from the center of the circle to the radius, with an angle of
%   angl degrees (note that angl=0 --> point along +x)
%
    g_out=[];
    if(nargin~=4)
        error('spinor() requires the use of four arguments: x, y, mass, angl.\n');
    else
        nsde=24;
        cang=0:2*pi/nsde:2*pi;
        x=varargin{1};
        y=varargin{2};
        mass=varargin{3};
        angl=varargin{4};
        
        try
            matr=[x(:); y(:); mass(:); angl(:)];
        catch
            error('all four inputs must be an array of the same size!\n');
        end

        crad=sqrt(mass);

        for il=1:numel(x)
            plot(crad(il).*cos(cang)+x(il),crad(il).*sin(cang)+y(il),'k','LineWidth',2); hold on;
            plot([x(il) x(il)+crad(il).*cos(angl(il))], ...
                 [y(il) y(il)+crad(il).*sin(angl(il))],'r','LineWidth',2);
        end
        hold off;

    end
end
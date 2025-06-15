function c = saffrongreen(m)
%SAFFRONGREEN Shades of saffron to green colormap
%   SAFFRONGREEN(M) returns an M-by-3 colormap array going from
%   deep green ([0.0157, 0.4157, 0.2196]) to white to saffron 
%   ([1.0, 0.4039, 0.1216]).
%
%   Example:
%       colormap(saffrongreen)

if nargin < 1, m = size(get(gcf,'colormap'),1); end

saffron = [1.0, 0.4039, 0.1216];
green   = [0.0157, 0.4157, 0.2196];
white   = [1, 1, 1];

if mod(m, 2) == 0
    m1 = m / 2;
    
    % Interpolate from green to white
    r1 = linspace(green(1), white(1), m1)';
    g1 = linspace(green(2), white(2), m1)';
    b1 = linspace(green(3), white(3), m1)';
    
    % Interpolate from white to saffron
    r2 = linspace(white(1), saffron(1), m1)';
    g2 = linspace(white(2), saffron(2), m1)';
    b2 = linspace(white(3), saffron(3), m1)';
    
    r = [r1; r2];
    g = [g1; g2];
    b = [b1; b2];
else
    m1 = floor(m / 2);
    
    % Interpolate from green to white
    r1 = linspace(green(1), white(1), m1)';
    g1 = linspace(green(2), white(2), m1)';
    b1 = linspace(green(3), white(3), m1)';
    
    % Center point is white
    r_mid = white(1);
    g_mid = white(2);
    b_mid = white(3);
    
    % Interpolate from white to saffron
    r2 = linspace(white(1), saffron(1), m1)';
    g2 = linspace(white(2), saffron(2), m1)';
    b2 = linspace(white(3), saffron(3), m1)';
    
    r = [r1; r_mid; r2];
    g = [g1; g_mid; g2];
    b = [b1; b_mid; b2];
end

c = [r g b];


function ma = get_ma_vector( ...
    m2,m3,m4,m5,m6, ...
    ag2x,ag2y,ag3x,ag3y,ag4,ag5x,ag5y,ag6x,ag6y, ...
    ddtheta3, ddtheta6, ...
    Ig3, Ig6)

% Matrix defined by equations underlined in blue in the assignment in order
ma = [ 
    m2*ag2x;
    m2*ag2y;
    0;
    m3*ag3x;
    m3*ag3y;
    Ig3*ddtheta3;
    m4*ag4;
    0;
    m5*ag5x;
    m5*ag5y;
    m6*ag6x;
    m6*ag6y;
    Ig6*ddtheta6
    ];

end
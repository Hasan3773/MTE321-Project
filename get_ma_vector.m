function ma = get_ma_vector(%m_i, ... % these are the examples of the possible input
        % ri ... % Only include the inputs that are necessary
        % theta_i ...
        % dtheta_i ...
        % ddtheta_i ...
        % ddr_i, ...
        % I_Gi);

% ENTER YOUR CODE HERE %

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
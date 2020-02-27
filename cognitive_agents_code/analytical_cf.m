% Analytical approach and "light model" for cognitive force. Suppose for
% the sake of this analysis, only friction and cognitive force for model
x0 = 0;
y0 = 0;
x1 = 1;
y1 = 2;

r0 = 5;
r1 = 1;

pi;

dif = sqrt( (x0-x1)^2 + (y0-y1)^2 );

dthetadx = 1 / sqrt( 1 - r1^2/dif^2  ) * (-r1) * (x0-x1) * (dif^(-3));
dthetady = 1 / sqrt( 1 - r1^2/dif^2  ) * (-r1) * (y0-y1) * (dif^(-3));

doadx = (x0-x1) / sqrt( (x0-x1)^2 + (y0-y1)^2  - r1^2);
doady = (y0-y1) / sqrt( (x0-x1)^2 + (y0-y1)^2  - r1^2);

fx = 0.5* 2* (log(2*pi) +1) * ( (r1^2-r0^2) *dthetadx + r1*doadx);
fy = 0.5* 2* (log(2*pi) +1) * ( (r1^2-r0^2) *dthetady + r1*doady);

cf = analytical_cog_force([x0,y0;x1,y1], 5, 1)





function cog_f = analytical_cog_force(agent_coor, cog_size, radius)
    cog_f = zeros( size(agent_coor,1), 2);
    
    xij = agent_coor(:,1) - agent_coor(:,1)' ;
    yij = agent_coor(:,2) - agent_coor(:,2)' ;
    rij2 =  xij.^2 + yij.^2  ;
%         if (rij ~= 0 && rij <= cog_size)
    dthetadx = -radius* (rij2 .^ (-1.5)) ./ sqrt(1-radius^2./rij2) .* xij;
    dthetady = -radius* (rij2 .^ (-1.5)) ./ sqrt(1-radius^2./rij2) .* yij;

    dthetadx(rij2 >= cog_size^2) = 0;
    dthetadx(rij2 == 0) =0;
    dthetady(rij2 >= cog_size^2) = 0;
    dthetady(rij2 == 0) =0;

    doadx = xij ./ sqrt( rij2 - radius^2);
    doady = yij ./ sqrt( rij2 - radius^2);
    doadx(rij2 >= cog_size^2) = 0;
    doadx(rij2 == 0) =0;
    doady(rij2 >= cog_size^2) = 0;
    doady(rij2 == 0) =0;

    fx = 0.5* 2* (log(2*pi) +1) * ( (radius^2 - cog_size^2) * sum(dthetadx,2) + radius * sum(doadx,2));
    fy = 0.5* 2* (log(2*pi) +1) * ( (radius^2 - cog_size^2) * sum(dthetady,2) + radius * sum(doady,2));
    cog_f = [fx fy]/(cog_size^2 * pi);
   
    
end


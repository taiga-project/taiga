function multirun_angled_det(r)
    d = dir(['../results/',r]);
    
    det_mid = [0.7 0.2 0];
    det_ang = [1 0]; %z/R tor/R   
    
    rad = [];
    z = [];
    tor = [];
    
    for i = 3:length(d)
        if d(i).isdir
            t_rad = load([d(i).folder,'/',d(i).name,'/t_rad.dat']);
            t_z = load([d(i).folder,'/',d(i).name,'/t_z.dat']);
            t_tor = load([d(i).folder,'/',d(i).name,'/t_tor.dat']);
        
            t_plane = (t_rad-det_mid(1)) + (t_z-det_mid(2))*det_ang(1) + (t_tor-det_mid(3))*det_ang(2);
        
        
            ix = nan(1,size(t_rad,2));
        
            for j = 1:size(t_rad,2)            
                ix(j) = find((t_plane(:,j)>0),1);
            end
            rad = [rad, mx_interp(t_rad,t_plane,ix)];            
            z   = [z,   mx_interp(t_z,  t_plane,ix)];
            tor = [tor, mx_interp(t_tor,t_plane,ix)];
            
        end
    end
    
    plane = (rad-det_mid(1)) + (z-det_mid(2))*det_ang(1) + (tor-det_mid(3))*det_ang(2);
    keyboard    

end


function I = mx_interp(A,B,ix);

    s = size(A);
    n = ix + s(1).*(0:s(2)-1);
    I = A(n-1) - ( B(n)  ) ./ ( B(n) - B(n-1) ) .* ( A(n) - A(n-1) );

end

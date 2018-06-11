function multirun_angled_det(r)
    d = dir(['../results/',r]);
    
    det_mid = [0.8 0.2 0];
    det_ang = [1 0]; %z/R tor/R   
    
    
    for i = 3:length(d)
        if d(i).isdir
            t_rad = load([d(i).folder,'/',d(i).name,'/t_rad.dat']);
            t_z = load([d(i).folder,'/',d(i).name,'/t_z.dat']);
            t_tor = load([d(i).folder,'/',d(i).name,'/t_tor.dat']);
        end
    end
    keyboard    

end

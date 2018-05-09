
    mainfolder = '../results';
    shotnumber = 'test_1';
    runnumber = '120';
    
    load([mainfolder,'/',shotnumber,'/',runnumber,'/t_rad.dat']);
    load([mainfolder,'/',shotnumber,'/',runnumber,'/t_z.dat']);
  %  load([mainfolder,'/',shotnumber,'/',runnumber,'/t_tor.dat']);
    
    t_rad = t_rad';
    t_z = t_z';
    t_tor = t_tor';
    
    L1=t_rad(:,2:end)<t_rad(:,1);
    L2=t_rad(:,1:end-1)<t_rad(:,1);
    [I1,I2] = find(L1.*L2);
    I = find(L1.*L2);
    
    zi = interp1([t_rad(I) t_rad(I+1)],[t_z(I) t_z(I+1)],t_rad(:,1));
    
    
    
    

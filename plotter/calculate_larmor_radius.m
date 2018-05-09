
    mainfolder = '../results';
    shotnumber = 'test_1';
    runnumber = '120';
    
    load([mainfolder,'/',shotnumber,'/',runnumber,'/t_rad.dat']);
    load([mainfolder,'/',shotnumber,'/',runnumber,'/t_z.dat']);
    load([mainfolder,'/',shotnumber,'/',runnumber,'/t_tor.dat']);
    
    L1=t_rad(2:end,:)<t_rad(1,:);
    L2=t_rad(1:end-1,:)<t_rad(1,:);
    I = find(L1.*L2);
    t_rad(I)
    
    
    
    
    

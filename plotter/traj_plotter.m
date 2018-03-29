function traj_plotter(varargin)


    mainfolder = '../results';
    
    
    if nargin >= 1
        shotnumber = varargin{1};
    else
        shotnumber = '11774_1000';
      
    end
    if nargin >= 2
        runnumber = varargin{2};
    else        
        runnumber = '0';
    end
    
    load([mainfolder,'/',shotnumber,'/',runnumber,'/t_rad.dat']);
    load([mainfolder,'/',shotnumber,'/',runnumber,'/t_z.dat']);
    load([mainfolder,'/',shotnumber,'/',runnumber,'/t_tor.dat']);
    
    figure
    subplot(2,2,1)
    plot(t_tor,t_z)
    xlabel('tor')
    ylabel('z')
    
    
    subplot(2,2,2)
    plot(t_rad,t_z)
    xlabel('R')
    ylabel('z')
    
    
    subplot(2,2,3)
    plot(t_tor,t_rad)
    xlabel('tor')
    ylabel('R')
    
    
    subplot(2,2,4)
    plot3(t_rad,t_z,t_tor)
    xlabel('R')
    ylabel('z')
    zlabel('tor')
    
end

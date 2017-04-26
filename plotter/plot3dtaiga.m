% Plot TAIGA trajectories


mainfolder = '../results';
fieldfolder = '../input/fieldGrid';
energy=60;
shotnumber = '11344_0285'
runnumber = '06May2016_140034'



load([mainfolder,'/',shotnumber,'/',runnumber,'/t_rad.dat'])
load([mainfolder,'/',shotnumber,'/',runnumber,'/t_z.dat'])
load([mainfolder,'/',shotnumber,'/',runnumber,'/t_tor.dat'])

load([fieldfolder,'/',shotnumber,'/',runnumber,'/psi2.dat'])
load([fieldfolder,'/',shotnumber,'/',runnumber,'/rcord.dat'])
load([fieldfolder,'/',shotnumber,'/',runnumber,'/zcord.dat'])

[R,Z]=meshgrid(rcord,zcord);

%contour(R,Z,psi2',30)
axis equal
hold on
%plot(t_rad,t_z,'r')

t_rad2 = sqrt(1+(t_tor./t_rad).^2).*t_rad;
t_rad
t_rad2
plot3(t_rad,t_z,t_tor,'r')

axis equal

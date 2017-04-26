% Plot TAIGA trajectories
close all
clear all
clc

%! detector position
detpos = 0.7089;

TORLIM = [-0.07, 0.03];
ZLIM = [0.18 0.28];%[0.12, 0.20];
TORLIM = [-0.04 0.04]
ZLIM = [0.20 0.28]

%TORLIM = [-0.03,-0.01]
%ZLIM = [0.24 0.25]
%mainfolder = '../../results';
mainfolder = '../../SOFT';
energy=60;
runnumber = '11347_0285';
%runnumber = '06May2016_140034'
%runnumber = '06May2016_141714'

load([mainfolder,'/',runnumber,'/rad.dat'])
load([mainfolder,'/',runnumber,'/z.dat'])
load([mainfolder,'/',runnumber,'/tor.dat'])

load([mainfolder,'/',runnumber,'/t_rad.dat'])
load([mainfolder,'/',runnumber,'/t_z.dat'])
load([mainfolder,'/',runnumber,'/t_tor.dat'])


s_rad = (t_rad(1,:));
s_z   = (t_z(1,:));
s_tor = (t_tor(1,:));

%ind = find(rad==rad(1,1));
ind = find(rad==detpos);

rions = [0.7,0.68,0.66,0.64];

ind2{1} = find(s_rad<rions(1)+1e-3 & s_rad>rions(1)-1e-3 );
ind2{2} = find(s_rad<rions(2)+1e-3 & s_rad>rions(2)-1e-3 );
ind2{3} = find(s_rad<rions(3)+1e-3 & s_rad>rions(3)-1e-3 );
ind2{4} = find(s_rad<rions(4)+1e-3 & s_rad>rions(4)-1e-3 );

rions2 = [0.71,0.69,0.67,0.65,0.63];

ind22{1} = find(s_rad<rions2(1)+1e-2 & s_rad>rions2(1)-1e-2 );
ind22{2} = find(s_rad<rions2(2)+1e-3 & s_rad>rions2(2)-1e-3 );
ind22{3} = find(s_rad<rions2(3)+1e-3 & s_rad>rions2(3)-1e-3 );
ind22{4} = find(s_rad<rions2(4)+1e-3 & s_rad>rions2(4)-1e-3 );


if false
figure
%for ploti=1:9
   % subplot(3,3,ploti)
    %contour(R,Z,psi2',30)
    hold on
    plot(tor(ind),z(ind),'r.')
    
    for i=1:4	
	    plot(tor(ind2{i}),z(ind2{i}),'ko','MarkerSize',0.5)
	    zs(i) = mean(z(ind2{i}));
	    ts(i) = mean(tor(ind2{i}));
	    zs2(i) = mean(z(ind22{i}));
	    ts2(i) = mean(tor(ind22{i}));
    end
    set(gca,'xdir','reverse')
    
ylim([0.04 0.16])
    axis equal
		title(['$R_\mathrm{detector} = ',num2str(detpos),'$'],'interpreter','latex','fontsize',14)

    xlabel('{$T$ (m)}','interpreter','latex','fontsize',14)
    ylabel('{$Z$ (m)}','interpreter','latex','fontsize',14)
    saveas(gcf,'total.pdf')
    saveas(gcf,'total.png')

end

figure
%for ploti=1:9
   % subplot(3,3,ploti)
    %contour(R,Z,psi2',30)
    %[values,categ] = 
    %hist3([z(ind)',tor(ind)'],[21,21]);
    %imagesc(categ{2},categ{1},values)
    %contourf(categ{2},categ{1},values)
    
   % [val,cate]=contourc(z(ind)',tor(ind)')
X=    [tor(ind)',z(ind)'];
rectangle('Position',[TORLIM(1) ZLIM(1) TORLIM(2)-TORLIM(1) ZLIM(2)-ZLIM(1)],'FaceColor','k')
hold on
       smoothhist2D(X,1,[1000, 1000],-0.1);
  set(gca,'xdir','reverse','ydir','normal')
    colorbar('YTick', [0, 64 ,128, 192, 256],'YTickLabel',{'0%','25%','50%','75%','100%'})

    %cb = colorbar('XTick', [0, 64 ,128, 192, 256],'XTickLabel',{'0%','25%','50%','75%','100%'})
%    set(gca,'YDir','normal')

		
xlim0=xlim;ylim0=ylim;
xlim0
ylim0
xlim0=TORLIM;
ylim0=ZLIM;
		%axis equal
xlim(xlim0)
ylim(ylim0)

if true
hold on
valy = ZLIM(1):0.01:ZLIM(2);
valx = [xlim0(2)+0.003,xlim0(2)-0.003];
[valX,valY]=meshgrid(valx,valy)
plot(valX',valY','w')

valy = ZLIM(1):0.005:ZLIM(2);
valx = [xlim0(2)+0.0015,xlim0(2)-0.0015];
[valX,valY]=meshgrid(valx,valy)
plot(valX',valY','w')

valx = TORLIM(1):0.01:TORLIM(2);
valy = [ylim0(1)-0.003,ylim0(1)+0.003];
[valX,valY]=meshgrid(valx,valy)
plot(valX,valY,'w')

valx = TORLIM(1):0.005:TORLIM(2);
valy = [ylim0(1)-0.0015,ylim0(1)+0.0015];
[valX,valY]=meshgrid(valx,valy)
plot(valX,valY,'w')

end
%ts = (ts - xlim0(2))/(xlim0(1)-xlim0(2));
%zs = (zs - ylim0(1))/(ylim0(2)-ylim0(1));


if false
for i=1%:4
	text(ts(i)+0.01/i+0.02,zs(i),['',num2str(rions(i)),' m \rightarrow'],'color','w')
	%annotation('textarrow', [ts(i),ts(i)],[zs(i)-0.01,zs(i)], 'String', ['',num2str(rions(i)),' m \rightarrow']', 'color', 'w');

	%	text(ts2(i)+0.01/i+0.008,zs2(i),['\rightarrow'],'color','w')
end
end

    ttli = [0 strfind(runnumber,'_') length(runnumber)+1];
    runshot = runnumber(ttli(1)+1:ttli(2)-1);
    
    if strcmp(runnumber (ttli(2)+1), 'm')
        if (ttli(2)+3==ttli(3))
            runangle = ['-',runnumber(ttli(2)+2),'.0'];            
        else
            runangle = ['-',runnumber(ttli(2)+2),'.',runnumber(ttli(2)+3:ttli(3)-1)];
        end
    else
        if (ttli(2)+2==ttli(3))
            if strcmp(runnumber(ttli(2)+1),'0')
                runangle = [runnumber(ttli(2)+1),'.0'];           
            else
                runangle = ['+',runnumber(ttli(2)+1),'.0'];
            end
        
        else
            runangle = ['+',runnumber(ttli(2)+1),'.',runnumber(ttli(2)+2:ttli(3)-1)];
        end
    end
	title(['$\#',runshot,...
        '~~R_\mathrm{detector} = ',num2str(detpos),'\mathrm{~m}~~(\varphi_\mathrm{in}=',...
        runangle,'^\circ) E=',num2str(energy),'~keV $'],'interpreter','latex','fontsize',14)
    xlabel('{$T$ (m)}','interpreter','latex','fontsize',14)
    ylabel('{$Z$ (m)}','interpreter','latex','fontsize',14)
    saveas(gcf,['plots/',runnumber,'_detector.pdf'])

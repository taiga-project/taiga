function detPlot(varargin)





    %TORLIM = [-0.03,-0.01]
    %ZLIM = [0.24 0.25]
    mainfolder = '../results';
    energy=60;
    %runnumber = '06May2016_141714'
    
    detector_plotted=false;
    particle_plotted=true;
    
    if nargin >= 1
        shotnumber = varargin{1};
    else
        shotnumber = '11344_1000'        
    end
    if nargin >= 2
        runnumber = varargin{2};
    else        
        runnumber = '06May2016_140034'
    end
    
    if nargin >= 3
        ZMID_DET = varargin{3};
    else
        ZMID_DET = 0.23;
    end
    
    if nargin >= 4
        TMID_DET = varargin{4};
    else
        TMID_DET = 0.00;
    end
    
    
    
    %! detector position
    detpos = 0.7089;

    TORLIM = [TMID_DET-0.04 TMID_DET+0.04];
    ZLIM = [ZMID_DET-0.04 ZMID_DET+0.04];
    

    TORDET = [-0.00425 -0.00225 -0.00175 -0.00075 -0.00025 0.00025 0.00075 0.00175 0.00225 0.00425];
    %                d        g      d       g       d      g      d      g      d


    ZDET = ZMID_DET+[-0.01075 -0.00575 -0.00525 -0.00025 0.00025 0.00525 0.00575 0.01075];
    %                        d        g        d        g       d       g        d   

    

    load([mainfolder,'/',shotnumber,'/',runnumber,'/rad.dat'])
    load([mainfolder,'/',shotnumber,'/',runnumber,'/z.dat'])
    load([mainfolder,'/',shotnumber,'/',runnumber,'/tor.dat'])

    load([mainfolder,'/',shotnumber,'/',runnumber,'/t_rad.dat'])
    load([mainfolder,'/',shotnumber,'/',runnumber,'/t_z.dat'])
    load([mainfolder,'/',shotnumber,'/',runnumber,'/t_tor.dat'])



    if  particle_plotted
        figure
        plot(tor(1:20),z(1:20),'.')
        xlabel('tor')
        ylabel('z')
    end


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
    rectangle('Position',[TORLIM(2) ZLIM(1) TORLIM(2)-TORLIM(1) ZLIM(2)-ZLIM(1)],'FaceColor','k','EdgeColor','none')
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




    if detector_plotted
    hold on


    TORDET = [TORLIM(1), TORDET, TORLIM(2)];
    ZDET = [ZLIM(1), ZDET, ZLIM(2)];

	
	
	    for i=1:2:length(ZDET)
		    rectangle('Position',[TORLIM(1) ZDET(i) TORLIM(2)-TORLIM(1) ZDET(i+1)-ZDET(i)],'FaceColor','k','EdgeColor','none')
	    end
	
	    for i=1:2:length(TORDET)
	    [TORDET(i) ZLIM(1) TORDET(i+1)-TORDET(i) ZLIM(2)-ZLIM(1)]
		    rectangle('Position',[TORDET(i) ZLIM(1) TORDET(i+1)-TORDET(i) ZLIM(2)-ZLIM(1)],'FaceColor','k','EdgeColor','none')
	    end
	
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

    if detector_plotted
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
        runangle,'^\circ) E=',num2str(energy),'~\mathrm{keV} $'],'interpreter','latex','fontsize',14)
    xlabel('{$T$ (m)}','interpreter','latex','fontsize',14)
    ylabel('{$Z$ (m)}','interpreter','latex','fontsize',14)
    mkdir('plots')
    mkdir(['plots/',shotnumber])
    saveas(gcf,['plots/',shotnumber,'/',runnumber,'_detector',num2str(ZMID_DET),'.pdf'])
end

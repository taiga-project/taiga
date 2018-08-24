function detPlot(varargin)


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
    
    
    if nargin >= 3
        energy = varargin{3};
    else        
        energy=60;
    end
    
    
    if nargin >= 4
        runangle = varargin{4};
    else        
        runangle=0;
    end

    if nargin >= 5
        DETECTOR_POS_ANGLE = str2num(varargin{5});
    else
        DETECTOR_POS_ANGLE = [0.7071 0.22 0.0 1.0 0.0] ;
    end
    
    if nargin >= 6
        scenario = varargin{6};
    else        
        scenario = '_all';
    end
    
    
    if nargin >= 7
        ZMID_DET = varargin{7};
    else
        ZMID_DET = 0.22;
    end
    
    if nargin >= 8
        TMID_DET = varargin{8};
    else
        TMID_DET = 0.00;
    end

    
    if nargin >= 9
        detector_plotted = varargin{9};
    else
        detector_plotted = false;
    end
    
    if nargin >= 10
        particle_plotted = varargin{10};
    else
        particle_plotted = false;
    end
    
    
    in.tokamak = 'compass';    
    in.folder = '../input/renate110';
    in.prof_folder = '../input/ionProf';
    in.shotNumber=shotnumber;    
    
    %! detector position    

    det_R = DETECTOR_POS_ANGLE(1); %0.684618;
    det_z = DETECTOR_POS_ANGLE(2);
    det_angle_Rz = DETECTOR_POS_ANGLE(4);
    det_tan_Rz = tan(det_angle_Rz/180*pi);

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





    s_rad = (t_rad(1,:));
    s_z   = (t_z(1,:));
    s_tor = (t_tor(1,:));

    %ind = find(rad==rad(1,1));
    
    %trafo
    
 %   rad = (z-det_z)+rad;
%    z = det_z + (z-det_z)*sqrt(2);
    
   ind = find(abs((rad-det_R)+det_tan_Rz*(z-det_z))<1e-4);
%    ind = find(rad==detpos);
    
    aaa=find(abs(tor(ind))<1);
    
    ind=ind(aaa);
    
    if strcmp(scenario,'_spx')	
        in.filepath = [mainfolder,'/',shotnumber,'/',runnumber,'/detector_',shotnumber,'_spx.pdf'];
  	in.spx = load(in.filepath);
        aaa=find(s_rad(ind)<in.spx-0.01 && s_rad(ind)>in.spx+0.01)		
        ind=ind(aaa);		
    end
	
    figure('visible', 'off')
    h = hist(s_rad,100);    
    hist(s_rad,100); 
    
    try
        ionrate = load([in.prof_folder,'/',shotnumber,'/ionrate.dat']);
        ionrad = load([in.prof_folder,'/',shotnumber,'/rad.dat']);
        
        hold on        
        plot(ionrad,ionrate/max(ionrate)*max(h),'r')
    end
    saveas(gcf,[mainfolder,'/',shotnumber,'/',runnumber,'/ionprof',scenario,'.pdf'])

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
    figure('visible', 'off')
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
		    title(['$z_\mathrm{detector} = ',num2str(det_z),'$'],'interpreter','latex','fontsize',14)

        xlabel('{$T$ (m)}','interpreter','latex','fontsize',14)
        ylabel('{$Z$ (m)}','interpreter','latex','fontsize',14)
        saveas(gcf,'total.pdf')
        saveas(gcf,'total.png')

    end

    figure('visible', 'off')
    %for ploti=1:9
       % subplot(3,3,ploti)
        %contour(R,Z,psi2',30)
        %[values,categ] = 
        %hist3([z(ind)',tor(ind)'],[21,21]);
        %imagesc(categ{2},categ{1},values)
        %contourf(categ{2},categ{1},values)
        
       % [val,cate]=contourc(z(ind)',tor(ind)')
       
    ind000=ind;
    a000=find(tor(ind000)<TORLIM(2));
    ind000=ind000(a000);
    a000=find(tor(ind000)>TORLIM(1));
    ind000=ind000(a000);
    a000=find(z(ind000)<ZLIM(2));
    ind000=ind000(a000);
    a000=find(z(ind000)>ZLIM(1));
    ind000=ind000(a000);
    
    X=    [tor(ind000)',z(ind000)'];
   % X
    %rectangle('Position',[TORLIM(2) ZLIM(1) TORLIM(2)-TORLIM(1) ZLIM(2)-ZLIM(1)],'FaceColor','k','EdgeColor','none')
    hold on
    X = [X; TORLIM(2) ZLIM(1); TORLIM(1) ZLIM(2)];
           smoothhist2D(X,1,[150, 150],-0.1);
      set(gca,'xdir','reverse','ydir','normal')
        colorbar('YTick', [0, 64 ,128, 192, 256],'YTickLabel',{'0%','25%','50%','75%','100%'})

        %cb = colorbar('XTick', [0, 64 ,128, 192, 256],'XTickLabel',{'0%','25%','50%','75%','100%'})
    %    set(gca,'YDir','normal')

		
    xlim0=xlim;ylim0=ylim;
    %xlim0
    %ylim0
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
    %runshot = runnumber(ttli(1)+1:ttli(2)-1);
    shotnumber_array = strsplit(shotnumber,'_');
    detector_angle_RZ = round(atan2(abs(rad(ind(1))-rad(ind(end))),abs(z(ind(1))-z(ind(end))))/pi*180,1);

    runshot = [shotnumber_array{1},'~~(t = ',shotnumber_array{2},'\mathrm{~s,~detector~angle:~}',num2str(detector_angle_RZ),'^\circ)'];

%	title(['$\#',runshot,...
%        '~~z_\mathrm{detector} = ',num2str(det_z),'\mathrm{~m}~~(\varphi_\mathrm{in}=',...
%        num2str(runangle),'^\circ) E=',num2str(energy),'~\mathrm{keV} $'],'interpreter','latex','fontsize',14)



    title(['$\#',runshot,' $'],'interpreter','latex','fontsize',14)
    xlabel('{$T$ (m)}','interpreter','latex','fontsize',14)
    ylabel('{$Z$ (m)}','interpreter','latex','fontsize',14)
    xlim0=xlim;
    xticks(xlim0(1):0.01:xlim0(2));
    yticks(ylim0(1):0.01:ylim0(2));
    
    mkdir('plots')
    mkdir(['plots/',shotnumber])
    %saveas(gcf,['plots/',shotnumber,'/',runnumber,'_detector',num2str(ZMID_DET),'.pdf'])
    disp(['Saved to ',mainfolder,'/',shotnumber,'/',runnumber,'/detector_',shotnumber,'.pdf'])
    saveas(gcf,[mainfolder,'/',shotnumber,'/',runnumber,'/detector_',shotnumber,'.pdf'])    
    saveas(gcf,[mainfolder,'/',shotnumber,'/',runnumber,'/detector_',shotnumber,'.fig'])


    xt = xticklabels;
    yt = yticklabels;
    
    for i=1:length(xt)
        xt{i} = num2str(str2num(xt{i})*100);
    end
    for i=1:length(yt)
        yt{i} = num2str(str2num(yt{i})*100);
    end
    
    xticklabels(xt)
    yticklabels(yt)
    xlabel('{vertical position (cm)}','interpreter','latex','fontsize',14)
    ylabel('{horizontal position (cm)}','interpreter','latex','fontsize',14)

    
    disp(['Saved to ',mainfolder,'/',shotnumber,'/',runnumber,'/detector_',shotnumber,'_cm.pdf'])
    saveas(gcf,[mainfolder,'/',shotnumber,'/',runnumber,'/detector_',shotnumber,'_cm.pdf'])    
    saveas(gcf,[mainfolder,'/',shotnumber,'/',runnumber,'/detector_',shotnumber,'_cm.fig'])
    
    
    figure('visible', 'off')
    plot(rad(ind),z(ind))
    xlabel('R [m]')
    ylabel('z [m]')
    
    saveas(gcf,[mainfolder,'/',shotnumber,'/',runnumber,'/RZ',scenario,'.pdf'])
end

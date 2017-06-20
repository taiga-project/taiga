% IDL output reader

function renate110_to_taiga(varargin)
    
    in.tokamak = 'compass';
    
    in.folder = '../input/renate110';
    out.folder = '../input/ionProf';
    out.plotfolder='../results';
    
	in.index.rad = 1;
	in.index.ionyeald = 3;
    in.ionrate=false;
    in.plot=false;

    if nargin == 0
        in.shotNumber = '11774';
        in.time=1099;
    elseif nargin == 1
        in.shotNumber = '11774'; 
        in.time = varargin{1};
    else    
        in.shotNumber = varargin{1};
        in.time = varargin{2};
    end    
    
    if nargin >=3    
        renate.angle = varargin{3};
	else
	    renate.angle = 0.0;
	end
    
    if nargin >=4    
        renate.rad.max = varargin{3};
	else
	    renate.rad.max = 0.78;
	end
    
    in.filepath = [in.folder,'/pop_',in.tokamak,in.shotNumber,'_',num2str(in.time),'.txt'];
    
  	data = load(in.filepath);

	mkdir(out.folder)
	foldername = ([out.folder,'/', in.shotNumber,'_',num2str(in.time),'/']);
	mkdir(foldername)
	
	out.rad      = renate.rad.max-data(:,in.index.rad);
	out.ionyeald = data(:,in.index.ionyeald);
	if in.ionrate
	    v = 1:length(ionyeald);
	    d = 1e-6;
	    out.ionrate = (interp1(v,out.ionyeald,v+d/2)-interp1(v,out.ionyeald,v-d/2) )/d;
	end
	
	dlmwrite([foldername, 'rad.dat']     , out.rad     , 'precision','%.16e','delimiter','\t');  
	dlmwrite([foldername, 'ionyeald.dat'], out.ionyeald, 'precision','%.16e','delimiter','\t'); 
    if in.ionrate
	    dlmwrite([foldername, 'ionrate.dat'], out.ionyeald, 'precision','%.16e','delimiter','\t'); 
    end	
	
	if in.plot
	    figure
	    plot(out.rad,1-out.ionyeald,'r','linewidth',2)
	    ylim([0 1])   
        title (['Ionisation yield from RENATE @ ', upper(in.tokamak), ' #', in.shotNumber, ' (', num2str(in.time),' s) angle: ',num2str(renate.angle),' degree'])
        xlabel('R [m]')
        ylabel('ionisation yield')
	    saveas(gcf,[out.plotfolder,'/', in.shotNumber,'_',num2str(in.time),'/renate_',num2str(renate.angle),'.pdf'])
	
        if in.ionrate		
	        figure
	        plot(out.rad,out.ionrate,'r','linewidth',2)
	        ylim([0 1])   
            title (['Ion population from RENATE @ ', upper(in.tokamak), ' #', in.shotNumber, ' (', num2str(in.time),' s) angle: ',num2str(renate.angle),' degree'])
            xlabel('R [m]')
            ylabel('ion population')
	        saveas(gcf,[out.plotfolder,'/', in.shotNumber,'_',num2str(in.time),'/renate_',num2str(renate.angle),'_pop.pdf'])
        end
    end
end

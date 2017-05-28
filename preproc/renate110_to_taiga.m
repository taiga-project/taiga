% IDL output reader

function renate110_to_taiga(varargin)

    in.tokamak = 'compass';
    
    in.folder = '../input/renate110';
    out.folder = '../input/ionProf';
    out.plotfolder='../results';
    
	in.index.rad = 1;
	in.index.ionyeald = 3;

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
	
	dlmwrite([foldername, 'rad.dat']     , out.rad     , 'precision','%.16e','delimiter','\t');  
	dlmwrite([foldername, 'ionyeald.dat'], out.ionyeald, 'precision','%.16e','delimiter','\t'); 
	ylim([0 1])   
    title (['Ionisation yield from RENATE @ ', upper(in.tokamak), ' #', in.shotNumber, ' (', num2str(in.time),' s) angle: ',num2str(renate.angle),' ^o'])
    xlabel('R [m]')
    ylabel('ionisation yield')
	saveas(gcf,[out.plotfolder,'/', in.shotNumber,'_',num2str(in.time),'/renate_',num2str(renate.angle),'.pdf'])

end

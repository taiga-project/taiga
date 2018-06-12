% IDL output reader

function pavel_profile_to_taiga(varargin)
    
    in.tokamak = 'compass';
    
    in.folder = '../input/ionProf/pavel';
    out.folder = '../input/ionProf';
    out.plotfolder='../results';
    
	in.index.rad = 1;
	in.index.ionyeald = 3;
    in.ionrate=true;
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
    

    
    in.filepath = [in.folder,'/', in.shotNumber,'_',num2str(in.time),'_ionization.mat'];
    
  	data = load(in.filepath);

	mkdir(out.folder)
	foldername = ([out.folder,'/', in.shotNumber,'_',num2str(in.time),'/']);
	mkdir(foldername)
	
	out.rad      = data.R';
	out.ionrate = data.ion_prob';
	
	out.ionyeald = 1-cumsum(out.ionrate)/max(cumsum(out.ionrate));
	
	dlmwrite([foldername, 'rad.dat']     , out.rad     , 'precision','%.16e','delimiter','\t');  
	dlmwrite([foldername, 'ionyeald.dat'], out.ionyeald, 'precision','%.16e','delimiter','\t'); 
    if in.ionrate
	    dlmwrite([foldername, 'ionrate.dat'], out.ionrate, 'precision','%.16e','delimiter','\t'); 
    end	
	

end

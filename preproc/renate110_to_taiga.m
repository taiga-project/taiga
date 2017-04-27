% IDL output reader

function renate110_to_taiga(varargin)

    in.tokamak = 'compass';
    
    in.folder = '../input/renate110';
    out.folder = '../input/ionProf';
    
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
	 

end



function old
	runshot = '11344';

	runmode.folder.renate='.';
	runmode.renate.rad.max=0.78
	runmode.renate.rad.min=0.6
	runmode.renate.number=400
	runmode.renate.index.r=1;
	runmode.renate.index.dens=3;
	runmode.renate.index.pop=3;

	files = dir([runmode.folder.renate,'/output/*.txt']);
		     
	%% read RENATE output
	data = load ([runmode.folder.renate,'/output/',files(length(files)).name]);

	%% ionisation distribution
	runmode.renate.distr = 1-cumsum(data(:,runmode.renate.index.pop))/sum(data(:,runmode.renate.index.pop));

	%% ionisation density function
	runmode.renate.profil = (data(1:end-1,runmode.renate.index.dens)-data(2:end,runmode.renate.index.dens))/(runmode.renate.rad.max-runmode.renate.rad.min)*(runmode.renate.number);
	runmode.renate.profil=[runmode.renate.profil;0];

	runmode.renate.profil2 = data(:,runmode.renate.index.pop);


	radial_coords=linspace(runmode.renate.rad.max,runmode.renate.rad.min,length(runmode.renate.distr));
	close all
	subplot(2,1,1)
	plot(radial_coords,runmode.renate.distr,'r','linewidth',2)


	title(['$\#',runshot,'$'],'interpreter','latex','fontsize',14)
	xlabel('{$R$ (m)}','interpreter','latex','fontsize',14)
	ylabel('{cumulative ionisation rate}','interpreter','latex','fontsize',14)

	%figure
	%hold on

	%plot(radial_coords,runmode.renate.profil,'r')

	subplot(2,1,2)
	plot(radial_coords,runmode.renate.profil2,'r','linewidth',2)

	title(['$\#',runshot,'$'],'interpreter','latex','fontsize',14)
	xlabel('{$R$ (m)}','interpreter','latex','fontsize',14)
	ylabel('{normalised ionisation rate}','interpreter','latex','fontsize',14)
	saveas(gcf,['plots/renate_',runshot,'.pdf'])
end

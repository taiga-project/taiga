% //! read hdf5 data files

function cdb_reader(varargin)
    startclean

    in.tokamak = 'compass';
    in.majorradius=0.56;  
    
    in.plot = false;
    
    in.folder = '../input/cdb';
    out.folder.renate = '../input/renate110';
    out.folder.grid = '../input/fieldGrid/';
    out.folder.spline = '../input/fieldSpl/';
    out.folder.plot = '../results/';
    

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
          
    out = getHiddenParameters(out);
    
    %   start    
    [efit, ts] = initEmpty;
    in = readTimes(in);
    [efit, out] = readEfitGrid(in, out, efit);
    
    efit = readMagneticAxis(in, efit);
    efit = readMagneticBoundary(in, efit);
    
    efit = readToroidalFlux(in, out, efit);
    efit = readPoloidalFlux(in, out, efit);
    
    efit = makeMagneticGrid (in, out, efit);
    out = normaliseFlux(in, out, efit);
    
    ts = readThomsonData(in, out, ts);
        
    out.efit.z      = linspace(min(ts.z),max(ts.z),200);
    out.efit.r      = ones(size(out.efit.z))*in.majorradius;
    
    out = fitProfilesNT(ts, out);
            
    saveRenateFlux(in, out);    
    saveRenateNT (in, out);
    
    saveMagneticGrid (in, out, efit);    
    saveMagneticSpline (in, out, efit);

    if in.plot
        plotProfilesNT (in, out, ts);        
        plotNormFlux (in, out, efit);
    end

end

function efit = makeMagneticGrid (in, out, efit)
    psi_RZ = efit.polflux;
    R_M = out.flux.r;
    Z_M = out.flux.z;
    dx=1e-6;
    psi_dR1 = interp2(R_M,Z_M,psi_RZ,R_M-dx/2,Z_M,'spline');
    psi_dR2 = interp2(R_M,Z_M,psi_RZ,R_M+dx/2,Z_M,'spline');
    psi_dZ1 = interp2(R_M,Z_M,psi_RZ,R_M,Z_M-dx/2,'spline');
    psi_dZ2 = interp2(R_M,Z_M,psi_RZ,R_M,Z_M+dx/2,'spline');

    psi_dR = (psi_dR2-psi_dR1)/dx;
    psi_dZ = (psi_dZ2-psi_dZ1)/dx;


    

    brad = -psi_dZ./R_M;%*pi;
    bz   =  psi_dR./R_M;%*pi;	
    efit.brad=brad';
    efit.bz=bz';
end

function saveMagneticGrid (in, out, efit)    
    
    foldername = ([out.folder.grid,'/', in.shotNumber,'_',num2str(in.time),'/']);
    mkdir(foldername)
    
    dlmwrite([foldername, 'rcord.dat'],efit.r, 'precision','%.16e','delimiter','\t');
    dlmwrite([foldername, 'zcord.dat'],efit.z, 'precision','%.16e','delimiter','\t');
    dlmwrite([foldername, 'brad.dat'],efit.brad, 'precision','%.16e','delimiter','\t');
    dlmwrite([foldername, 'bz.dat'],efit.bz, 'precision','%.16e','delimiter','\t');
    dlmwrite([foldername, 'btor.dat'],efit.btor, 'precision','%.16e','delimiter','\t');    
    dlmwrite([foldername, 'psi2.dat'],efit.polflux, 'precision','%.16e','delimiter','\t');    
    
    % a btor szar
    
end


function savePlot (in, out, plotname)
    mkdir(out.folder.plot)
    foldername = [out.folder.plot,'/', in.shotNumber,'_',num2str(in.time),'/'];
    mkdir(foldername)
    saveas(gcf,[foldername, plotname, '.pdf'])
    close
end

function saveMagneticSpline (in, out, efit)
      
	complist = {'brad','bz','btor'};
    for i=1:length(complist)
        comp=complist{i};
        sp = csapi({efit.r,efit.z},efit.(comp));

        i11=1:sp.pieces(1);
        i12=i11+1*sp.pieces(1);
        i13=i11+2*sp.pieces(1);
        i14=i11+3*sp.pieces(1);

        i21=1:sp.pieces(2);
        i22=i21+1*sp.pieces(2);
        i23=i21+2*sp.pieces(2);
        i24=i21+3*sp.pieces(2);


        c=reshape(sp.coefs,size(sp.coefs,2),[]);
        b = sp.breaks{1};
        b2 = sp.breaks{2};

	    mkdir(out.folder.spline)
	    foldername = [out.folder.spline,'/', in.shotNumber,'_',num2str(in.time),'/'];
	    mkdir(foldername)
        save([foldername,'/r.spline'],'-ascii','-tabs','-double','b')
        save([foldername,'/z.spline'],'-ascii','-tabs','-double','b2')
        disp('Spline calculated')        
        
        c11 = c(i11,i21);    c12 = c(i11,i22);    c13 = c(i11,i23);    c14 = c(i11,i24);
        c21 = c(i12,i21);    c22 = c(i12,i22);    c23 = c(i12,i23);    c24 = c(i12,i24);
        c31 = c(i13,i21);    c32 = c(i13,i22);    c33 = c(i13,i23);    c34 = c(i13,i24);
        c41 = c(i14,i21);    c42 = c(i14,i22);    c43 = c(i14,i23);    c44 = c(i14,i24);
        
        for fi1 = 1:4
        	for fi2 = 1:4        	
		        save([foldername,'/',comp,'.spl',num2str(fi1),num2str(fi2)],'-ascii','-tabs','-double',['c',num2str(fi1),num2str(fi2)])
            end
        end
        disp(comp)
        disp('Spline saved')
    end
end


function out = getHiddenParameters(out)
    out.nt.zeff = 2.0;
    out.nt.psi  = 0:.01:2;
    out.nt.psi_in  = 0:.01:1.25;
end

function [efit, ts] = initEmpty
    efit = '';
    ts = '';
end

function in = readTimes(in)

    % EFIT
    in.source = 'efitxx';
    in.hdf5flag = '/time';
    in.index.efit=find(readDataFile(in)>=in.time/1000,1);
    
    % TS 
    in.source = 'ts';
    in.hdf5flag = 'TS_record_time';
    in.index.ts=find(readDataFile(in)>=in.time,1);
    
end

function ts = readThomsonData(in, out, ts)

    in.source = 'ts';
    % TS density profile (from TS)
    in.hdf5flag = 'ne';
    ts.densityRaw = readVectorData(in);
    in.hdf5flag = 'ne_err';
    ts.densityRawErr = readVectorData(in);
        
    % TS coordinates
    in.hdf5flag = 'TS_z_axis';
    ts.zRaw = readDataFile(in);
    
    % TS temperature profile (from TS)
    in.hdf5flag = 'Te';
    ts.temperatureRaw = readVectorData(in);
    in.hdf5flag = 'Te_err';
    ts.temperatureRawErr = readVectorData(in);

    %value to renate110
    ts.densityRaw        = ts.densityRaw        / 1e19;
    ts.densityRawErr     = ts.densityRawErr     / 1e19;
    ts.temperatureRaw    = ts.temperatureRaw    / 1000;
    ts.temperatureRawErr = ts.temperatureRawErr / 1000;

    notnanindex      = (~isnan(ts.densityRaw));
    ts.density       = [ts.densityRaw(notnanindex)];        % 1e19 m-3
    ts.densityErr    = [ts.densityRawErr(notnanindex)];     % 1e19 m-3
    ts.temperature   = [ts.temperatureRaw(notnanindex)];    % keV
    ts.temperatureErr= [ts.temperatureRawErr(notnanindex)]; % keV
    ts.z             = [ts.zRaw(notnanindex)];
    ts.r             = ones(size(ts.z))*in.majorradius;
    ts.psi           = interp2(out.flux.r,out.flux.z,out.flux.normPolFlux,ts.r,ts.z);
end


function efit = readToroidalFlux(in, out, efit)
    in.source = 'efitxx';
    
    %B*Rgeom
    in.hdf5flag = '/output/globalParameters/bvacRgeom';
    efit.bvacrgeom = readScalarData(in);
    
    in.hdf5flag = '/output/fluxFunctionProfiles/rBphi';
    efit.rbtor = readVectorData(in);
    
    efit.btor = -kron(efit.rbtor ./ efit.r, ones(1,length(efit.z)));

end

function [efit, out] = readEfitGrid(in, out, efit)

    in.source = 'efitxx';
    
    % EFIT major radius (R)
    in.hdf5flag = '/output/profiles2D/r';
    efit.r       = readVectorData(in);
    
    % EFIT vertical position (Z)
    in.hdf5flag = '/output/profiles2D/z';
    efit.z       = readVectorData(in);
    
    % size of R and Z
    efit.R_size    = length(efit.r);
    efit.Z_size    = length(efit.z);
    
    
    % RZ meshgrid
    [out.flux.r,out.flux.z] = meshgrid(efit.r,efit.z);

end

function efit = readPoloidalFlux(in, out, efit)

    in.source = 'efitxx';
    % EFIT poloidal flux
    in.hdf5flag = '/output/profiles2D/poloidalFlux';
    efit.polflux    = readMatrixData(in);        

    
    
end



function efit = readMagneticAxis(in, efit)

    in.source = 'efitxx';
    % EFIT magnetic axis
    in.hdf5flag = '/output/globalParameters/magneticAxisR';
    efit.magnax.r     = readScalarData(in);
    in.hdf5flag = '/output/globalParameters/magneticAxisZ';
    efit.magnax.z     = readScalarData(in);
end


function efit = readMagneticBoundary(in, efit)
    
    in.source = 'efitxx';   
    in.hdf5flag = '/output/separatrixGeometry/boundaryClassification';
    efit.boundary     = readScalarData(in);
            
    % EFIT x point
    in.hdf5flag = '/output/separatrixGeometry/xpointCoordsR';
    efit.spx.r     = readVectorData(in);
    in.hdf5flag = '/output/separatrixGeometry/xpointCoordsZ';
    efit.spx.z     = readVectorData(in);
    
    % EFIT limiter
    in.hdf5flag = '/output/separatrixGeometry/limiterCoordsR';
    efit.limiter.r     = readScalarData(in);
    in.hdf5flag = '/output/separatrixGeometry/limiterCoordsZ';
    efit.limiter.z     = readScalarData(in);
end
    
    

function out = normaliseFlux(in, out, efit)
    
    in.source = 'efitxx';
    in.hdf5flag = '/output/singleFluxSurface/poloidalFlux';
    polflux1dEFIT      = readScalarData(in);
    
    in.hdf5flag = '/output/singleFluxSurface/normalizedPoloidalFlux';
    polflux1dfactorEFIT   = readScalarData(in);
    
    
    polfluxMagnAx = interp2(out.flux.r,out.flux.z,efit.polflux,efit.magnax.r,efit.magnax.z);
    
    % Flux normalisation
    % f(x) = (F(x) - F0)/(Fq - F0) * q
    out.flux.normPolFlux = (efit.polflux - polfluxMagnAx) / (polflux1dEFIT - polfluxMagnAx) * polflux1dfactorEFIT;
end


% nt_compass_mx matrix for RENATE 1.1.0
function saveRenateNT(in, out)
    
    for i=1:length(out.nt.density)
        nt_compass_mx(i,:) = [out.nt.psi(i), out.nt.temperature(i), out.nt.density(i), out.nt.zeff, 5.0];
    end
    
    if exist(out.folder.renate) ~= 7
        mkdir(out.folder.renate);
    end
    
    filename = [out.folder.renate,'/nt_',in.tokamak,in.shotNumber,'_',int2str(in.time),'.txt'];    
    disp(['Writing density profile into ',filename]);
    dlmwrite(filename, nt_compass_mx, 'precision','%.4E','delimiter','\t');
end

% flux_compass_mx matrix for RENATE 1.1.0
function saveRenateFlux(in, out)
    for i=1:size(out.flux.r,1)*size(out.flux.r,2)
        flux_compass_mx(i,:) = [out.flux.r(i), out.flux.z(i), out.flux.normPolFlux(i)];
    end
    
    if exist(out.folder.renate) ~= 7
        mkdir(out.folder.renate);
    end
    
    filename = [out.folder.renate,'/flux_',in.tokamak,in.shotNumber,'_',int2str(in.time),'.txt'];
    disp(['Writing density profile into ',filename]);
    dlmwrite(filename, flux_compass_mx, 'precision','%.4E','delimiter','\t');
end


function ts = profileHack(ts,out)
	density = ts.density;
	densityErr = ts.densityErr;
	temperature = ts.temperature;
	temperatureErr = ts.temperatureErr;
	psi_in = ts.psi;
	psi_out=[];
	Te_out=[];
	ne_out=[];
	l = length(out.flux.r);
	l2=100;
	for i = 1:l
		psi_out=[psi_out,psi_in(i)*ones(1,l2)];		
		ne_out=[Te_out,random('norm',density(i),densityErr(i),l2)];
		Te_out=[Te_out,random('norm',temperature(i),temperatureErr(i),l2)];
	end
	ts.psi=psi_out;
	ts.temperature=Te_out;
	ts.density=ne_out;
	

end

function out = fitProfilesNT(ts, out)
    
    ts = profileHack(ts,out);
    
    in = find(ts.psi<max(out.nt.psi_in) & ts.psi > min(out.nt.psi_in));
    
    try
        out.nt.temperature = interp1(ts.psi(in), smooth(ts.psi(in), ts.temperature(in), 0.3, 'rlowess'), out.nt.psi,'linear','extrap');   
        out.nt.density     = interp1(ts.psi(in), smooth(ts.psi(in), ts.density(in), 0.3, 'rlowess'), out.nt.psi,'linear','extrap');     
        disp('  Data smoothed by Fitting Curve Toolbox.') 
    catch e
        try
            out.nt.temperature = interp1(ts.psi(in), smoothdata(ts.psi(in), ts.temperature(in), 0.3, 'rlowess'), out.nt.psi,'linear','extrap');   
            out.nt.density     = interp1(ts.psi(in), smoothdata(ts.psi(in), ts.density(in), 0.3, 'rlowess'), out.nt.psi,'linear','extrap');      
            disp('  Fitting Curve Toolbox is not installed, but Matlab natively supports data smooth.') 
        catch e
            out.nt.temperature = interp1(ts.psi(in), ts.temperature(in), out.nt.psi,'linear','extrap');   
            out.nt.density     = interp1(ts.psi(in), ts.density(in), out.nt.psi,'linear','extrap'); 
            disp('  WARNING: Data not smoothed! Please install Fitting Curve Toolbox! ')
        end        
    end
    
    out.nt.temperature = out.nt.temperature .* (out.nt.temperature > 0);
    
    out.nt.density     = out.nt.density .* (out.nt.density > 0);
end

function plotProfilesNT (in, out, ts);
    l = length(out.nt.psi_in);

    figure
    hold on
    title (['Electron density profile @ ', upper(in.tokamak), ' #', in.shotNumber, ' (', num2str(in.time),' s)'])
    plot(ts.psi,ts.density,'x')
    plot(out.nt.psi,out.nt.density,'m')
    plot(out.nt.psi_in,out.nt.density(1:l),'r','linewidth',2)
    xlabel('\psi')
    ylabel('n_e [10^{19} m^{-3}]')
    savePlot (in, out, 'ne')    
    
    figure
    hold on
    title (['Electron density profile @ ', upper(in.tokamak), ' #', in.shotNumber, ' (', num2str(in.time),' s)'])
    errorbar(ts.psi,ts.density,ts.densityErr)
    plot(out.nt.psi,out.nt.density,'m')
    plot(out.nt.psi_in,out.nt.density(1:l),'r','linewidth',2)
    xlabel('\psi')
    ylabel('n_e [10^{19} m^{-3}]')
    ylim([0,1.1*max(ts.density)])
    savePlot (in, out, 'ne_err')

    figure
    hold on
    title (['Electron temperature profile @ ', upper(in.tokamak), ' #', in.shotNumber, ' (', num2str(in.time),' s)'])
    plot(ts.psi,ts.temperature,'x')
    plot(out.nt.psi,out.nt.temperature,'m')
    plot(out.nt.psi_in,out.nt.temperature(1:l),'r','linewidth',2)
    xlabel('\psi')
    ylabel('T_e [keV]')
    savePlot (in, out, 'Te')
    
    figure
    hold on
    title (['Electron temperature profile @ ', upper(in.tokamak), ' #', in.shotNumber, ' (', num2str(in.time),' s)'])
    errorbar(ts.psi,ts.temperature,ts.temperatureErr)
    plot(out.nt.psi,out.nt.temperature,'m')
    plot(out.nt.psi_in,out.nt.temperature(1:l),'r','linewidth',2)
    xlabel('\psi')
    ylabel('T_e [keV]')
    ylim([0,1.1*max(ts.temperature)])
    savePlot (in, out, 'Te_err')

end

function plotNormFlux (in, out, efit)
    warning('off', 'MATLAB:HandleGraphics:noJVM')

    figure
    hold on
    title(['Normalised flux @ ', upper(in.tokamak), ' #', in.shotNumber, ' (', num2str(in.time),' s)'])
    
    contourStep = 0.005;
    contourMin = 0;
    contourMax = 1;
    clevel = contourMin-contourStep/2:contourStep:contourMax+contourStep/2;
    
    contourf(out.flux.r,out.flux.z,out.flux.normPolFlux,clevel,'color','none')
    contour(out.flux.r,out.flux.z,out.flux.normPolFlux,0:0.1:1,'k')
    caxis(clevel([1 end])-contourStep); % fix up the color so it spans across my desired levels.
    
    if strcmp(efit.boundary, 'xpoint')
        plot(efit.spx.r,efit.spx.z,'wx')
    else        
        plot([min(efit.r),max(efit.r)],[efit.limiter.z,efit.limiter.z],'w')
        plot([efit.limiter.r,efit.limiter.r],[min(efit.z),max(efit.z)],'w')
        plot(efit.limiter.r,efit.limiter.z,'wo')
    end
    
    
    xlabel('R [m]')
    ylabel('z [m]')
       
    colorbar
    axis equal
    
    savePlot(in, out, 'psi')
    
end

function  outputData = readScalarData(in)
    rawData = readDataFile(in);
    outputData = rawData(getTimeIndex(in));
    % disp(['Reading ',in.hdf5flag, '\t 0d(', num2str(size(outputData)), ')'])
end


function  outputData = readVectorData(in)
    rawData = readDataFile(in);
    outputData = rawData(:,getTimeIndex(in));
    % disp(['Reading ',in.hdf5flag, '\t 1d(', num2str(size(outputData)), ')'])
end

function  outputData = readMatrixData(in)
    rawData = readDataFile(in);
    outputData = rawData(:,:,getTimeIndex(in));
    % disp(['Reading ',in.hdf5flag, '\t 2d(', num2str(size(outputData)), ')'])
end

function outputin = readDataFile(in)
    switch in.source
        case 'efitxx'
            filename = [in.folder,'/',in.shotNumber,'/EFITXX/EFITXX.1.h5'];
            outputin = h5read(filename, in.hdf5flag);
        case 'ts'
            filename = [in.folder,'/',in.shotNumber,'/THOMSON/',in.hdf5flag,'.1.h5'];
            outputin = h5read(filename, ['/', in.hdf5flag]);
    end
end

function index = getTimeIndex(in)
    switch in.source
        case 'efitxx'
            index = in.index.efit;
        case 'ts'
            index = in.index.ts;
        otherwise
            disp('ERROR: Corrupt time index')
    end
end


function startclean
    clear all
    clc
    close all
end

% //! read hdf5 data files

function cdb_reader(varargin)
    startclean
    add_csapi_if_its_not_installed

    in.tokamak = 'compass';

    in.plot = true;
    in.renate = true;
    in.folder = '../input/cdb';
    out.folder.renate = '../input/renate110';
    out.folder.grid = '../input/fieldGrid/';
    out.folder.spline = '../input/fieldSpl/';
    out.folder.plot = '../results/';

    %% input
    if nargin == 0
        in.shotNumber = '11774';
        in.time = 1099;
    elseif nargin == 1
        in.shotNumber = '11774';
        in.time = varargin{1};
    else
        in.shotNumber = varargin{1};
        in.time = varargin{2};
    end

    if (nargin >=4 && ~isempty( varargin{4}) )
        in.make_field.e.loaded = varargin{3};
        in.make_field.e.equation = varargin{4};
    else
        in.make_field.e.loaded = 0;
    end

    if (strcmp(in.shotNumber,'test') && nargin >=5 && ~isempty( varargin{5}) )
        in.make_field.b.equation = varargin{5};
        in.make_field.b.loaded = 1;
    else
        in.make_field.b.loaded = 0;
    end

    %% start
    out = getHiddenParameters(out);

    efit = initEmpty;

    if in.make_field.b.loaded
        efit = makeMagneticGrid (in, out, efit);
    else
        in = readTimes(in);
        [efit, out] = readEfitGrid(in, out, efit);
        
        efit = readMagneticAxis(in, efit);
        efit = readMagneticBoundary(in, efit);
        
        efit = readPoloidalFlux(in, out, efit);
        efit = readToroidalFlux(in, out, efit);

        efit = calcMagneticGrid (in, out, efit);
    end

    saveMagneticGrid (in, out, efit);
    saveMagneticSpline (in, out, efit);
    savePsiProfile (in, out, efit);

    if in.make_field.e.loaded
        efit = makeElectricGrid (in, out, efit);
        saveElectricGrid (in, out, efit);
        saveElectricSpline (in, out, efit);
    end
end

function efit = calcMagneticGrid (in, out, efit)
    psi_RZ = efit.polflux';
    R_M = out.flux.r;
    Z_M = out.flux.z;
    dx = 1e-6;
    psi_dR1 = interp2(R_M,Z_M, psi_RZ, R_M-dx/2,Z_M, 'spline');
    psi_dR2 = interp2(R_M,Z_M, psi_RZ, R_M+dx/2,Z_M, 'spline');
    psi_dZ1 = interp2(R_M,Z_M, psi_RZ, R_M,Z_M-dx/2, 'spline');
    psi_dZ2 = interp2(R_M,Z_M, psi_RZ, R_M,Z_M+dx/2, 'spline');

    psi_dR = (psi_dR2-psi_dR1)/dx;
    psi_dZ = (psi_dZ2-psi_dZ1)/dx;

    brad = -psi_dZ./R_M;
    bz   =  psi_dR./R_M;
    efit.brad = brad';
    efit.bz = bz';
end

function efit = makeMagneticGrid (in, out, efit)
    e = lower(in.make_field.b.equation);
    e = regexprep(e,'*','.*');
    e = regexprep(e,'/','./');
    e = regexprep(e,'\^','.\^');

    s = [50,50];
    efit.r = linspace(0.5,1.5,s(1));
    efit.z = linspace(-0.5,0.5,s(2));
    efit.polflux = 'test';
    r = efit.r;
    z = efit.z;

    brad = zeros(s);
    bz = zeros(s);
    btor = zeros(s);

    try
        disp([e,';'])
        eval([e,';'])
    catch
        disp(['ERROR: magnetic_field_value is invalid (',in.make_field.b.equation,')'])
    end

    efit.brad = brad .* ones(s); efit.brad = efit.brad';
    efit.bz   = bz   .* ones(s); efit.bz   = efit.bz';
    efit.btor = btor .* ones(s); efit.btor = efit.btor';
end

function efit = makeElectricGrid (in, out, efit)
    e = lower(in.make_field.e.equation);
    e = regexprep(e,'*','.*');
    e = regexprep(e,'/','./');
    e = regexprep(e,'\^','.\^');

    r = efit.r;
    z = efit.z;
    s = size(efit.brad);

    erad=zeros(s);
    ez=zeros(s);
    etor=zeros(s);

    try
        disp([e,';'])
        eval([e,';'])
    catch
        disp(['ERROR: electric_field_value is invalid (',in.make_field.e.equation,')'])
        %keyboard
    end

    efit.erad = erad .* ones(s); efit.erad = efit.erad';
    efit.ez   = ez   .* ones(s); efit.ez   = efit.ez';
    efit.etor = etor .* ones(s); efit.etor = efit.etor';
end

function saveMagneticGrid (in, out, efit) 
    complist = {'r',    'z',    'brad','bz','btor','polflux', 'psi_n'};
    filelist = {'rcord','zcord','brad','bz','btor','psi2', 'psi_n'};
    saveGrid (in, out, efit, complist, filelist);
end

function saveElectricGrid (in, out, efit) 
    complist = {'erad','ez','etor'};
    filelist = complist;
    saveGrid (in, out, efit, complist, filelist);
end

function saveGrid (in, out, efit, complist, filelist)
    foldername = ([out.folder.grid,'/', in.shotNumber,'_',num2str(in.time),'/']);
    mkdir(foldername)

    for i = 1:length(complist)
        dlmwrite([foldername, filelist{i},'.dat'],efit.(complist{i}), 'precision','%.16e','delimiter','\t');
    end
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
    saveSplineCoeffs (in, out, efit,complist);
end

function savePsiProfile (in, out, efit)
    complist = {'polflux', 'psi_n'};
    saveSplineCoeffs (in, out, efit,complist);
end

function saveElectricSpline (in, out, efit)
    complist = {'erad','ez','etor'};
    saveSplineCoeffs (in, out, efit,complist);
end

function saveSplineCoeffs (in, out, efit,complist)
    for i=1:length(complist)
        comp=complist{i};
        try
            sp = csapi({efit.r,efit.z},efit.(comp));

            i11 = 1:sp.pieces(1);
            i12 = i11+1*sp.pieces(1);
            i13 = i11+2*sp.pieces(1);
            i14 = i11+3*sp.pieces(1);

            i21 = 1:sp.pieces(2);
            i22 = i21+1*sp.pieces(2);
            i23 = i21+2*sp.pieces(2);
            i24 = i21+3*sp.pieces(2);

            c = reshape(sp.coefs,size(sp.coefs,2),[]);
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
            disp(['Spline saved to ',foldername])
        catch
            error(['Error in spline fitting: \n sizeof ',comp,': ',num2str(size(efit.(comp))),' must be ',num2str(length(efit.r)),' x ',num2str(length(efit.z))])
        end
    end
end

function out = getHiddenParameters(out)
    out.nt.zeff = 2.0;
    out.nt.psi  = 0:.01:2;
    out.nt.psi_in  = 0:.01:1.25;
end

function efit = initEmpty
    efit = '';
end

function in = readTimes(in)
    % EFIT
    in.source = 'efitxx';
    in.hdf5flag = '/time';
    in.index.efit=find(readDataFile(in)>=in.time/1000,1);
end

function efit = readToroidalFlux(in, out, efit)
    in.source = 'efitxx';

    %B*Rgeom
    in.hdf5flag = '/output/globalParameters/bvacRgeom';
    efit.bvacrgeom = readScalarData(in);

    in.hdf5flag = '/output/fluxFunctionProfiles/rBphi';
    rbtor = readVectorData(in);
    
    in.hdf5flag = '/output/fluxFunctionProfiles/poloidalFlux';
    polflux = readVectorData(in);

    in.hdf5flag = '/output/fluxFunctionProfiles/normalizedPoloidalFlux';
    normpolflux = readDataFile(in);

    psi_RZ = efit.polflux';
    rbtor = interp1(polflux, rbtor, psi_RZ, 'spline', rbtor(end));
    btor = -rbtor./out.flux.r;
    efit.btor =  btor';

    efit.psi_n = interp1(polflux, normpolflux, efit.polflux, 'spline');
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
    efit.polflux    = readMatrixData(in)';
    
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

function  outputData = readScalarData(in)
    rawData = readDataFile(in);
    outputData = rawData(getTimeIndex(in));
end

function  outputData = readVectorData(in)
    rawData = readDataFile(in);
    outputData = rawData(:,getTimeIndex(in));
end

function  outputData = readMatrixData(in)
    rawData = readDataFile(in);
    outputData = rawData(:,:,getTimeIndex(in));
end

function outputin = readDataFile(in)
    switch in.source
        case 'efitxx'
            filename = [in.folder,'/',in.shotNumber,'/EFITXX/EFITXX.1.h5'];
            outputin = h5read(filename, in.hdf5flag);
    end
end

function index = getTimeIndex(in)
    switch in.source
        case 'efitxx'
            index = in.index.efit;
        otherwise
            disp('ERROR: Corrupt time index')
    end
end

function startclean
    clear all
    clc
    close all
end

function add_csapi_if_its_not_installed
    if ~exist('csapi')
        p = '/home/maradi/public/splines';
        addpath(p)
        disp(['CSAPI loaded from ',p])
    end
end

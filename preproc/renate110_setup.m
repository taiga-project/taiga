% //! read hdf5 data files

function renate110_setup(varargin)
    startclean

    in.tokamak = 'compass';
    in.majorradius = 0.56;

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

    ts = initEmpty;

    if (in.renate && ~in.make_field.b.loaded)
        in = readTimes(in);
        out = efitManager(in, out);
        ts = readThomsonData(in, out, ts);
        %out = fitProfilesNT(ts, out);
        saveRenateFlux(in, out);
        %saveRenateNT (in, out);
        if in.plot
            %plotProfilesNT (in, out, ts);
        end
    end
end

function out = efitManager(in, out)
    in.source = 'efitxx';
    % EFIT major radius (R)
    in.hdf5flag = '/output/profiles2D/r';
    efit.r       = readVectorData(in);
    % EFIT vertical position (Z)
    in.hdf5flag = '/output/profiles2D/z';
    efit.z       = readVectorData(in);
    [out.flux.r,out.flux.z] = meshgrid(efit.r,efit.z);

    % EFIT poloidal flux
    in.hdf5flag = '/output/profiles2D/poloidalFlux';
    polflux2D    = readMatrixData(in)';
        
    in.hdf5flag = '/output/singleFluxSurface/poloidalFlux';
    polflux1dEFIT      = readScalarData(in);

    in.hdf5flag = '/output/singleFluxSurface/normalizedPoloidalFlux';
    polflux1dfactorEFIT   = readScalarData(in);

    % EFIT magnetic axis
    in.hdf5flag = '/output/globalParameters/magneticAxisR';
    magnax_r     = readScalarData(in);
    in.hdf5flag = '/output/globalParameters/magneticAxisZ';
    magnax_z     = readScalarData(in);

    polfluxMagnAx = interp2(out.flux.r,out.flux.z,polflux2D,magnax_r,magnax_z);
    % Flux normalisation
    % f(x) = (F(x) - F0)/(Fq - F0) * q
    out.flux.normPolFlux = (polflux2D - polfluxMagnAx) / (polflux1dEFIT - polfluxMagnAx) * polflux1dfactorEFIT;
end

function savePlot (in, out, plotname)
    mkdir(out.folder.plot)
    foldername = [out.folder.plot,'/', in.shotNumber,'_',num2str(in.time),'/'];
    mkdir(foldername)
    saveas(gcf,[foldername, plotname, '.pdf'])
    close
end

function out = getHiddenParameters(out)
    out.nt.zeff = 2.0;
    out.nt.psi  = 0:.01:2;
    out.nt.psi_in  = 0:.01:1.25;
end

function ts= initEmpty
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
    in.index.ts = find(readDataFile(in)>=in.time,1);
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

    notnanindex       = (~isnan(ts.densityRaw));
    ts.density        = [ts.densityRaw(notnanindex)];        % 1e19 m-3
    ts.densityErr     = [ts.densityRawErr(notnanindex)];     % 1e19 m-3
    ts.temperature    = [ts.temperatureRaw(notnanindex)];    % keV
    ts.temperatureErr = [ts.temperatureRawErr(notnanindex)]; % keV
    ts.z              = [ts.zRaw(notnanindex)];
    ts.r              = ones(size(ts.z))*in.majorradius;
    ts.psi            = interp2(out.flux.r,out.flux.z,out.flux.normPolFlux,ts.r,ts.z);
end

% nt_compass_mx matrix for RENATE 1.1.0
function saveRenateNT(in, out)
    for i = 1:length(out.nt.density)
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
    for i = 1:size(out.flux.r,1)*size(out.flux.r,2)
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

function out = fitProfilesNT_old(ts, out)
    try
        ts = profileHack(ts,out);
    catch e
    end

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

function out = fitProfilesNT(ts, out)
    try
        stefanikova_l_mode = @(a,r) a(1)*exp(-(sqrt(r)/a(2)).^a(3));
        T_init = [800,0.69,3];
        T_fit = nlinfit(ts.psi,ts.temperature,stefanikova_l_mode,T_init);
        out.nt.temperature = stefanikova_l_mode(T_fit,ts.psi);

        n_init = [4e19,0.69,4.6];
        n_fit = nlinfit(ts.psi,ts.density,stefanikova_l_mode,n_init);
        out.nt.density = stefanikova_l_mode(n_fit,ts.psi);
    end
end

function plotProfilesNT (in, out, ts)
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

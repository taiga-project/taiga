% //! read hdf5 data files

function cdb_reader(varargin)
    startclean

    in.tokamak = 'compass';
    in.majorradius=0.56;  
    
    in.plot = false;
    
    in.folder = '../input/renate110';
    out.folder = '../dataio/data';
    

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
        in.ionNumber = varargin{3};
    else
        in.ionNumber = 100000;
    end
         
    out = readRenateProfile(in);
    
    out = getIonPosition(out);
            
    %saveInputIons(in, out);  
    %saveRenateProfile(in, out);  
    

    if in.plot
        plotInputIons (in, out);    
    end

end

function out = readRenateProfile(in)
    renate_raw = load([in.folder,'/pop_compass',in.shotNumber,'_',in.time,'.txt']);
    out.pos = renate_raw(1,:);
    out.rate = renate_raw(3,:);
end

function out = getIonPosition(out)
    

end

function saveRenateProfile (in, out)    
    
    foldername = ([out.folder.grid,'/', in.shotNumber,'_',num2str(in.time),'/']);
    mkdir(foldername)
    
    filename = [out.folder,'/pos.dat'];   
    disp(['Writing ionisation position into ',filename]);     
    dlmwrite(filename, out.pos, 'precision','%.4E','delimiter','\t');
    
    filename = [out.folder,'/prof.dat'];   
    disp(['Writing ionisation profile into ',filename]);     
    dlmwrite(filename, out.rate, 'precision','%.4E','delimiter','\t');
    
end

function saveInputIons (in, out)    
    
    foldername = ([out.folder.grid,'/', in.shotNumber,'_',num2str(in.time),'/']);
    mkdir(foldername)
    
    filename = [out.folder,'/rad.dat'];   
    disp(['Writing ion coordinates into ',filename]);     
    dlmwrite(filename, out.ions, 'precision','%.4E','delimiter','\t');
    
end




function plotRenateProfile (in, out)

    figure
    hold on
    title(['Ionisation profile @ ', upper(in.tokamak), ' #', in.shotNumber, ' (', num2str(in.time),' s)'])


end

function startclean
    clear all
    clc
    close all
end

% //! read hdf5 data files


clear all
clc
close

shotNumber = '11774'; 		% should be a STRING! not a number !!!
dataNames = {'R_geom_axis','B_vac_R_geom','psi_RZ','R','Z'};
folderName = '../input/shot';

%t_i = 200;
timestamp=1099

dx = 1e-6;



dataStruct=cell(5,1);
timeStruct=cell(5,1);
for file_i=1:5

	dataName = dataNames{file_i}

	filename = [folderName,'/',shotNumber,'/',dataName,'-',shotNumber,'.h5']

	%unix(['ls -l ',filename])

	h5info(filename)
	h5disp(filename)

	data1 = h5read(filename,'/values');
	time1 = h5read(filename,'/time');
	%dataStruct{file_i}
	dataStruct{file_i} = data1;
	timeStruct{file_i} = time1;
	size(data1)
end

s2 = size(data1);
Nt = s2(end);

% naming data
R_geom_D  = dataStruct{1}; % vector
BR_geom_D = dataStruct{2}; %
psi_RZ_D  = dataStruct{3};
R_D       = dataStruct{4};
Z_D       = dataStruct{5};


% get t_i timeset from dataset
t_i=find(timeStruct{1}>=timestamp,1)


R_geom    = R_geom_D(t_i);
BR_geom   = BR_geom_D(t_i);
r_v       = R_D(:,t_i);
R_size    = length(r_v);
z_v       = Z_D(:,t_i);
Z_size    = length(z_v);
psi_RZ    = psi_RZ_D(:,:,t_i);



% toroidal component of magnetic field
btor_v = -BR_geom*R_geom./(r_v) ; 
z_o = ones(1,Z_size);
r_o = ones(1,R_size);
rz_o = ones(R_size,Z_size);
btor = kron(btor_v,z_o);

%plot(r_v,btor_v)


psi_RZ	 = psi_RZ';

[R_M,Z_M]=meshgrid(r_v,z_v);
%R_M
%Z_M

%contour(R_M,Z_M,psi_RZ',30)
%axis equal	


psi_dR1 = interp2(R_M,Z_M,psi_RZ,R_M-dx/2,Z_M,'spline');
psi_dR2 = interp2(R_M,Z_M,psi_RZ,R_M+dx/2,Z_M,'spline');
psi_dZ1 = interp2(R_M,Z_M,psi_RZ,R_M,Z_M-dx/2,'spline');
psi_dZ2 = interp2(R_M,Z_M,psi_RZ,R_M,Z_M+dx/2,'spline');

psi_dR = (psi_dR2-psi_dR1)/dx;
psi_dZ = (psi_dZ2-psi_dZ1)/dx;


brad = -psi_dZ./R_M;%*pi;
bz   =  psi_dR./R_M;%*pi;	
brad=brad';
bz=bz';

%R_M
%Z_M

if false

cmin=-0.2;
cmax=0.2;

figure
subplot(1,2,1)
contour(R_M',Z_M',brad,linspace(cmin,cmax,20))
caxis([ cmin cmax ])
colorbar


subplot(1,2,2)
contour(R_M',Z_M',bz,linspace(cmin,cmax,20))
caxis([ cmin cmax ])
colorbar
error
end

if false
contour(R_M,Z_M,psi_RZ')
axis equal
xlim([0.3 0.8])
saveas(gcf,'field_psi.svg')

end

if false
cd '../input/fieldGrid'
mkdir([shotNumber,'_',num2str(timestamp),'_old'])
cd([shotNumber,'_',num2str(timestamp),'_old'])
save('rcord.dat','r_v','-ascii','-double','-tabs')
save('zcord.dat','z_v','-ascii','-double','-tabs')
save('brad.dat','brad','-ascii','-double','-tabs')
save('bz.dat','bz','-ascii','-double','-tabs')
save('psi2.dat','psi_RZ','-ascii','-double','-tabs')
save('btor.dat','btor','-ascii','-double','-tabs')
end


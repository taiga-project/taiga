% //! read magnetic field and make spline coefficients


function fieldSpl

    close all
    clear all
    clc
    
    shotnumber = '11347';
        
    inDir  = ['../input/fieldGrid/',shotnumber];
    outDir = ['../input/fieldSpl/',shotnumber];
	filenames = {'brad','bz','btor'}; 
	for i=1:3
		filename=filenames{i};
		if true
			[XS,YS,ZS]=peaks(100);

			x = load([inDir,'/rcord.dat'],'-ascii');
			y = load([inDir,'/zcord.dat'],'-ascii');
	
			Z = load([inDir,'/',filename,'.dat'],'-ascii');
			[X,Y]=meshgrid(x,y);
		
		end
		
		main(filename,inDir,outDir)
    
    end
      
end


function main(filename,inDir,outDir)
    sample=rand(1,2);
    sample(1)=sample(1)*.45+0.2906;
    sample(2)=sample(2)*1.14-0.665;


    x = load([inDir,'/rcord.dat'],'-ascii');
    y = load([inDir,'/zcord.dat'],'-ascii');
	
    Z = load([inDir,'/',filename,'.dat'],'-ascii');

    sp = csapi({x,y},Z);

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

	mkdir(outDir)
    save([outDir,'/r.spline'],'-ascii','-tabs','-double','b')
    save([outDir,'/z.spline'],'-ascii','-tabs','-double','b2')
    disp('READY')

    c11 = c(i11,i21);    c12 = c(i11,i22);    c13 = c(i11,i23);    c14 = c(i11,i24);
    c21 = c(i12,i21);    c22 = c(i12,i22);    c23 = c(i12,i23);    c24 = c(i12,i24);
    c31 = c(i13,i21);    c32 = c(i13,i22);    c33 = c(i13,i23);    c34 = c(i13,i24);
    c41 = c(i14,i21);    c42 = c(i14,i22);    c43 = c(i14,i23);    c44 = c(i14,i24);
    
    for fi1 = 1:4
    	for fi2 = 1:4
    	
    		locVar = genvarname(['c',num2str(fi1),num2str(fi2)]);
		    save([outDir,'/',filename,'.spl',num2str(fi1),num2str(fi2)],'-ascii','-tabs','-double',locVar)
        end
    end

    if false
	    r=rand(1,length(b(2:end)));
	    dx = ((b(2:end)-b(1:end-1)).*r)';

	    x2 = x(1:end-1)+dx';


	    r=rand(1,length(b2(2:end)));
	    dy = ((b2(2:end)-b2(1:end-1)).*r)';

	    y2 = y(1:end-1)+dy';


	    [dX,dY]=meshgrid(dx,dy);

	    [X2,Y2]=meshgrid(x2,y2);


	    dX=dX';
	    dY=dY';

	   Z20= c11.*dX.^3.*dY.^3 + c12.*dX.^3.*dY.^2 + c13.*dX.^3.*dY + c14.*dX.^3 + ...
		c21.*dX.^2.*dY.^3 + c22.*dX.^3.*dY.^2 + c23.*dX.^2.*dY + c24.*dX.^2 + ...
		c31.*dX   .*dY.^3 + c32.*dX   .*dY.^2 + c33.*dX   .*dY + c34.*dX    + ...
		c41       .*dY.^3 + c42       .*dY.^2 + c43       .*dY + c44;

	    Z2= c11.*dX.^3.*dY.^3 + c21.*dX.^3.*dY.^2 + c31.*dX.^3.*dY + c41.*dX.^3 + ...
		c12.*dX.^2.*dY.^3 + c22.*dX.^3.*dY.^2 + c32.*dX.^2.*dY + c42.*dX.^2 + ...
		c13.*dX   .*dY.^3 + c23.*dX   .*dY.^2 + c33.*dX   .*dY + c43.*dX    + ...
		c14       .*dY.^3 + c24       .*dY.^2 + c34       .*dY + c44;



	    bs1 = find((b <= sample(1)),1,'last');
	    bs2 = find((b2 <= sample(2)),1,'last');
	    
	    if bs1==length(b)
	       bs1=bs1-1; 
	    end

	    if bs2==length(b2)
	       bs2=bs2-1; 
	    end

	    dsx = sample(1)-b(bs1);
	    dsy = sample(2)-b2(bs2);


	 try
	% EZ A JÓ
	    sample2(3) =c11(bs1,bs2)*dsx^3*dsy^3 + c12(bs1,bs2)*dsx^3*dsy^2 + c13(bs1,bs2)*dsx^3*dsy + c14(bs1,bs2)*dsx^3 + ...
		        c21(bs1,bs2)*dsx^2*dsy^3 + c22(bs1,bs2)*dsx^3*dsy^2 + c23(bs1,bs2)*dsx^2*dsy + c24(bs1,bs2)*dsx^2 + ...
		        c31(bs1,bs2)*dsx  *dsy^3 + c32(bs1,bs2)*dsx  *dsy^2 + c33(bs1,bs2)*dsx  *dsy + c34(bs1,bs2)*dsx    + ...
		        c41(bs1,bs2)      *dsy^3 + c42(bs1,bs2)      *dsy^2 + c43(bs1,bs2)      *dsy + c44(bs1,bs2);
	 end
 
    
   end
    
    
end


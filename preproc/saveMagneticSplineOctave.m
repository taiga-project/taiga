function saveMagneticSplineOctave (comp,folder_grid,folder_spl)

efit_r = load([folder_grid, 'rcord.dat'])
efit_z = load([folder_grid, 'rcord.dat'])
efit_comp = load([folder_grid, comp, '.dat'])
      
        sp = csapi({efit_r,efit_z},efit_comp);

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

        save([folder_spl,'/r.spline'],'-ascii','-tabs','-double','b')
        save([folder_spl,'/z.spline'],'-ascii','-tabs','-double','b2')
        disp('Spline calculated by OCTAVE')        
        
        c11 = c(i11,i21);    c12 = c(i11,i22);    c13 = c(i11,i23);    c14 = c(i11,i24);
        c21 = c(i12,i21);    c22 = c(i12,i22);    c23 = c(i12,i23);    c24 = c(i12,i24);
        c31 = c(i13,i21);    c32 = c(i13,i22);    c33 = c(i13,i23);    c34 = c(i13,i24);
        c41 = c(i14,i21);    c42 = c(i14,i22);    c43 = c(i14,i23);    c44 = c(i14,i24);
        
        for fi1 = 1:4
        	for fi2 = 1:4        	
		        save([folder_spl,'/',comp,'.spl',num2str(fi1),num2str(fi2)],'-ascii','-tabs','-double',['c',num2str(fi1),num2str(fi2)])
            end
        end
        disp(comp)
        disp(['Spline saved to ',folder_spl])
   
end

clear all
    mainfolder = '../results';
    shotnumber = 'test_1';    runnumber = '120';
    shotnumber = 'test_2';    runnumber = '128';
    shotnumber = 'test_5';    runnumber = '129';
    shotnumber = 'test_6';    runnumber = '130';
    
    
    load([mainfolder,'/',shotnumber,'/',runnumber,'/t_rad.dat']);
    load([mainfolder,'/',shotnumber,'/',runnumber,'/t_z.dat']);
    load([mainfolder,'/',shotnumber,'/',runnumber,'/t_tor.dat']);
    

    t_rad = t_rad';
    t_z = t_z';
    t_tor = t_tor';
    
if false
    t_rad = t_rad([1,1],:);
    t_z = t_z([1,1],:);
    t_tor = t_tor([1,1],:);
end



    L1=t_rad(:,2:end)>t_rad(:,1);
    L2=t_rad(:,1:end-1)<t_rad(:,1);
    [I1,I2] = find(L1.*L2);
    I = find(L1.*L2);
    
    zi = zeros(size(I));
ri = zi;
    ind = cell(max(I1),1);

for i = 1:length(I)    
    ind{I1(i)}=[ind{I1(i)},i];  
    ri(i) = t_rad(I1(i),1);   
end
   
s1 = size(t_rad,1) 
    
zi = t_z(I) + (t_z(I+s1) - t_z(I)) ./ (t_rad(I+s1) - t_rad(I)) .* (t_rad(I+s1) - ri) ;



    L1r=t_rad(:,2:end)<t_rad(:,1);
    L2r=t_rad(:,1:end-1)>t_rad(:,1);
    [I1r,I2r] = find(L1r.*L2r);
    Ir = find(L1r.*L2r);
    
    zir = zeros(size(Ir)); rir=zir;
    indr = cell(max(I1r),1);

for i = 1:length(Ir)    
    indr{I1r(i)}=[indr{I1r(i)},i];    
    rir(i) = t_rad(I1r(i),1);   
end
    
zir = t_z(Ir) + (t_z(Ir+s1) - t_z(Ir)) ./ (t_rad(Ir+s1) - t_rad(Ir)).* (t_rad(Ir+s1) - rir) ;




    av = nan(max(I1),1);
    dif= nan(max(I1),1);
    avs = nan(max(I1),1);
    difs= nan(max(I1),1);
for i=1:length(indr)
    r = zir(indr{i});
    l = zi(ind{i});
    r_i = I2r(indr{i});
    l_i = I2(ind{i});
    
    m = min([length(l),length(r)]);
    ix = [l(1:m),r(1:m)]';
    ix = ix(:);


av(i) = mean([diff(r)./diff(r_i);diff(l)./diff(l_i)]); 
avs(i) = std([(r(2:end)-r(1:end-1))./(r_i(2:end)-r_i(1:end-1));(l(2:end)-l(1:end-1))./(l_i(2:end)-l_i(1:end-1))]); 
dif(i) = mean((l(1:m)-r(1:m)));
difs(i) = std((l(1:m)-r(1:m))); 
end


disp(['Larmor radius: ',num2str(mean(dif*100),'%.4f'), ' +/- ',num2str(std(difs*100),'%.4f'),' cm'])

disp(['Av. drift: ',num2str(mean(av*1e6),'%.4f'), ' +/- ',num2str(std(avs*1e6),'%.4f'),' um/step'])


B = 1;
dB = 0.002  ;
dt = 1e-9;
E = 10

xgradB = mean(60 * 1000 ./ (B +dB*t_rad(:,1)).^ 2 *dB *dt);
xE = mean(60 * 1000 /B*E *dt);
disp(['expected gradB drift: ',num2str(xgradB*1e6,'%.4f'),' um/step'])
disp(['expected E field drift: ',num2str(xE*1e6,'%.4f'),' um/step'])

if true
    close all
    figure
    plot(t_rad(1,:),t_z(1,:))
    hold on 
    plot(t_rad(1,1)*ones(size(ind{1})),zi(ind{1}),'rx')
%plot(t_rad(I(indr{1})),t_z(I(indr{1})),'rx')
 %   plot(t_rad(1,1)*ones(size(indr{1})),zir(indr{1}),'gx')
end


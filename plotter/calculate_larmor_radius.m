clear all
    mainfolder = '../results';
    shotnumber = 'test_1';    runnumber = '120';
    shotnumber = 'test_2';    runnumber = '128';
    shotnumber = 'test_5';    runnumber = '129';
    
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
    ind = cell(max(I1),1);

for i = 1:length(I)    
    ind{I1(i)}=[ind{I1(i)},i];    
end
    
    
zi = t_z(I) + (t_z(I+1) - t_z(I)) ./ (t_rad(I+1) - t_rad(I)) ;

if true
    close all
    figure
    plot(t_rad(1,:),t_z(1,:))
    hold on 
    plot(t_rad(1,1)*ones(size(ind{1})),zi(ind{1}),'rx')
end



    L1r=t_rad(:,2:end)<t_rad(:,1);
    L2r=t_rad(:,1:end-1)>t_rad(:,1);
    [I1r,I2r] = find(L1r.*L2r);
    Ir = find(L1r.*L2r);
    
    zir = zeros(size(Ir));
    indr = cell(max(I1r),1);

for i = 1:length(Ir)    
    indr{I1r(i)}=[ind{I1r(i)},i];    
end
    
zir = t_z(Ir) + (t_z(Ir+1) - t_z(Ir)) ./ (t_rad(Ir+1) - t_rad(Ir)) ;




    av = nan(max(I1),1);
    dif= nan(max(I1),1);
    avs = nan(max(I1),1);
    difs= nan(max(I1),1);
for i=1:length(indr)
    r = zir(indr{i});
    l = zi(ind{i});
    m = min([length(l),length(r)]);
    ix = [l(1:m),r(1:m)]';
    ix = ix(:);

av(i) = mean([r(2:end)-r(1:end-1);l(2:end)-l(1:end-1)]); %mean(ix(2:end)-ix(1:end-1));
avs(i) = std([r(2:end)-r(1:end-1);l(2:end)-l(1:end-1)]);
dif(i) = mean(l(1:m)-r(1:m));
difs(i) = std(l(1:m)-r(1:m));    
end


disp(['Larmor radius: ',num2str(mean(dif*100),'%.4f'), ' +/- ',num2str(std(difs*100),'%.4f'),' cm'])

disp(['drift: ',num2str(mean(av*100),'%.4f'), ' +/- ',num2str(std(avs*100),'%.4f'),' cm'])


B = 1;
dB = 0.01;
dt = 1e-9;
t_step=2000;

xgradB = 60 * 1000 / B ^ 2 *dB *dt*t_step;
disp(['expected gradB drift: ',num2str(xgradB*100,'%.4f'),' cm'])



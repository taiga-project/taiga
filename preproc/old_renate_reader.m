nt_compass30866_106=load('../../../renate/Tags/1.0-stable/data/nt_compass30866_106.txt');
phi=(nt_compass30866_106(:,1));
te=(nt_compass30866_106(:,2));
ne=(nt_compass30866_106(:,3));

out_indices = find(phi>1&te>0&ne>0);

phi_out=phi(out_indices);
te_out=te(out_indices);
ne_out=ne(out_indices);

close all
figure
plot(phi_out,te_out)

figure
plot(phi_out,ne_out)

figure
loglog(phi_out,te_out)

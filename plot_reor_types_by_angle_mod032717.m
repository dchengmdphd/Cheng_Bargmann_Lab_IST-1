%plots different angles of reorientations on same plot
%make sure a linkedtracks is loaded...


reor_type = 'all';

[reor_0_44_minus1,~] = events_per_animal_twindows_byangle_mod032717(Tracks,1,reor_type,3, 0, 44.9);
avg_reor_0_44_minus1 = nanmean(reor_0_44_minus1,1);

[reor_45_89_minus1,~] = events_per_animal_twindows_byangle_mod032717(Tracks,1,reor_type,3, 45, 89);
avg_reor_45_89_minus1 = nanmean(reor_45_89_minus1,1);

[reor_90_134_minus1,~] = events_per_animal_twindows_byangle_mod032717(Tracks,1,reor_type,3, 90, 134);
avg_reor_90_134_minus1 = nanmean(reor_90_134_minus1,1);

[reor_135_180_minus1,~] = events_per_animal_twindows_byangle_mod032717(Tracks,1,reor_type,3, 135, 180);
avg_reor_135_180_minus1 = nanmean(reor_135_180_minus1,1);

figure();
hold on;
plot(avg_reor_0_44_minus1,'r-');
plot(avg_reor_45_89_minus1,'g-');
plot(avg_reor_90_134_minus1,'b-');
plot(avg_reor_135_180_minus1,'m-');
title('Dus track example: allnonUps angle at index ONE BEFORE reori');
xlabel('time (min)');
ylabel('frequency')';

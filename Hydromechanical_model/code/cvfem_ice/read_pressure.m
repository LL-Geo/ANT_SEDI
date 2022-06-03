sigx_t_top=sigx_t(2811:3091,:);

sigy_t_top=sigy_t(2811:3091,:);

sigz_t_top=sigz_t(2811:3091,:);

sigxy_t_top=sigxy_t(2811:3091,:);

sigkk_t_top=sig_kk_t(2811:3091,:);


figure()
subplot(2,2,1);
imagesc(sigx_t_top)
title('sigx_t_top')
colormap('jet')
subplot(2,2,2);
imagesc(sigy_t_top)
colormap('jet')
title('sigy_t_top')

subplot(2,2,3);
imagesc(sigz_t_top)
colormap('jet')
title('sigz_t_top')


subplot(2,2,4);
imagesc(sigkk_t_top)
colormap('jet')
title('sigkk_t_top')

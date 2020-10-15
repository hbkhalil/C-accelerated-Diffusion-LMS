%% compilation file 

if ismac || isunix
    mex -O CFLAGS="\$CFLAGS -std=c99" ATC_mex.c
    mex -O CFLAGS="\$CFLAGS -std=c99" compressed_diffusion.c
    mex -O CFLAGS="\$CFLAGS -std=c99" doubly_compressed_diffusion.c
    mex -O CFLAGS="\$CFLAGS -std=c99" ATC_partial_model.c
    mex -O CFLAGS="\$CFLAGS -std=c99" ATC_RCD.c
    mex -O CFLAGS="\$CFLAGS -std=c99" compressive_diff_ATC.c
    mex -O CFLAGS="\$CFLAGS -std=c99" ATC_En.c
    mex -O CFLAGS="\$CFLAGS -std=c99" compressed_diffusion_En.c
    mex -O CFLAGS="\$CFLAGS -std=c99" doubly_compressed_diffusion_En.c
    mex -O CFLAGS="\$CFLAGS -std=c99" ATC_RMT.c
    mex -O CFLAGS="\$CFLAGS -std=c99" ATC_RMT_theo.c
elseif ispc
    mex  ATC_mex.c
    mex  compressed_diffusion.c
    mex  doubly_compressed_diffusion.c
    mex  ATC_partial_model.c
    mex  ATC_RCD.c
    mex  compressive_diff_ATC.c
else
    disp('platform not suported')
end



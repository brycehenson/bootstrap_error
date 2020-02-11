function fit_params=test_gauss_fit_bootstrap_error_fit_wrapper(x)
    verbose=5;
    shdat=smooth_hist(x);
    norm_count_rate=shdat.count_rate.smooth/numel(x);
    model_fun=@(c,x) c(1)+c(2).*exp(-(1/2).*((x - c(3))./c(4)).^2);
    mdl=fitnlm(shdat.bin.centers,norm_count_rate,model_fun,[1,1,1,1]);
    fit_params=mdl.Coefficients.Estimate;
    if verbose>3
        stfig('gauss fit diagnostics');
        clf
        plot(shdat.bin.centers,shdat.count_rate.smooth)
        hold on
        predict_x=col_vec(linspace(shdat.bin.centers(1),shdat.bin.centers(end),1e3))
        predict_out=predict(mdl,predict_x);
        plot(shdat.bin.centers,shdat.count_rate.smooth)
    end
end
function [img_rec,rmse] = nlcg_new(kdata,it_para, w_para, f_para, tv_para)
%-----------------------------------------------------------------------
%
% 
%-------------------------------------------------------------------------

    img_rec = 0.*kdata;
    rmse = [];

    % iteration parameters
    maxit = it_para.maxiter;
    gradtol = it_para.gradtol;
    alpha = it_para.alpha;
    beta = it_para.beta;
    t_maxiter = it_para.tmax;

    g = grad(img_rec,kdata,w_para,f_para,tv_para);
    dm = -g;
    grd = g(:)'*g(:);
    rmse = [rmse grd];
    fprintf("initial residual = %.2f\n",grd);
    
    for i = 1:maxit
        
        % stop condition
        rmse = [rmse g(:)'*g(:)/grd];
        if rmse(i+1) < gradtol
            break;
        end
        fprintf('%d-th iteration // gradient = %.6f\n', i,rmse(i+1));
    
        % backtraking line-search
        t = 1;
        t_iter = 0;
        f_obj = obj(img_rec+t*dm,kdata,w_para,f_para,tv_para);
        f_line = obj(img_rec,kdata,w_para,f_para,tv_para) - alpha*t*abs(g(:)'*dm(:));
        while ((f_obj > f_line) && (t_iter < t_maxiter))
            t = beta*t;
            f_obj = obj(img_rec+t*dm,kdata,w_para,f_para,tv_para);
            f_line = obj(img_rec,kdata,w_para,f_para,tv_para) - alpha*t*abs(g(:)'*dm(:));
            t_iter = t_iter + 1;
        end
        
        img_rec = img_rec + t*dm;
        g_nxt = grad(img_rec,kdata,w_para,f_para,tv_para);
        gamma = (g_nxt(:)'*g_nxt(:)) / (g(:)'*g(:));
        dm = -g_nxt + gamma*dm;
        
        g = g_nxt;
    end

end


function res = obj(img,kdata,w_para,f_para,tv_para)
    % data consistency term
    FT = f_para.operator;
    tmp = FT*img - kdata;
    tmp_f = tmp.*conj(tmp);
%     tmp_f = fft2c(img).*f_para.mask - kdata;
%     tmp_f = tmp_f.*conj(tmp_f);
    
    % sparsity term
    tmp_w = 0;
    W = w_para.operator;
    if w_para.weight ~= 0
        wav_img = W*img;
        wav_img = sqrt(wav_img.*conj(wav_img)+w_para.smoothingfct);
        tmp_w = wavshow(abs(wav_img),w_para.idx,w_para.filter);
%         [wav_img,widx] = wav2c(img,w_para.lvl,w_para.filter);
%         wav_img = sqrt(wav_img.*conj(wav_img)+w_para.smoothingfct);
%         tmp_w = wavshow(abs(wav_img),widx,w_para.filter);
    end
    
    % TV term
    tmp_tv = 0;
    TV = tv_para.operator;
    if tv_para.weight ~= 0
        tv_img = TV*img;
        tmp_tv = sqrt(tv_img.*conj(tv_img)+tv_para.smoothingfct);
%         tv_img = D(img);
%         tmp_tv = sqrt(tv_img.*conj(tv_img)+tv_para.smoothingfct);
%         tmp_tv = sqrt(tv_img.*conj(tv_img));
    end
    
    res = sum(tmp_f(:)) + w_para.weight*sum(tmp_w(:)) + tv_para.weight*sum(tmp_tv(:));
end

function df = grad(img,kdata,w_para,f_para,tv_para)
    % data consistency term
    gradF = gObj(img,kdata,f_para);
    
    % sparsity term
    gradW = 0;
    if w_para.weight ~= 0
        gradW = gWav(img,w_para);
    end
    
    % TV term
    gradTV = 0;
    if tv_para.weight ~= 0
        gradTV = gTV(img,tv_para);
    end
    
    df = gradF + w_para.weight.*gradW + tv_para.weight.*gradTV;
end

function res = gObj(img,kdata,f_para)
% 2*F'(F(m)-y)
    FT = f_para.operator;
    tmp = FT*img - kdata;
    res = FT'*tmp;
%     tmp = fft2c(img).*f_para.mask - kdata;
%     res = ifft2c(tmp.*f_para.mask);
    res = 2*res;
end

function res = gWav(img,w_para)
% Phi'*W-1*Phi(m)
    W = w_para.operator;
    
    wav_img = W*img;
    tmp = wav_img./sqrt(wav_img .* conj(wav_img) + w_para.smoothingfct);
    res = W'*tmp;
    
%     [wav_img,widx] = wav2c(img,w_para.lvl,w_para.filter);
%     tempDeno = sqrt(wav_img .* conj(wav_img) + w_para.smoothingfct);
%     derImg = wav_img./tempDeno;
%     res = iwav2c(derImg,widx,w_para.filter);
end

function res = gTV(img,tv_para)
    TV = tv_para.operator;
    tmp = TV*img;
    tmp = tmp./sqrt(tmp.*conj(tmp) + tv_para.smoothingfct);
    res = TV'*tmp;
    
%     Dx = D(img);
%     G = Dx./sqrt(Dx.*conj(Dx) + tv_para.smoothingfct);
%     res = invD(G,size(img));
end


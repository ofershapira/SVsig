function [mix_model,opt_alpha] = mix_model_param( mfull, model1, model2, events, bins, CHR )
% finds the alpha's of the mix model pij = alpha*model1 + (1-alpha)*
% alpha's can be determined by a subset of tiles
% the best model is found by minimizing the BIC

annot_tiles=tiles_annot('length',events,bins,CHR);


% setting up some needed variables
num_param=length(annot_tiles(1,1,:));
alpha=0.5*ones(num_param,1);

if issymmetric(model1),
    model1=triu(model1);
    model1(eye(size(model1))==1)=diag(model1)/2;
end

if issymmetric(model2),
    model2=triu(model2);
    model2(eye(size(model2))==1)=diag(model2)/2;
end

if issymmetric(mfull),
    mfull=triu(mfull);
    mfull(eye(size(mfull))==1)=diag(mfull)/2;
end

log_fac(1)=0;
for c1=1:max(mfull(:));
    log_fac(c1+1)=sum(log(1:c1));
end

nume=sum(mfull(:));
alpha=rand(num_param,1);

% optimize over model
options = optimoptions('fmincon','DiffMinChange',1e-6,'TolFun',1e-1,'TolX',1e-6);
[opt_alpha, f_bic] = fmincon(@mix_optim_fun,alpha,[],[],[],[],zeros(num_param,1),ones(num_param,1),[],options);
mix_model=zeros(size(mfull));
for c1=1:num_param,
    mix_model(annot_tiles(:,:,c1)) =  mix_model(annot_tiles(:,:,c1)) + opt_alpha(c1)*model1(annot_tiles(:,:,c1))+(1-opt_alpha(c1))*model2(annot_tiles(:,:,c1));
end
mix_model=mix_model/sum(mix_model(:));

    function BIC = mix_optim_fun(alpha)

        % the mix model probability function
        mix_model=zeros(size(mfull));
        for c1=1:num_param,
            mix_model(annot_tiles(:,:,c1)) =  mix_model(annot_tiles(:,:,c1)) + alpha(c1)*model1(annot_tiles(:,:,c1))+(1-alpha(c1))*model2(annot_tiles(:,:,c1));
        end
        mix_model=mix_model/sum(mix_model(:));
        
        % the log likelihood function
        nnz_idc=mix_model>0&mfull>0;
%        nnz_idc=mix_model>0;
        sLij = sum(sum(mfull(nnz_idc).*log(nume*mix_model(nnz_idc))-nume*mix_model(nnz_idc)-log_fac(mfull(nnz_idc)+1)'));

        % the BIC value
        BIC = -2*sLij+log(nume)*num_param

    end
        
end

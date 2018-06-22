%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% do_irspca_main(EEG) for demonstartion
%
% Deaprtment of Brain and Cognitive Engineering, Korea University 
% Brain Signal Processing Laboraty,BSPL
%
% updated 07/25/2014
%
% Any suggestions or errors, please contact us, hyunchul_kim@korea.ac.kr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [rXbar] = do_irspca_main(sig,dim,depth,loc_vec_ui32,singpdB,thkval,stop_var,chnnel_info,chnnel_idx,fs,fpt,flg_verbose,max_depth,max_nkurt)

FOI.gamma=[20; fs/2]; window=1;

x=sig;

nsmp=length(x)-dim+1;
tdim=length(x);

if flg_verbose == 1;
disp('data matrix is being generated');
end

%% Segmentation
X = flipud(toeplitz(x(dim:-1:1),x(dim:end)));
x=[];

%% Reconstruction
if flg_verbose ==1;
disp('rsPCA...............................');
end

rmX = X- mean(X,2)*ones(1,size(X,2));
[V D]=eig((rmX*rmX')./length(X));

[eigenval,index] = sort(diag(D));

index=rot90(rot90(index)); Dval=rot90(rot90(eigenval))';
D = Dval; V = V(:,index); normalD=D./sum(D);

coeff = (V'*X)'; V_t = V';

rmX =[]; X =[];
reconXbar = zeros(dim,tdim);
Xbar = zeros(nsmp,dim);

if flg_verbose ==1;
tic
end
parfor ii=1:dim
    reconXbar(ii,:) = subfunc_recon(loc_vec_ui32, coeff(:,ii), V_t(ii,:));
end

ids = 1:dim-1;
reconXbar(:,ids) = reconXbar(:,ids)./(ones(dim,1)*ids);

ids = dim:tdim-dim;
reconXbar(:,ids) = reconXbar(:,ids)/dim;

ids = tdim-dim+1:tdim;
reconXbar(:,ids) = reconXbar(:,ids)./(ones(dim,1)*abs(ids-tdim-1));
coeff =[]; x=[];

if flg_verbose ==1;
disp('rsPCA is done!');
toc;
end

%% Estimation of mean variance from each eigenvectors
kval = kurtosis(reconXbar(:,dim*2:tdim-dim*2)')-3;

%% Identification of the single-peak EV
[candidate psigpk]=determination_EVs(V,fpt,fs,FOI,window,singpdB,flg_verbose);

%% Extraction of features according to krutosis and freqeucny information
consigpk = intersect(candidate, psigpk);
ikval = find(kval<thkval);
idx=intersect(consigpk,ikval);

idx_num = 1:dim;
diffidx=setdiff(idx_num,psigpk);
iter_idx=diffidx;

nidx=length(iter_idx);

if flg_verbose ==1
disp('Corresponding eigenvectors have multiple-peaks');
disp(sprintf('%d / %d ', length(iter_idx),dim));
end

%%
h = waitbar(0, sprintf('%s (ch#%02d) is being processed ...',chnnel_info,chnnel_idx));
steps = length(iter_idx);

if isempty(iter_idx)
    
    ROIidx = setdiff(idx_num,idx);
    rXbar = sum(reconXbar(ROIidx,:));
    
    if flg_verbose==1
    disp('end');
    end
    
else
    depth = depth + 1;
    for i=1:nidx
        nsig = reconXbar(iter_idx(i),:);
        tgVar = normalD(iter_idx(i));
        
        if (stop_var <= tgVar && max_nkurt>= kval(iter_idx(i)) && depth<=max_depth)
            if flg_verbose ==1
            disp(sprintf('%s (ch#%02d): depth %02d: %02d eigenvector proceeds another recursion of the rsPCA ...',chnnel_info,chnnel_idx,depth,iter_idx(i)));
            end
            [out_rXbar] = do_irspca_recursion_func(nsig,dim,depth,loc_vec_ui32,singpdB,thkval,stop_var,chnnel_info,chnnel_idx,fs,fpt,tgVar,flg_verbose,max_depth,max_nkurt);
            reconXbar(iter_idx(i),:)= out_rXbar;
        else
            continue;
        end
        
        waitbar(i/steps);
    end
    
    ROIidx = setdiff(idx_num,idx);
    rXbar = sum(reconXbar(ROIidx,:));
end

close(h);

end

function reconXbar = subfunc_recon(loc_vec_ui32, coeff, V_t)
Xbar = coeff*V_t;
reconXbar= accumarray(loc_vec_ui32,Xbar(:));
end


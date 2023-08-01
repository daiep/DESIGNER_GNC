function [b0, dt, md, rd, ad, fa, fe, v11, v21, v31] = dti_fit_gnc_v3(dwi, grad, mask)
% DTI fitting with Gradient Nonlinearity Correction (GNC)
% (c) Erpeng Dai, Stanford University

    %% check grad order
%     bval = grad(:, 4);
%     order = floor(log(abs(max(bval)+1))./log(10));
%     if order >= 2
%         grad(:, 4) = grad(:, 4)/1000;
%         bval = grad(:, 4);
%     end
    
    % daiep 
    fprintf('... dti_fit_gnc\n');
    bval = grad(:, 4, :, :, :);
    order = floor(log(abs(max(bval(:))+1))./log(10));
    if order >= 2 % w/ GNC, bval: [ndwi, nx*ny*nz];
        grad(:, 4, :, :, :) = grad(:, 4, :, :, :)/1e3;
        bval = grad(:, 4, :, :, :);
    end
    
    %% parameter checks 
    dwi = double(dwi);
    dwi(dwi<=0)=eps;
    [x, y, z, ndwis] = size(dwi);
    
    % daiep 
    ind00=find(bval(:, 1)==0);
    b00=dwi(:, :, :, ind00);
    
    for x0=1:x
        for y0=1:y
            for z0=1:z
                b00tmp0=squeeze(b00(x0, y0, z0, :));
                b00tmp1=b00tmp0(find(b00tmp0>2*eps));
                if isempty(b00tmp1)
                    b00mean(x0, y0, z0)=0;
                else
                    b00mean(x0, y0, z0)=mean(b00tmp1);
                end
            end
        end
    end
    clear b00tmp0 b00tmp1 b00
    
    if ~exist('grad','var') || size(grad,1) ~= ndwis || size(grad,2) ~= 4
        error('');
    end
    grad = double(grad);
%     normgrad = sqrt(sum(grad(:, 1:3).^2, 2)); normgrad(normgrad == 0) = 1;
%     grad(:, 1:3) = grad(:, 1:3)./repmat(normgrad, [1 3]);
%     grad(isnan(grad)) = 0;
%     bval = grad(:, 4);
    % daiep 
    normgrad = sqrt(sum(grad(:, 1:3, :, :, :).^2, 2)); normgrad(normgrad == 0) = 1;
    grad(:, 1:3, :, :, :) = grad(:, 1:3, :, :, :)./repmat(normgrad, [1 3]);
    grad(isnan(grad)) = 0;
    bval = grad(:, 4, :, :, :);

    if ~exist('mask','var') || isempty(mask)
        mask = true(x, y, z);
    end
    
    dwi = custom_vectorize(dwi, mask);
    % daiep 
    b00mean = custom_vectorize(b00mean, mask);
    
    nvoxels = size(dwi,2);
    
    %% tensor fit
    [D_ind, D_cnt] = createTensorOrder(2);
    
%     bS = ones(ndwis, 1);
%     bD = D_cnt(ones(ndwis, 1), :).*grad(:,D_ind(:, 1)).*grad(:,D_ind(:, 2));
    
%     b = [bS, -bval(:, ones(1, 6)).*bD];
% 
%     % unconstrained LLS fit
%     dt = b\log(dwi);
%     w = exp(b*dt);
    
    % daiep 
    s0=size(bval);
    bval=bval(:, :, :);
    bS = ones([ndwis, s0(2:end)]);
    bD = D_cnt(ones(ndwis, 1), :).*grad(:,D_ind(:, 1), :, :, :).*grad(:,D_ind(:, 2), :, :, :);
    ind0=find(mask>0);
    clear grad normgrad
    
    if size(s0, 2)<=2 % w/o GNC bval: [ndwi, 1]
        b = [bS, -repmat(bval, [1 6]).*bD];
        bval = repmat(bval, [1 1  nvoxels]);
        b = repmat(b, [1 1  nvoxels]);
        % daiep  better mitigate noise bias similar to FSL dtifit
        parfor i=1:nvoxels
            ind10 = find(dwi(:, i)./b00mean(:, i) > exp(-bval(:, :, 1)*5));
            dt0(:, i) = (b(ind10, :, i)+eps)\log(dwi(ind10, i)+eps); % The sizes of dt and w keep unchanged
        end
        w = exp(b(:, :, 1)*dt0);
    else
        bval=bval(:, :, ind0); % w/ GNC, bval: [ndwi, 1, nx*ny*nz];
        bS=bS(:, :, ind0);
        bD=bD(:, :, ind0);
        b = cat(2, bS, -repmat(bval, [1 6]).*bD);
        % daiep: better mitigate noise bias similar to FSL dtifit
        parfor i=1:nvoxels
            ind10 = find(dwi(:, i)./b00mean(:, i) > exp(-bval(:, :, i)*5));
            dt0(:, i) = (b(ind10, :, i)+eps)\log(dwi(ind10, i)+eps); % The sizes of dt and w keep unchanged
            w(:, i) = exp(b(:, :, i)*dt0(:, i));
        end
    end
    % daiep
    ind20=find(isnan(dt0(2, :)) | isinf(dt0(2, :)));
    dt0(:, ind20)=[];
    w(:, ind20)=[];
    dwi(:, ind20)=[];
    b00mean(:, ind20)=[];
    mask(ind0(ind20))=0;
    ind0(ind20)=[];
    nvoxels = size(dwi,2);

    % WLLS fit initialized with LLS   
    parfor i = 1:nvoxels
        wi = diag(w(:,i)); 
        logdwii = log(dwi(:,i)+eps);
        % dt(:,i) = (wi*b)\(wi*logdwii);
        
        % daiep: better mitigate noise bias similar to FSL dtifit
        ind10 = find(dwi(:, i)./b00mean(:, i) > exp(-bval(:, :, i)*5));
        dt0(:,i) = (wi(ind10, ind10)*b(ind10, :, i)+eps)\(wi(ind10, ind10)*logdwii(ind10, :)+eps);
    end
    clear dwi
    
    % daiep
    ind20=find(isnan(dt0(2, :)) | isinf(dt0(2, :)));
    dt0(:, ind20)=[];
    mask(ind0(ind20))=0;
    ind0(ind20)=[];
    
    b0 = exp(dt0(1,:));
    dt = dt0(2:7, :);
    
    b0 = custom_vectorize(b0, mask);
    
    DT = reshape([dt(1,:); dt(2,:); dt(3,:);...
        dt(2,:); dt(4,:); dt(5,:);...
        dt(3,:); dt(5,:); dt(6,:)],[3,3,size(dt,2)]);

    parfor i=1:size(dt,2)
    % for i=1:size(dt,2)
        [v,l] = eig(DT(:,:,i));
        l = diag(l); 
        [l, idx] = sort(l,'descend'); 
        v = v(:, idx);
        
        v1(i,:)=v(:,1);
        v2(i,:)=v(:,2);
        v3(i,:)=v(:,3);
        l1(i)=abs(l(1));
        l2(i)=abs(l(2));
        l3(i)=abs(l(3));
    end
    dt = custom_vectorize(dt, mask);
    md = custom_vectorize((l1+l2+l3)./3, mask);
    rd = custom_vectorize((l2+l3)./2, mask);
    ad = custom_vectorize(l1, mask);
    fa = custom_vectorize(sqrt(1/2).*sqrt((l1-l2).^2+(l2-l3).^2+(l3-l1).^2)./sqrt(l1.^2+l2.^2+l3.^2+eps), mask);
    % daiep
    v11 = custom_vectorize(v1', mask);
    v21 = custom_vectorize(v2', mask);
    v31 = custom_vectorize(v3', mask);
    
    for i = 1:3
        ev(:,:,:,i) = custom_vectorize(v1(:,i)', mask);
    end
    fe = cat(4,fa.*ev(:,:,:,1),fa.*ev(:,:,:,2),fa.*ev(:,:,:,3));

end

function [X, cnt] = createTensorOrder(order)
    X = nchoosek(kron([1, 2, 3], ones(1, order)), order);
    X = unique(X, 'rows');
    for i = 1:size(X, 1)
        cnt(i) = factorial(order) / factorial(nnz(X(i, :) ==1))/ factorial(nnz(X(i, :) ==2))/ factorial(nnz(X(i, :) ==3));
    end

end


function [s, mask] = custom_vectorize(S, mask)
    if nargin == 1
       mask = ~isnan(S(:,:,:,1));
    end
    if ismatrix(S)
        n = size(S, 1);
        [x, y, z] = size(mask);
        % s = NaN([x, y, z, n], 'like', S);
        % daiep 
        s = zeros([x, y, z, n], 'like', S);
        for i = 1:n
            % tmp = NaN(x, y, z, 'like', S);
            % daiep 
            tmp = zeros(x, y, z, 'like', S);
            tmp(mask(:)) = S(i, :);
            s(:,:,:,i) = tmp;
        end
    else
        for i = 1:size(S, 4)
           Si = S(:,:,:,i);
           s(i, :) = Si(mask(:));
        end
    end
end


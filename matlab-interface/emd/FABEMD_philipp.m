% -------------function emdc=FABEMD(signal,parameter,fid,iy)---------------
% Purpose :
% compute the empirical mode decomposition of a 2D signal and save some 
% computation infos in a file fid
% Input: 
% signal: 2D matrix (m*n) of real doubles
% parameter: structure, contains the parameters for the DNS and the FABEMD
% parameter.MNAI= maximum number of allowed iterations
% parameter.type: string, indicates the windowing width
% - 'HD': highest distance
% - 'LD': lowest distance
% - 'hD': the min of the maximal distances in upper and lower envelopes
% Output :
% emdc= empirical mode decomposition components, matrix(m*n*(N+1))

function emdc=FABEMD(signal,parameter,verbose)

% verification of the valitidy of the input parameter
if (strcmp(parameter.type,'HD')==false && strcmp(parameter.type,'LD')==false && strcmp(parameter.type,'hD')==false)
    error('error in the Input argument type. Allowable values: HD or LD or hD')
end

if strcmp(parameter.type,'LD')
    ws=@(a,b) 2*floor(min(min(a),min(b))/2)+1;
elseif strcmp(parameter.type,'hD')
    ws=@(a,b) 2*floor(min(max(a),max(b))/2)+1;
elseif strcmp(parameter.type,'HD')
    ws=@(a,b) 2*floor(max(max(a),max(b))/2)+1;
end
disp('preallocating emdc');
fprintf('running on %i cores \n',maxNumCompThreads);
emdc=repmat(zeros(size(signal)),[1,1,parameter.n+1]);
disp('preallocating done');

% 1) initialisation of the residual: res=signal
residual=signal;
w_old=0;
for i=1:parameter.n
    res=residual;
    j=1;
    
    while j<=parameter.MNAI
	disp('finding min/max');
        tstart_tot=tic;
        tstart=tic;
        % 2) upper and lower envelopes

        indexf=imregionalmax(res');
        [r,c]=find(indexf==1);
        Xmax = [c r];
        indexf=imregionalmin(res');
        [r,c]=find(indexf==1);
        Xmin = [c r];
        t1=toc(tstart);
	disp('finding min/max done');
        
	disp('distance check');

        tstart=tic;
        [~,d_max]=knnsearch(Xmax,Xmax,'k',2);
        [~,d_min]=knnsearch(Xmin,Xmin,'k',2);t2=toc(tstart);
	disp('distance check done');
        
        if (size(d_min,1)>=2 && size(d_max,1)>=2)
            w = ws(d_max(:,2),d_min(:,2));
            
            if w<=w_old
                w=2*floor(1.5*w_old/2)+1; 
            end
            

            tstart2=tic;
            if w<floor(min(size(signal))/2)
                disp('envelope calculation');
                tstart=tic;
                xLayer = (w-1)/2;       
                env_max = ordfilt2(res,w*w,ones(w,w),'symmetric');
                env_min = ordfilt2(res,1,ones(w,w),'symmetric');
                t3=toc(tstart);
                disp('bp1')                
                tstart=tic;
                h = fspecial('average',w);
                env_max_pad = padarray(env_max,[xLayer  xLayer],'circular');
                env_max = conv2(env_max_pad,h,'same');
                env_max = env_max(xLayer+1:end-xLayer, xLayer+1:end-xLayer);

%                 env_max = conv2(env_max,h,'same');

% 
                h = fspecial('average',w);
                env_min_pad = padarray(env_min,[xLayer xLayer],'circular');
                env_min = conv2(env_min_pad,h,'same');
                env_min = env_min(xLayer+1:end-xLayer, xLayer+1:end-xLayer);
                t4=toc(tstart);
                disp('envelope calculation done');

%                 env_max=smooth(env_max,w);
%                 tstart=tic;env_min=smooth(env_min,w);tsmin=toc(tstart);
                % 3) computation of the mean of the two envelopes env_mean
                disp('mean env and new res');

                env_mean=(env_max+env_min)/2;
                % 4) computation of the new residual : res=res-env_mean
                res=res-env_mean;
                disp('new res obtained');

                %                 if verbose
%                     tstart=tic;[Xmax,~]=maximum_local_BC(res);t5=toc(tstart5);
%                     tstart=tic;[Xmin,~]=minimum_local_BC(res);t6=toc(tstart);
%                     mean=sum(sum(res))/(size(res,1)*size(res,2));
%                     disp(iy,i,j,w,size(Xmin,1),size(Xmax,1),toc(tstart_tot),mean,t2,t4,tfilter,tsmin);
%                 end
                j=j+1;
            else
                emdc(:,:,i)=res;
                return;
            end
        else 
            emdc(:,:,i)=res;
            return;
        end
        ttot=toc(tstart_tot);
        disp(sprintf('%i: t1=%0.5g t2=%0.5g t3=%0.5g t4=%0.5g tt=%0.5g w=%i',i,t1,t2,t3,t4,t1+t2+t3+t4,w));
    end
    % 6) a new residual is defined: res=res-IMF. Repeat from 2).
    emdc(:,:,i)=res;
    residual=residual-res;
    w_old=w;
end

emdc(:,:,i+1)=residual;

end

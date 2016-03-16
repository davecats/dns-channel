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

function emdc=FABEMD(signal,parameter,fid,iy)

% verification of the valitidy of the input parameter
if (strcmp(parameter.type,'HD')==false && strcmp(parameter.type,'LD')==false && strcmp(parameter.type,'hD')==false)
    error('error in the Input argument type. Allowable values: HD or LD or hD')
end

emdc=repmat(zeros(size(signal)),[1,1,parameter.n+1]);

% 1) initialisation of the residual: res=signal
residual=signal;
w_old=0;
for i=1:parameter.n
    res=residual;
    j=1;
    
    while j<=parameter.MNAI
        tstart_tot=tic;
        % 2) upper and lower envelopes
        [Xmax,Valmax]=maximum_local_BC(res); 
        tstart=tic;[Xmin,Valmin]=minimum_local_BC(res);tmin=floor(toc(tstart));
        [k_max,d_max]=knnsearch(Xmax,Xmax,'k',2); 
        tstart=tic;[k_min,d_min]=knnsearch(Xmin,Xmin,'k',2);tdmin=floor(toc(tstart));
        
        if (size(d_min,1)>=2 && size(d_max,1)>=2)
            if parameter.type=='LD'
                w=2*floor(min(min(d_max(:,2)),min(d_min(:,2)))/2)+1;
            elseif parameter.type=='hD'
                w=2*floor(min(max(d_max(:,2)),max(d_min(:,2)))/2)+1;
            elseif parameter.type=='HD'
                w=2*floor(max(max(d_max(:,2)),max(d_min(:,2)))/2)+1;
            end
            
            if w<=w_old
                w=2*floor(1.5*w_old/2)+1; 
            end
            
            if w<floor(min(size(signal))/2)
                tstart=tic;[env_max,env_min]=extrema_filter(res,w);tfilter=floor(toc(tstart));
                env_max=smooth(env_max,w);
                tstart=tic;env_min=smooth(env_min,w);tsmin=floor(toc(tstart));
                % 3) computation of the mean of the two envelopes env_mean
                env_mean=(env_max+env_min)/2;
                % 4) computation of the new residual : res=res-env_mean
                res=res-env_mean;
                [Xmax,~]=maximum_local_BC(res);
                [Xmin,~]=minimum_local_BC(res);
                mean=sum(sum(res))/(size(res,1)*size(res,2));
                %fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%.3f\t%d\t%d\t%d\t%d\n',iy,i,j,w,size(Xmin,1),size(Xmax,1),toc(tstart_tot),mean,tmin,tdmin,tfilter,tsmin);
                j=j+1;
            else
                emdc(:,:,i)=res;
                return;
            end
            
        else 
            emdc(:,:,i)=res;
            return;
        end
        
    end
    % 6) a new residual is defined: res=res-IMF. Repeat from 2).
    emdc(:,:,i)=res;
    residual=residual-res;
    w_old=w;
end

emdc(:,:,i+1)=residual;
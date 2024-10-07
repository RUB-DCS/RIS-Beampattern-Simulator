function [configs] = maximizeRIS_CB(params)
    channelDir  = params.chanDir;
    channelVec1 = params.chan1;
    channelVec2 = params.chan2;
        
    h_eff_mat= channelVec1.*channelVec2;

    C_all = linspace(0,2*pi,3600);
    C = C_all([(0:params.phaseResolution-1)*floor(3600/params.phaseResolution)]+1);

    switch params.direct_channel_property
        case 'blocked'

            phases = [];

            vector=C;
            matrix=(angle(h_eff_mat));

            % Replicate the vector along the third dimension to match the size of the matrix
            vector_replicated = repmat(vector, size(matrix, 1), size(matrix, 2), 1);
            
            % Perform element-wise subtraction
            phases = permute(repmat(vector', 1, size(matrix, 1), size(matrix, 2)) - permute(matrix, [3, 1, 2]),[2, 3, 1]);
            
            h_test_mat = h_eff_mat.*exp(1i*phases);

       
            optPhase_quant=[];
            optPhase_quant = angle(exp(1i*phases));
            optPhase_quant(abs(optPhase_quant)<pi/2)=0;
            optPhase_quant(abs(optPhase_quant)>=pi/2)=pi;
                            
          
            %include db loss
            scalMult = 0.*phases+1;
            scalMult(find(optPhase_quant)) = db2pow(-3);
            h_test_quant= (h_eff_mat.*scalMult).*exp(1i*optPhase_quant);
           



            if ~isfield(params,'TxSweepPoints')
                if isfield(params,'RxSplitPos')
                  
%                     conf = permute(reshape(permute(optPhase_quant,[1,3,2]),64,48,2),[2,1,3]);
%                     figure
%                     imagesc(conf(:,:,1))
%                     figure
%                     imagesc(conf(:,:,2))

                    %select left/right half of RIS    
                    selMat1 = zeros(params.N_surfaces_Y*params.N_y,params.N_surfaces_X*params.N_x);
                    selMat2 = zeros(params.N_surfaces_Y*params.N_y,params.N_surfaces_X*params.N_x);
                    selMat1(:,1:ceil(params.N_surfaces_X*params.N_x/2)) = 1;
                    selMat2(:,floor(params.N_surfaces_X*params.N_x/2):end) = 1;
                    selMat1 = reshape(selMat1.',[],1);
                    selMat2 = reshape(selMat2.',[],1);



                else
                   [bestPhaseDirVals,bestPhaseDirInd] = max(abs(squeeze(sum(h_test_quant,1))),[],1);
                   %bestPhaseDirInd = 4;
                   optPhase_quant_out = squeeze(optPhase_quant(:,:,bestPhaseDirInd));
                end
                  % h_test_quant_out(:,k) = h_test_quant(:,k,bestPhaseDirInd(k));
            else
                [bestPhaseDirVals,bestPhaseDirInd] = max(abs(squeeze(sum(h_test_quant,1))),[],2);
                for k = 1:length(params.TxSweepPoints)
                   optPhase_quant_out(:,k) = optPhase_quant(:,k,bestPhaseDirInd(k));
                  % h_test_quant_out(:,k) = h_test_quant(:,k,bestPhaseDirInd(k));
                end
            end
        case 'available'
            optPhase = angle(channel_Dir)-(angle(channelVec1)+angle(channelVec2)).';
            h_test2=channel_Dir+channelVec1*diag(exp(1i*optPhase))*channelVec2.';
            power2=abs(h_test2)^2;

            optPhase_quant = angle(exp(1i*optPhase));
            optPhase_quant(abs(optPhase_quant)<pi/2)=0;
            optPhase_quant(abs(optPhase_quant)>=pi/2)=pi;

            h_test_quant=channel_Dir+channelVec1*diag(exp(1i*optPhase_quant))*channelVec2.';
            power_quant=abs(h_test_quant)^2;
    end
    
    configs.phaseOffsets =  C;
    configs.cont  = phases;
    configs.quant = optPhase_quant_out;
    configs.bestPhaseInd = bestPhaseDirInd;
    
%     gains.quantChan = h_test_quant;
%     gains.linEq = gains_linEq;
end


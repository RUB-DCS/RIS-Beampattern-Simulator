function [configs,gains] = maximizeRIS(params,channelParams)
    channelDir  = channelParams.dir;
    channelVec1 = channelParams.chan1;
    channelVec2 = channelParams.chan2;
        
    
    %%tiling
    indVecTile = channelParams.indVecTile;
    cloneIndVec = channelParams.indVecTile;
    tileSize = channelParams.tileSize;
    
    h_eff_tiled = channelVec1.*channelVec2;
%     if tileSize==1
%        save('h_eff','h_eff')
%     end
%     for k = 1 : (16*16)/tileSize
%       h_eff_tiled(k) = sum(h_eff(cloneIndVec(1:tileSize)));
%       cloneIndVec(1:tileSize)=[];
%     end
    %%end tiling

    h_eff = h_eff_tiled;

    switch params.direct_channel_property
        case 'blocked'
            C_all = linspace(0,2*pi,3600);
            C = C_all([(0:params.polarPoints-1)*floor(3600/params.polarPoints)]+1);
            for c = C
                c_ind = find(c==C);
                
                phases = [];
                phases_tiled = c - (angle(h_eff_tiled)).';
                
               % phases(indVecTile,1) = kron(phases_tiled,ones(tileSize,1));
                phases = kron(phases_tiled,ones(tileSize,1));


                
                h_test = h_eff*exp(1i*phases);
                
               
                
                optPhase_quant=[];
                optPhase_quant = angle(exp(1i*phases));
                optPhase_quant(abs(optPhase_quant)<pi/2)=0;
                optPhase_quant(abs(optPhase_quant)>=pi/2)=pi;
                                
                gains_cont(c_ind)=abs(h_test)^2;
                optPhases{c_ind} = phases;
                
                scalMult = ones(length(h_eff),1);
                scalMult(find(optPhase_quant)) = db2pow(-3);
                %scalMult(1:256) = (1+0.5.*sin([1:256]./256)) .*exp(1i*2*pi*linspace(0,1,256));
                h_test_quant= (h_eff.*scalMult.')*exp(1i*optPhase_quant);
                
               % h_test_quant= h_eff*exp(1i*optPhase_quant);
                gains_linEq(c_ind)= h_test_quant;
                gains_quant(c_ind)=abs(h_test_quant)^2;
                optQuantPhases{c_ind} = optPhase_quant;   
            end     
            
%             for index = 1:256
%                optPhase_quant=optPhase_quant.*0;
%                optPhase_quant(index)=pi;
%                scalMult = ones(length(h_eff),1);
%                scalMult(find(optPhase_quant)) = db2pow(-3);
%                h_test_quant= (h_eff.*scalMult.')*exp(1i*optPhase_quant);
%                gains_linEq(c_ind+index)= h_test_quant;
%                gains_quant(c_ind+index)=abs(h_test_quant)^2;
%                 optQuantPhases{c_ind+index} = optPhase_quant;  
%             end
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
    configs.cont  = optPhases;
    configs.quant = optQuantPhases;
    
    gains.phaseOffsets = C;
    gains.channelParams  = channelParams;
    gains.cont  = gains_cont;
    gains.quant = gains_quant;
    gains.linEq = gains_linEq;
end

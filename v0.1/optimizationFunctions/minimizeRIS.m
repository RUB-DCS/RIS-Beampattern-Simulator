function [configs,gains] = minimizeRIS(params,channelParams)
    channelDir  = channelParams.dir;
    channelVec1 = channelParams.chan1;
    channelVec2 = channelParams.chan2;
        
    %%tiling
    indVecTile = channelParams.indVecTile;
    cloneIndVec = channelParams.indVecTile;
    tileSize = channelParams.tileSize;
    
    h_eff = channelVec1.*channelVec2;

    for k = 1 : (16*16)/tileSize
      h_eff_tiled(k) = sum(h_eff(cloneIndVec(1:tileSize)));
      cloneIndVec(1:tileSize)=[];
    end
    %%end tiling
    
    h_eff_tiled = h_eff_tiled*10^10;

    kappa=0.0000001;
    gamma=0.05;
    chanval=[];
    switch params.direct_channel_property
        case 'blocked'
            
            v_tild = real(exp(1i*2*pi*rand((16*16)/tileSize,1)));
            v = (v_tild.*exp(1i*2*pi*rand((16*16)/tileSize,1)));
            count = 0;
            oldChan = 0;
            newChan = 1;
            kappa2=1000;
            optPhase_quant = angle(v_tild);
            optPhase_quant(abs(optPhase_quant)<pi/2)=0;
            optPhase_quant(abs(optPhase_quant)>=pi/2)=pi;
             for k = 1:(256/tileSize)
                cumChan2(k) = h_eff(1:k)*exp(1i*(optPhase_quant(1:k)));
             end
            
            fMin = figure(2);   
            fMin.Position = [1300,500,500 500];
            
            best.v = [];
            best.chan=1;
            while ((abs(oldChan-newChan)>0.001 || ~(16*16-1 < sum(abs(v)) && sum(abs(v))<16*16+1))) || kappa<10^3
            
            oldChan = cumChan2(end);
            count = count+1;
            if kappa<10^8
               kappa=kappa*1.5;
            else
              break
            end
            cvx_solver mosek
            cvx_begin quiet
          
            variable v((16*16)/tileSize,1)        %phase shifters
            
            minimize(  0.0000001*pow_abs(h_eff_tiled*v,2)+ 0.9999999*kappa /((16*16)/tileSize) * sum( real (kappa*conj(v_tild).*v + pow_abs(v_tild,2) ))  ) ;
                 
            subject to
              % pow_pos(v)>=1
                abs(v) <= 1;
            cvx_end   

               
            v_tild(v<0)=-gamma;% = v_tild + min(gamma,1)*(v-v_tild);%;exp(1i*angle(v));
            v_tild(v>0)=gamma;
            gamma=min(gamma+0.05,1);
            
            phases(indVecTile,1) = kron(v_tild,ones(tileSize,1));
            phases_cont(indVecTile,1) = kron(v,ones(tileSize,1));
            
%             for k = 1:256
%                 cumChan(k) = h_eff(1:k)*exp(1i*(optPhase_quant(1:k)));
%              end
            
            
%             figure(1) 
%             polarplot(cumChan)
%             
            optPhase_quant = angle(phases);
            optPhase_quant(abs(optPhase_quant)<pi/2)=0;
            optPhase_quant(abs(optPhase_quant)>=pi/2)=pi;
           % hold on
             for k = 1:(256/tileSize)
                cumChan2(k) = h_eff(1:k)*exp(1i*(optPhase_quant(1:k)));
             end
            
            try
             % polarplot(cumChan2)
              %hold on
            catch
%                pause(0.5)
%                polarplot(cumChan2)
            end
               
             newChan = cumChan2(end);
             pause(0.01)
             
             channelVec1 = channelParams.chan1;
            channelVec2 = channelParams.chan2;
    
            h_test=channelVec1*diag(phases_cont)*channelVec2.';  
            gains_cont=abs(h_test)^2;
            
            h_test=channelVec1*diag(exp(1i*optPhase_quant))*channelVec2.';  
            gains_quant=abs(h_test)^2;
            
             if best.chan>gains_quant
                  best.chan = gains_quant;
                  best.v = optPhase_quant;
             end
             chanval=[chanval,gains_quant];
             if mod(count,5)==0
             plot(pow2db(chanval))
             end
             %hold on
            end
             
%             if 16*16-0.1 < sum(abs(v)) || sum(abs(v))<16*16+0.1))
%                disp('+ min')
%             end
            
            if sum(abs(v))<((16*16)/tileSize)-1
               disp('- failed min')
            end
            
            channelVec1 = channelParams.chan1;
            channelVec2 = channelParams.chan2;
    
            h_test=channelVec1*diag(phases_cont)*channelVec2.';  
            gains_cont=abs(h_test)^2;
            
            h_test=channelVec1*diag(phases)*channelVec2.';  
            gains_quant=abs(h_test)^2;

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
    
    optPhase_quant=best.v;
    h_test=channelVec1*diag(exp(1i*optPhase_quant))*channelVec2.'; 
    
    gains_quant=abs(h_test)^2;
    
    
    optPhases{1} = v;
    optQuantPhases{1} = optPhase_quant;

    configs.phaseOffsets =  [];
    configs.cont  = optPhases;
    configs.quant = optQuantPhases;
    
    gains.phaseOffsets = [];
    gains.cont  = gains_cont;
    gains.quant = gains_quant;
end


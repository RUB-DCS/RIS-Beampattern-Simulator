function [configs,gains] = maxquantRIS(params,channelParams)
    channelDir  = channelParams.dir;
    channelVec1 = channelParams.chan1;
    channelVec2 = channelParams.chan2;
    
    channelVec1=channelVec1.*10^3;
    channelVec2=channelVec2.*10^3;
    
    kappa=0.01;
    switch params.direct_channel_property
        case 'blocked'
            
            v_tild = exp(1i*2*pi*ones(16*16,1));
            v = v_tild*exp(1i*pi);
            count = 0;
            oldChan = 0;
            newChan = 1;
            
            optPhase_quant = angle(v_tild);
            optPhase_quant(abs(optPhase_quant)<pi/2)=0;
            optPhase_quant(abs(optPhase_quant)>=pi/2)=pi;
             for k = 1:256
                cumChan2(k) = channelVec1(1:k)*diag(channelVec2(1:k))*exp(1i*(optPhase_quant(1:k)));
             end
            
            
            
            while abs(oldChan-newChan)>0.001
            
            oldChan = cumChan2(end);
            count = count+1;
            
            cvx_solver 
            cvx_begin 
          
%             variable R(1,K) nonnegative   %rates
%             variable tp(1,K) nonnegative  %sinrs
            variable v(16*16,1) complex       %phase shifters
            
            maximize( abs((channelVec1*diag(channelVec2)*v)' + 2*kappa*sum(real(2*v_tild.*(v-v_tild))) ) % fairness  - penalty
                 
            subject to
            
               abs(v)<=1;
                
            cvx_end    
            
            v_tild=v;
            
            
            for k = 1:256
                cumChan(k) = channelVec1(1:k)*diag(channelVec2(1:k))*v(1:k);
            end
            
            
%             figure(1) 
%             polarplot(cumChan)
%             
            optPhase_quant = angle(v);
            optPhase_quant(abs(optPhase_quant)<pi/2)=0;
            optPhase_quant(abs(optPhase_quant)>=pi/2)=pi;
           % hold on
             for k = 1:256
                cumChan2(k) = channelVec1(1:k)*diag(channelVec2(1:k))*exp(1i*(optPhase_quant(1:k)));
             end
%              polarplot(cumChan2)
    
             newChan = cumChan2(end);
             
            end
            
            h_test=channelVec1*diag(v)*channelVec2.';  
            gains_cont=abs(h_test)^2;
            
            h_test=channelVec1*diag(v)*channelVec2.';  
            gains_quant=abs(cumChan2(end))^2;
           
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
    
    configs.phaseOffsets =  [];
    configs.cont  = v;
    configs.quant = optPhase_quant;
    
    gains.phaseOffsets = [];
    gains.cont  = gains_cont;
    gains.quant = gains_quant;
end


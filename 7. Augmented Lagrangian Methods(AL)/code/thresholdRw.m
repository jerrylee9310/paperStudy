function [v] = thresholdRw(Rw, params, mu)

Operator = params.Operator;
lambdaTV = params.lambda.tv;
lambdaWavelet = params.lambda.wavelet;

% dataTV = Rw(:,:,1:2);
% dataWav = Rw(:,:,3:2+params.Wavelet.dnum); 
% 
% tmpTV = sos(dataTV);
% tmpWav = abs(dataWav);
% 
% % TV
% dTV = tmpTV - lambdaTV/mu;
% dTV = (dTV + abs(dTV))/ 2;
% ddTV = dataTV./tmpTV.*dTV;
% 
% 
% % wavelet
% dWav = tmpWav - lambdaWavelet/mu;
% dWav = (dWav + abs(dWav))/ 2;
% ddWav = dataWav./tmpWav.*dWav;
% 
% 
% v = cat(3,ddTV,ddWav);
% v(isnan(v)) = 0; % Ensure v does not contain NaN obtained from 0/0 operations


%% Perform thresholding
switch(Operator)
    case{'FD'} % Finite difference
        dataTV = Rw;
        tmpTV = sos(dataTV);
        % TV
        dTV = tmpTV - lambdaTV/mu;
        dTV = (dTV + abs(dTV))/ 2;
        ddTV = dataTV./tmpTV.*dTV;
        
        v = ddTV;

    case{'W'} % Wavelets
        dataWav = Rw;
        tmpWav = abs(dataWav);
        
        % wavelet
        dWav = tmpWav - lambdaWavelet/mu;
        dWav = (dWav + abs(dWav))/ 2;
        ddWav = dataWav./tmpWav.*dWav;
        
        v = ddWav;

    case{'WFD'}
        dataTV = Rw(:,:,1:2);
        dataWav = Rw(:,:,3:2+params.Wavelet.dnum); 

        tmpTV = sos(dataTV);
        tmpWav = abs(dataWav);

        % TV
        dTV = tmpTV - lambdaTV/mu;
        dTV = (dTV + abs(dTV))/ 2;
        ddTV = dataTV./tmpTV.*dTV;


        % wavelet
        dWav = tmpWav - lambdaWavelet/mu;
        dWav = (dWav + abs(dWav))/ 2;
        ddWav = dataWav./tmpWav.*dWav;
        
        v = cat(3,ddTV,ddWav);
end

v(isnan(v)) = 0;



% %% mycode
% % wavelet
% tolwav = lambdaWavelet/mu;
% dWav = dataWav./abs(dataWav).*max(abs(dataWav)-tolwav,0);
% 
% % TV
% toltv = lambdaTV/mu;
% dTV = dataTV./sos(dataTV).*max(sos(dataTV)-toltv,0);
% 
% v = cat(3,dTV,dWav);
% v(isnan(v)) = 0;
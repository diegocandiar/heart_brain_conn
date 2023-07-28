function [conn, mrkv] = arx_model(EEG_comp, Fs, timearx, wind)

% Author: Diego Candia-Rivera (diego.candia.r@ug.uchile.cl)
%%
[Nch,Nt] = size(EEG_comp);

% Optional scaling for model convergence
EEG_comp(1:Nch-2, :) = sqrt(EEG_comp(1:Nch-2, :));
EEG_comp(Nch-1, :) = EEG_comp(Nch-1, :) / std(EEG_comp(Nch-1, :));
EEG_comp(Nch, :) = EEG_comp(Nch, :) / std(EEG_comp(Nch, :));

conn = zeros(Nch,Nch,length(timearx));
mrkv = zeros(Nch,Nch,length(timearx));

parfor ch = 1:Nch % Consider use parallel computing here
    [conn(ch,:,:), mrkv(ch,:,:)] = arx_ch(EEG_comp, ch, wind*Fs);
end

end

function [conn_ch, mrkv_ch] = arx_ch(EEG_ch, ch1, window)
    [Nch,Nt] = size(EEG_ch);
    EEG_ch1 = EEG_ch(ch1,:);
    conn_ch = ones(Nch, Nt-window);
    mrkv_ch = zeros(Nch, Nt-window);

    for ch2 = 1 : Nch
        if ch1 == ch2
            continue;
        end
        EEG_ch2 = EEG_ch(ch2,:);
        for i = 1 : Nt-window
            arx_data = iddata(EEG_ch1(i:i+window)', EEG_ch2(i:i+window)',1); 
            model_eegP = arx(arx_data,[1 1 1]); 
            conn_ch(ch2,i) = model_eegP.B(2); 
            mrkv_ch(ch2,i) = model_eegP.A(2);         
        end
    end
end
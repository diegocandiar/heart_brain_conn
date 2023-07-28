function MIC = compute_coupling(conn,t_conn, hrv, t_hrv)


[nch, nch2, nt] = size(conn);
t1 = max([ t_conn(1) t_hrv(1) ]);
t2 = min([ t_conn(end) t_hrv(end) ]);

st = 1;
t_f = t1:st:t2;
    
sig2_hrv = interp1(t_hrv, hrv, t_f, 'spline');
MIC = zeros(nch,nch);
for ch1 = 1 : nch
    for ch2 = 1 : nch
        if ch1 == ch2
            continue;
        else
            X = conn(ch1, ch2, :);
            sig1 = interp1(t_conn, X, t_f, 'spline');
            sig1 = zscore(sig1); 
            sig2_hrv = zscore(sig2_hrv); 
            sig2_hrv = zscore(sig2_hrv);

            stats = MICfunction(sig1, sig2_hrv);
            MIC(ch1,ch2) = stats(1);
        end
    end
end

end
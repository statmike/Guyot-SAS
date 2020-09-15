/* NOTES:

This code is from 

    CHANGES FROM ORIGINAL:
        - removed write_path variables and actions
    
    FINISH UP:
        header for this file
        readme.md file with examples

*/

proc iml;

	start guyot_method(x, y, t, r, tot_events=.);

        /* 44-58 */
        /* check for error and warning states */
            if max(y) > 1 then do;
                print "Survival (y coordinates) can't be greater than 1. Please scale your y to be in [0, 1] range.";
                abort;
            end;
            if min(x) > 0 then do;
                print "Adding coordinate (x=0, y=1) to your data.";
                x = 0 || x;
                y = 1 || y;
            end;
            if min(t) > 0 then print "Warning: Number at risk at time 0 is needed. The calculation might not work correctly otherwise.";
            if max(t) > max(x) then do;
                print "Warning: Times for number at risk (t) where t > x (x = coordinates of the curve) are not used.";
                r = r[,loc(t<max(x))];
                t = t[,loc(t<max(x))];
            end;

        /* 59-60 */
        /* define variables */
            t_S = x;
            S = y;
            n_risk = r;
            t_risk = t;

        /* 62-76 */
        /* Need to grab `lower` and `upper`, to "match" vector x against vector t
            See Guyot Table 2 for visual explanation of this */
            lower = j(1, ncol(t));
				do i = 1 to ncol(lower);
					lower[i] = min(loc(x >= t[i]));
				end;
			upper = j(1, ncol(t));
				temp = t[1,2:ncol(t)] || max(x)+1;
				do i = 1 to ncol(upper);
					upper[i] = max(loc(x < temp[i]));
				end;
				*free temp;

            /* Initialise vectors */
                n_int = ncol(n_risk);
                n_t = upper[n_int];
                
                n_censor = j(1,n_int-1,0);
                n_hat = j(1,n_t,n_risk[1]+1);
                cen = j(1,n_t,0);
                d = j(1,n_t,0);
                KM_hat = j(1,n_t,1);
                last_i = j(1,n_int,1);
                sumdL = 0;

        /* 77-118 */
            if n_int > 1 then do;
                /* loop time intervals 1 to n_int-1 */
                do i = 1 to n_int-1;
                    /* start with approximation of number of censored on interval i */
                        n_censor[i] = round(n_risk[i]*S[lower[i+1]]/S[lower[i]]- n_risk[i+1]);
                    /* adjust total number of censored until n_hat = n_risk at start of interval i+1 */
                        do while ((n_hat[lower[i+1]] > n_risk[i+1]) | ((n_hat[lower[i+1]] < n_risk[i+1])  &  (n_censor[i] > 0)));
                            if n_censor[i] <= 0 then do;
                                cen[lower[i]:upper[i]] = 0;
                                n_censor[i] = 0;
                            end;
                            if n_censor[i] > 0 then do;
                                cen_t = j(1,n_censor[i],0);
                                do j = 1 to n_censor[i];
                                    cen_t[j] = t_S[lower[i]] + j * (t_S[lower[i+1]] - t_S[lower[i]]) / (n_censor[i]+1);
                                end;
                                /* Distribute censored observations evenly over time */
                                        b = bin(cen_t,t_S[lower[i]:lower[(i+1)]]);
                                        call tabulate(levels, freq, b);
                                        do f = 1 to ncol(levels);
                                            cen[(lower[i]+f-1)] = freq[levels[f]];
                                        end;                            
                            end;
                            /* find number of events and number at risk on each interval to agree with KM estimates read from curve */
                                n_hat[lower[i]] = n_risk[i];
                                last = last_i[i];
                                do k = lower[i] to upper[i];
                                    if (i = 1 & k = lower[i]) then do;
                                        d[k] = 0;
                                        KM_hat[k] = 1;
                                    end;
                                    else do;
                                        d[k] = round(n_hat[k] * (1 - (S[k] / KM_hat[last])));
                                        KM_hat[k] = KM_hat[last] * (1 - (d[k] / n_hat[k]));
                                    end;
                                    n_hat[k+1] = n_hat[k] - d[k] - cen[k];
                                    if d[k] ^= 0 then last = k;                                    
                                end;
                                n_censor[i] = n_censor[i] + (n_hat[lower[i+1]] - n_risk[i+1]);
                        end;
                    /* store the results in last_i and n_hat */
                        if n_hat[lower[i+1]] < n_risk[i+1] then n_risk[i+1] = n_hat[lower[i+1]];
                        last_i[i+1] = last;
                end;
            end;

        /* 119-138 */
        /* Time interval n_int */
            if n_int > 1 then do;
                /* Assume same censor rate as average over previous time intervals. */
                n_censor = n_censor || min(round(sum(n_censor[1:(n_int-1)]) * (t_S[upper[n_int]] - t_S[lower[n_int]]) / (t_S[upper[(n_int-1)]] - t_S[lower[1]])), n_risk[n_int]);
            end;
            if n_int = 1 then n_censor = n_censor || 0;
            if n_censor[n_int] <= 0 then do;
                cen[lower[n_int]:(upper[n_int]-1)] = 0;
                n_censor[n_int] = 0;
            end;
            if n_censor[n_int] > 0 then do;
                cen_t = j(1,n_censor[n_int],0);
                do j = 1 to n_censor[n_int];
                    cen_t[j] = t_S[lower[n_int]] + j * (t_S[upper[n_int]] - t_S[lower[n_int]]) / (n_censor[n_int] + 1);
                end;
                /* Distribute censored observations evenly over time */
                    b = bin(cen_t,t_S[lower[n_int]:upper[n_int]]);
                    call tabulate(levels, freq, b);
                    do f = 1 to ncol(levels);
                        cen[(lower[n_int]+f-1)] = freq[levels[f]];
                    end;
            end;

        /* 139-153 */
        /* Find number of events and number at risk on each interval to agree with KM estimates read from curves */
            n_hat[lower[n_int]] = n_risk[n_int];
            last = last_i[n_int];
            do k = lower[n_int] to upper[n_int];
                if KM_hat[last] ^= 0 then d[k] = round(n_hat[k] * (1 - (S[k] / KM_hat[last])));
                else d[k] = 0;
                /* condition used here avoid div 0 - original R code makes div/0 NAN, SAS does also but gives a warning */
                if n_hat[k] ^=0 then KM_hat[k] = KM_hat[last] * (1 - (d[k] / n_hat[k]));
                else KM_hat[k]=.;
				if ncol(n_hat) = k then n_hat = n_hat || (n_hat[k] - d[k] - cen[k]);
				else n_hat[k+1] = n_hat[k] - d[k] - cen[k];
                /* the number at risk cannot be negative ... */
                if n_hat[k+1] < 0 then do;
                    n_hat[k+1] = 0;
                    cen[k] = n_hat[k] - d[k];
                end;
                if d[k] ^= 0 then last = k;
            end;

        /* 154-204 */
        /* If total number of events reported (input variable tot_events=), adjust number censored so that total number of events agrees. */
            if tot_events then do;
                if n_int > 1 then do;
                    sumdL = sum(d[1:upper[(n_int-1)]]);
                    /* If total number of events already too big, then set events and censoring = 0 on all further time intervals */
                    if (sumdL >= tot_events) then do;
                        d[lower[n_int]:upper[n_int]] = j(1,(upper[n_int] - lower[n_int] + 1),0);
                        cen[lower[n_int]:(upper[n_int])] = j(1,(upper[n_int] - lower[n_int] + 1),0);
                        n_hat[(lower[n_int]+1):(upper[n_int]+1)] = j(1,(upper[n_int] + 1 - lower[n_int]),n_risk[n_int]);
                    end;
                end;
                /* Otherwise, adjust number censored to give correct total number events */
                if (sumdL < tot_events) | (n_int = 1) then do;
                    sumd = sum(d[1:upper[n_int]]);
                    do while ((sumd > tot_events) | ((sumd < tot_events) & (n_censor[n_int] > 0)));
                        n_censor[n_int] = n_censor[n_int] + (sumd - tot_events);
                        if n_censor[n_int] <= 0 then do;
                            cen[lower[n_int]:(upper[n_int]-1)] = 0;
                            n_censor[n_int] = 0;
                        end;
                        if n_censor[n_int] > 0 then do;
                            cen_t = j(1,n_censor[n_int],0);
                            do j = 1 to n_censor[n_int];
                                cen_t[j] = t_S[lower[n_int]] + j * (t_S[upper[n_int]] - t_S[lower[n_int]]) / (n_censor[n_int] + 1);
                            end;
                            /* Distribute censored observations evenly over time */
                                b = bin(cen_t,t_S[lower[n_int]:upper[n_int]]);
                                call tabulate(levels, freq, b);
                                do f = 1 to ncol(levels);
                                    cen[(lower[n_int]+f-1)] = freq[levels[f]];
                                end;
                        end;
                        n_hat[lower[n_int]] = n_risk[n_int];
                        last = last_i[n_int];
                        do k = lower[n_int] to upper[n_int];
                            d[k] = round(n_hat[k] * (1 - (S[k] / KM_hat[last])));
                            KM_hat[k] = KM_hat[last] * (1 - (d[k] / n_hat[k]));
                            if k ^= upper[n_int] then do;
                                n_hat[k+1] = n_hat[k] - d[k] - cen[k];
                                /* Number at risk cannot be negative */
                                if n_hat[k+1] < 0 then do;
                                    n_hat[k+1] = 0;
                                    cen[k] = n_hat[k] - d[k];
                                end;
                            end;
                            if (d[k] ^= 0) then last = k;
                        end;
                        sumd = sum(d[1:upper[n_int]]);
                    end;
                end;
            end;
			/* note that n_hat is 1,101 but when requesting n_hat[1:n_t] it becomes 101,1 */
            wt1 = t(t_S) || n_hat[1:n_t] || t(d) || t(cen);

        /* 206-234 */
        /* Form IPD */
        /* Initialise vectors */
            t_IPD = j(1,n_risk[1],t_S[n_t]);
            event_IPD = j(1,n_risk[1],0);
        /* Write event time and event indicator (=1) for each event, as separate row in t_IPD and event_IPD */
            k=1;
            do j = 1 to n_t;
                if d[j] ^= 0 then do;
                    t_IPD[k:(k+d[j]-1)] = j(1,abs(d[j]),t_S[j]);
                    event_IPD[k:(k+d[j]-1)] = j(1,abs(d[j]),1);
                    k = k + d[j];
                end;
            end;
        /* Write censor time and event indicator (=0) for each censor, as separate row in t_IPD and event_IPD */
            do j = 1 to (n_t-1);
                if cen[j] ^= 0 then do;
                    t_IPD[k:(k+cen[j]-1)] = j(1,cen[j],((t_S[j]+t_S[j+1])/2));
                    event_IPD[k:(k+cen[j]-1)] = j(1,cen[j],0);
                    k = k + cen[j];
                end;
            end;
            IPD = t(t_IPD) || t(event_IPD);

        return(IPD);

	finish guyot_method;


    /* Example of Use */
        ipd_curve1 = guyot_method(0:100, exp((-1/50)*(0:100)), {0 10 50}, {1000 800 250});
        ipd_curve2 = guyot_method(0:100, exp((-1/30)*(0:100)), {0 10 50}, {1000 700 100});

        colnames={"time" "event" "treatment"};
            treat1 = j(nrow(ipd_curve1),1,'active');
            treat2 = j(nrow(ipd_curve2),1,'chemo');
        create ipd from ipd_curve1 treat1 [colname=colnames];
            append from ipd_curve1 treat1;
            append from ipd_curve2 treat2;
        close ipd;

        submit;
            ods graphics on;
                proc lifetest data=ipd plots=s(cl nocensor);
                    time time * event(0);
                    strata treatment;
                run;
            ods graphics off;
        endsubmit;

quit;



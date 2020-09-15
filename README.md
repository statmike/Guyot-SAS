# Guyot-SAS

---

## Reconstructing individual-level data from published Kaplan-Meier survival curves using the method proposed by Guyot, Patricia et al (2012)

---




Outline:
- The main idea
- Key article cited
- other methods
- overview of methods link
- what has already be done: original, optimized R
- What this is: SAS version following the optimized R version
- IML overview
- DS overview??

### Example
The following code runs the guyot_method module with created input data.  The resulting individual patient data (ipd) is the constructed into a SAS dataset with columns <time, event, treatment>.  The SAS dataset is used to fit a survival model with PROC LIFETEST (SAS/STAT) and the resulting Kaplan-Meier fit can be viewed with the included plot below.

```
PROC IML;
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
QUIT;
```

<p align="center">
  <img src="example.png" />
</p>


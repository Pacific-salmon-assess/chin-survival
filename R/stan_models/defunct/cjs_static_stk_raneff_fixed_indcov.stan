// Stan Code for CJS Model w/ Random Effects 
// Based on swallows example

data {
	int<lower = 2> K; // capture events
	int<lower = 0> I; // number of individuals
	int<lower = 0, upper = 1> CH[I, K]; // CH[i,k]: individual i captured at K
	int<lower = 0> nstock; // number of stocks
	int<lower = 0, upper = nstock> stock[I];// index of group variable
	vector[I] exp_z; // explanatory ind covariate (e.g. size), z-trans.
	//int<lower = 1,upper = 4> year[I]; // year effect (removed)
	//vector[K] agec; // covariate for probability of detection (removed)
}

transformed data {
	int<lower = 0, upper = K + 1> last[I]; // last[i]: ind i last capture
	last = rep_array(0,I);
	for (i in 1:I) {
		for (k in 1:K) {
			if (CH[i, k] == 1) {
				if (k > last[i]) last[i] = k;
			}
		}
	}
}

parameters {
	real a[K - 1]; // intercept of phi
	real a1; // coef of phi
	real<lower = 0> sigmaphi; // between stock sd in logit(phi)
	real<lower = 0> sigmap; // between stock sd in logit(p)
	real stkeffphi[nstock]; // stock effects for phi
	real stkeffp[nstock]; // stock effects for p
}

transformed parameters {
	real<lower = 0, upper = 1> p[I, K]; // capture probability
	real<lower = 0, upper = 1> phi[I, K - 1]; // survival probability
	real<lower = 0, upper = 1> chi[I, K + 1]; // probability that an individual  
									   // is never recap. after its last cap.
	{
		int k;
		for(ii in 1:I) {
			for(tt in 1:(K - 1)) {
				// linear predictor with random effect for phi:
				// add fixed and random effects here
				phi[ii,tt] = inv_logit(a[tt] + (a1 * exp_z[ii]) + 
					(sigmaphi * stkeffphi[stock[ii]]));
			}
		}
		for(i in 1:I) {
			// linear predictor with random effect
			// for p: add fixed and random effects here
			p[i, 1] = 1; // first occasion is marking occasion
			for(kk in 2:K) {
				p[i, kk] = inv_logit(sigmap * stkeffp[stock[i]]);
			}

			// probability that an individual is never recaptured after its
			// last capture
			chi[i, K + 1] = 1.0;
			k = K;
			while (k > 1) {
				chi[i, k] = (1 - phi[i, k-1]) + phi[i, k-1] * (1 - p[i, k]) * 
					chi[i, k+1];
				k = k - 1;
			}
			chi[i,1] = (1 - p[i, 1]) * chi[i, 2];
		}
	}
}

model {
	// priors
	for(v in 1:(K-1)) {
		a[v] ~ normal(0, 5);
	}
	a1 ~ normal(0, 5);
	sigmaphi ~ student_t(2, 0, 1);
	sigmap ~ student_t(2, 0, 1);
	
	// random effects
	for(g in 1:nstock) {
		stkeffphi[g] ~ normal(0, 1);
		stkeffp[g] ~ normal(0,1);
	}
	
	// likelihood
	for (i in 1:I) {
		if (last[i] > 0) {
			for (k in 1:last[i]) {
				if(k > 1) 1 ~ bernoulli(phi[i, k-1]);
				CH[i,k] ~ bernoulli(p[i, k]);
			}
		}
		// deprecated function; replace w/ below 
		// increment_log_prob(log(chi[i, last[i] + 1]));
		target += log(chi[i, last[i] + 1]);
	}
}

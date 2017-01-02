#include <R.h>
#include <Rmath.h>

int convert_2D_indices_to_1D (int i, int j, int* nrow, int* ncol) {
	return(j* (*nrow) + i);
}

/* disp = size parameter for R's negative binomial */

void test (int* x, double* size, double* mu, double* out) {
	GetRNGstate();
	*out = dnbinom_mu(*x, *size, *mu, FALSE);
	PutRNGstate();
}

/* M3Drop */
double calc_prob_m3d (int obs, double mu, double K, double coeff1, double coeff2) {
	if (obs == 0) {
		return(log(1.0 - mu/(mu+K)));
	}
	double cv2 = exp(coeff1+coeff2*log(mu)/log(10));
	double v = cv2*mu*mu;
	if (v <= mu) {v = 1.01*mu;}
	double disp = (mu*mu)/(v-mu);

	GetRNGstate();
	double p = dnbinom_mu(obs, disp, mu, FALSE)
	PutRNGstate();
	if (p < 10^-100) {
		p = 10^-100;
	}
	return(p);	
}

/*double cv2_2_disp = function (double cv2, double mu) {
	double v = cv2*mu*mu;
	if (v <= mu) {v = 1.01*mu;}
	return((mu*mu)/(v-mu));
}*/

void loglikehood_m3d (int* counts, double* mus, int* groups, double* group_factors, double* disps, int* nc, int* ng, int* n_group, double* disp_slope, double* pvalues_out) {
	
	int i, j;
	for (j = 0; j < *ng; j++) {
		double p_null = 0;
		double p_test = 0;
		for (i = 0; i < *nc; i++) {
			int group = groups[i];
			int coords = convert_2D_indices_to_1D(j, i, ng, nc);
			double this_mu = mus[coords];
			double group_mu = this_mu*group_factors[coords];
			int this_obs = counts[coords];
			double this_disp = disps[j];

			/* Shift Disp */
			double new_intercept = log(this_disp)-disp_slope*log(this_mu);
			double new_disp = new_intercept + disp_slope*log(group_mu);
			double shifted_disp = exp(new_disp);

			GetRNGstate();
			double p = dnbinom_mu(this_obs, this_disp, this_mu, TRUE);
			double p2 = dnbinom_mu(this_obs, shifted_disp, group_mu, TRUE);
			PutRNGstate();
			p_null = p_null+p;
			p_test = p_test+p2;
		}
		double D = -2*(p_null-p_test);
		int df = *n_group -1;
		GetRNGstate();
		double p = pchisq(D, df, FALSE, FALSE);
		PutRNGstate();
		pvalues_out[j] = p;
	}
}{

}
/* NBumi */
double calc_prob_nbumi () {

}
/*double shift_disp (double mu, double disp, double mu_i, double slope) {
	double new_intercept = log(disp)-slope*log(mu);
	double new_disp = new_intercept + slope*log(mu_i);
	return(exp(new_disp));
}*/
double var_2_dip () {
}

void loglikehood_nbumi (int* counts, double* mus, int* groups, double* group_factors, double* disps, int* nc, int* ng, int* n_group, double* disp_slope, double* pvalues_out) {
	
	int i, j;
	for (j = 0; j < *ng; j++) {
		double p_null = 0;
		double p_test = 0;
		for (i = 0; i < *nc; i++) {
			int group = groups[i];
			int coords = convert_2D_indices_to_1D(j, i, ng, nc);
			double this_mu = mus[coords];
			double group_mu = this_mu*group_factors[convert_2D_indices_to_1D(j, group, ng, n_group)];
			int this_obs = counts[coords];
			double this_disp = disps[j];

			/* Shift Disp */
			double new_intercept = log(this_disp)-disp_slope*log(this_mu);
			double new_disp = new_intercept + disp_slope*log(group_mu);
			double shifted_disp = exp(new_disp);

			GetRNGstate();
			double p = dnbinom_mu(this_obs, this_disp, this_mu, TRUE);
			double p2 = dnbinom_mu(this_obs, shifted_disp, group_mu, TRUE);
			PutRNGstate();
			p_null = p_null+p;
			p_test = p_test+p2;
		}
		double D = -2*(p_null-p_test);
		int df = *n_group -1;
		GetRNGstate();
		double p = pchisq(D, df, FALSE, FALSE);
		PutRNGstate();
		pvalues_out[j] = p;
	}
}

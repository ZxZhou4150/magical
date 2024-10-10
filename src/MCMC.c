#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>

// 1
void TF_activity_T_sampling(double **A, double **A_sample, int *ATAC_Cell_Sample_vector, 
                            double **R, double **R_sample, int *RNA_Cell_Sample_vector, 
                            double **T_sample, double **T_A, double **T_R, 
                            double *T_prior_mean, double **T_prior_var, 
                            double **B, int **B_state, double sigma_A_noise, 
                            int P, int G, int M, int S) {
    int *TF_index = (int*)malloc(M * sizeof(int));
    
    // Randomly permute indices (shuffling)
    for (int i = 0; i < M; i++) {
        TF_index[i] = i;
    }
    for (int i = M - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = TF_index[i];
        TF_index[i] = TF_index[j];
        TF_index[j] = temp;
    }

    for (int i = 0; i < M; i++) {
        int m = TF_index[i];
        int B_sum = 0;
        for (int j = 0; j < P; j++) {
            B_sum += B_state[j][m];
        }

        if (B_sum > 0) {
            double temp_var = 0.0;
            for (int j = 0; j < P; j++) {
                temp_var += B[j][m] * B[j][m];
            }
            temp_var = temp_var * T_prior_var[m][0] / B_sum + sigma_A_noise;

            double mean_T[S];
            double variance_T[S];

            for (int s = 0; s < S; s++) {
                mean_T[s] = 0.0;
                for (int j = 0; j < P; j++) {
                    mean_T[s] += B[j][m] * (A_sample[j][s] - B[j][s] * T_sample[m][s] + B[j][m] * T_sample[m][s]);
                }
                mean_T[s] = (mean_T[s] / B_sum + T_prior_mean[m] * sigma_A_noise) / temp_var;
                variance_T[s] = T_prior_var[m][s] * sigma_A_noise / temp_var;
            }

            for (int s = 0; s < S; s++) {
                double aa = rnorm(0.0, 1.0);
                if (aa - 3 > 0) {
                    aa = 3;
                }
                if (aa + 3 < 0) {
                    aa = -3;
                }

                T_sample[m][s] = aa * sqrt(fabs(variance_T[s])) + mean_T[s];

                for (int idx = 0; idx < P; idx++) {
                    if (ATAC_Cell_Sample_vector[idx] == s) {
                        T_A[m][idx] = rnorm(0.0, 1.0) * sqrt(fabs(variance_T[s])) + T_sample[m][s];
                    }
                    if (RNA_Cell_Sample_vector[idx] == s) {
                        T_R[m][idx] = rnorm(0.0, 1.0) * sqrt(fabs(variance_T[s])) + T_sample[m][s];
                    }
                }
            }
        }
    }
    free(TF_index);
}

// 2
void TF_peak_binding_B_sampling(double **A, double **A_sample, int *ATAC_Cell_Sample_vector, 
                                double **T_sample, double **T_A, double **TFA_T_sample, 
                                double **B, int **B_state, double *B_prior_mean, 
                                double *B_prior_var, double sigma_A_noise, 
                                int P, int G, int M, int S) {
    double **T_sample_matrix = (double **)malloc(M * sizeof(double *));
    for (int i = 0; i < M; i++) {
        T_sample_matrix[i] = (double *)malloc(S * sizeof(double));
    }

    // Filling T_sample
    for (int s = 0; s < S; s++) {
        int count = 0;
        double *totake_sum = (double *)calloc(M, sizeof(double));

        for (int p = 0; p < P; p++) {
            if (ATAC_Cell_Sample_vector[p] == s) {
                count++;
                for (int m = 0; m < M; m++) {
                    totake_sum[m] += T_A[m][p];
                }
            }
        }

        for (int m = 0; m < M; m++) {
            if (count == 1) {
                T_sample_matrix[m][s] = TFA_T_sample[m][s];
            } else if (count > 1) {
                T_sample_matrix[m][s] = totake_sum[m] / count;
            }
        }
        free(totake_sum);
    }

    // Sampling TF indices
    int *TF_index = (int *)malloc(M * sizeof(int));
    for (int i = 0; i < M; i++) {
        TF_index[i] = i;
    }
    for (int i = M - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = TF_index[i];
        TF_index[i] = TF_index[j];
        TF_index[j] = temp;
    }

    // Main loop
    for (int i = 0; i < M; i++) {
        int m = TF_index[i];
        double temp_var = 0.0;
        double mean_B[P];
        double variance_B = 0.0;

        // Calculate temp_var
        for (int s = 0; s < S; s++) {
            temp_var += TFA_T_sample[m][s] * TFA_T_sample[m][s];
        }
        temp_var = temp_var * B_prior_var[m] / S + sigma_A_noise;

        // Calculate mean_B
        for (int p = 0; p < P; p++) {
            double dot_product = 0.0;
            for (int s = 0; s < S; s++) {
                dot_product += (A_sample[p][s] - B[p][m] * TFA_T_sample[m][s]) * TFA_T_sample[m][s];
            }
            mean_B[p] = (dot_product * B_prior_var[m] / S + B_prior_mean[p] * sigma_A_noise) / temp_var;
        }

        variance_B = B_prior_var[m] * sigma_A_noise / temp_var;

        double bb[P];
        for (int p = 0; p < P; p++) {
            bb[p] = rnorm(0.0, 1.0);
            if (bb[p] - 3 > 0) bb[p] = 3;
            if (bb[p] + 3 < 0) bb[p] = -3;

            B[p][m] = (bb[p] * sqrt(fabs(variance_B)) + mean_B[p]) * B_state[p][m];
        }
    }

    // Clean up
    for (int i = 0; i < M; i++) {
        free(T_sample_matrix[i]);
    }
    free(T_sample_matrix);
    free(TF_index);
}

// 3
void Peak_gene_looping_L_sampling(double **R, double **R_sample, int *RNA_Cell_Sample_vector, 
                                  double **TFA_T_sample, double **TFA_T_R, 
                                  double **B, double **L, int **L_state, 
                                  double *L_prior_mean, double L_prior_var, 
                                  double sigma_R_noise, int P, int G, int M, int S) {
    double **T_sample = (double **)malloc(M * sizeof(double *));
    for (int i = 0; i < M; i++) {
        T_sample[i] = (double *)malloc(S * sizeof(double));
    }

    // Calculate T_sample
    for (int s = 0; s < S; s++) {
        int count = 0;
        double *totake_sum = (double *)calloc(M, sizeof(double));

        for (int p = 0; p < P; p++) {
            if (RNA_Cell_Sample_vector[p] == s) {
                count++;
                for (int m = 0; m < M; m++) {
                    totake_sum[m] += TFA_T_R[m][p];
                }
            }
        }

        for (int m = 0; m < M; m++) {
            if (count == 1) {
                T_sample[m][s] = TFA_T_sample[m][s];
            } else if (count > 1) {
                T_sample[m][s] = totake_sum[m] / count;
            }
        }
        free(totake_sum);
    }

    // Estimate A
    double **A_estimate = (double **)malloc(P * sizeof(double *));
    for (int i = 0; i < P; i++) {
        A_estimate[i] = (double *)malloc(S * sizeof(double));
        for (int s = 0; s < S; s++) {
            A_estimate[i][s] = 0.0;
            for (int m = 0; m < M; m++) {
                A_estimate[i][s] += B[i][m] * TFA_T_sample[m][s];
            }
        }
    }

    // Randomly permute indices for peaks
    int *Peak_index = (int *)malloc(P * sizeof(int));
    for (int i = 0; i < P; i++) {
        Peak_index[i] = i;
    }
    for (int i = P - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = Peak_index[i];
        Peak_index[i] = Peak_index[j];
        Peak_index[j] = temp;
    }

    // Main loop
    for (int i = 0; i < P; i++) {
        int f = Peak_index[i];
        double temp_var = 0.0;
        double mean_L[G];
        double variance_L = 0.0;

        // Calculate temp_var
        for (int s = 0; s < S; s++) {
            temp_var += A_estimate[f][s] * A_estimate[f][s];
        }
        temp_var = temp_var * L_prior_var / S + sigma_R_noise;

        // Calculate mean_L
        for (int g = 0; g < G; g++) {
            double dot_product = 0.0;
            for (int s = 0; s < S; s++) {
                dot_product += A_estimate[f][s] * (R_sample[g][s] - L[f][g] * A_estimate[f][s]);
            }
            mean_L[g] = (dot_product * L_prior_var / S + L_prior_mean[f] * sigma_R_noise) / temp_var;
        }

        variance_L = L_prior_var * sigma_R_noise / temp_var;

        double ll[G];
        for (int g = 0; g < G; g++) {
            ll[g] = rnorm(0.0, 1.0);
            if (ll[g] - 3 > 0) ll[g] = 3;
            if (ll[g] + 3 < 0) ll[g] = -3;

            L[f][g] = (ll[g] * sqrt(fabs(variance_L)) + mean_L[g]) * L_state[f][g];
        }
    }

    // Clean up
    for (int i = 0; i < M; i++) {
        free(T_sample[i]);
    }
    free(T_sample);

    for (int i = 0; i < P; i++) {
        free(A_estimate[i]);
    }
    free(A_estimate);
    free(Peak_index);
}

// 4
void TF_peak_binary_binding_B_state_sampling(double **A, double **A_sample, int *ATAC_Cell_Sample_vector, 
                                             double **TFA_T_A, double **TFA_T_sample, 
                                             double **B, int **B_state, double **B_prior_mean, 
                                             double *B_prior_var, double **B_prior_prob, 
                                             double sigma_A_noise, int P, int G, int M, int S) {
    double **T_sample = (double **)malloc(M * sizeof(double *));
    for (int i = 0; i < M; i++) {
        T_sample[i] = (double *)malloc(S * sizeof(double));
    }

    // Calculate T_sample
    for (int s = 0; s < S; s++) {
        int count = 0;
        double *totake_sum = (double *)calloc(M, sizeof(double));

        for (int p = 0; p < P; p++) {
            if (ATAC_Cell_Sample_vector[p] == s) {
                count++;
                for (int m = 0; m < M; m++) {
                    totake_sum[m] += TFA_T_A[m][p];
                }
            }
        }

        for (int m = 0; m < M; m++) {
            if (count == 1) {
                T_sample[m][s] = TFA_T_A[m][s];
            } else if (count > 1) {
                T_sample[m][s] = totake_sum[m] / count;
            }
        }
        free(totake_sum);
    }

    // Sampling TF and Peak indices
    int *TF_index = (int *)malloc(M * sizeof(int));
    for (int i = 0; i < M; i++) {
        TF_index[i] = i;
    }
    for (int i = M - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = TF_index[i];
        TF_index[i] = TF_index[j];
        TF_index[j] = temp;
    }

    int *Peak_index = (int *)malloc(P * sizeof(int));
    for (int i = 0; i < P; i++) {
        Peak_index[i] = i;
    }
    for (int i = P - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = Peak_index[i];
        Peak_index[i] = Peak_index[j];
        Peak_index[j] = temp;
    }

    // Main loop
    for (int i = 0; i < P; i++) {
        int f = Peak_index[i];
        double *temp = (double *)calloc(S, sizeof(double));

        for (int s = 0; s < S; s++) {
            for (int m = 0; m < M; m++) {
                temp[s] += B[f][m] * TFA_T_sample[m][s];
            }
        }

        for (int j = 0; j < M; j++) {
            int m = TF_index[j];
            double temp_var = 0.0;
            double mean_B = 0.0;
            double variance_B = 0.0;
            double P1 = 0.0;
            double threshold_c = 0.0;

            for (int s = 0; s < S; s++) {
                temp_var += TFA_T_sample[m][s] * TFA_T_sample[m][s];
            }
            temp_var = temp_var * B_prior_var[m] / S + sigma_A_noise;

            if (B_state[f][m] > 0) {
                for (int s = 0; s < S; s++) {
                    mean_B += (A_sample[f][s] - temp[s] + B[f][m] * TFA_T_sample[m][s]) * TFA_T_sample[m][s];
                }
                mean_B = (mean_B * B_prior_var[m] / S + B_prior_mean[f][m] * sigma_A_noise) / temp_var;

                variance_B = B_prior_var[m] * sigma_A_noise / temp_var;

                double post_b1 = exp(-(B[f][m] - mean_B) * (B[f][m] - mean_B) / (2 * variance_B)) * (B_prior_prob[f][m] + 0.25) + 1e-6;
                double post_b0 = exp(-mean_B * mean_B / (2 * variance_B)) * (1 - B_prior_prob[f][m] + 0.25) + 1e-6;

                P1 = post_b1 / (post_b1 + post_b0);
                threshold_c = runif(0.0, 1.0);

                if (isnan(P1) || isnan(P1)) {
                    P1 = 0.5;
                }

                if (P1 < threshold_c) {
                    B[f][m] = 0;
                    B_state[f][m] = 0;
                } else {
                    B[f][m] = B[f][m];
                    B_state[f][m] = 1;
                }
            }

            if (B_state[f][m] == 0 && B_prior_prob[f][m] > 0) {
                for (int s = 0; s < S; s++) {
                    mean_B += (A_sample[f][s] - temp[s]) * TFA_T_sample[m][s];
                }
                mean_B = (mean_B * B_prior_var[m] / S + B_prior_mean[f][m] * sigma_A_noise) / temp_var;

                variance_B = B_prior_var[m] * sigma_A_noise / temp_var;

                double bb = rnorm(0.0, 1.0);
                if (bb - 3 > 0) bb = 3;
                if (bb + 3 < 0) bb = -3;

                double B_temp = bb * sqrt(fabs(variance_B)) + mean_B;

                double post_b1 = exp(-(bb * bb) / 2) * (B_prior_prob[f][m] + 0.25) + 1e-6;
                double post_b0 = exp(-mean_B * mean_B / (2 * variance_B)) * (1 - B_prior_prob[f][m] + 0.25) + 1e-6;

                P1 = post_b1 / (post_b1 + post_b0);
                threshold_c = runif(0.0, 1.0);

                if (isnan(P1) || isnan(P1)) {
                    P1 = 0.5;
                }

                if (P1 < threshold_c) {
                    B[f][m] = 0;
                    B_state[f][m] = 0;
                } else {
                    B[f][m] = B_temp;
                    B_state[f][m] = 1;
                }
            }
        }
        free(temp);
    }

    // Clean up
    for (int i = 0; i < M; i++) {
        free(T_sample[i]);
    }
    free(T_sample);
    free(TF_index);
    free(Peak_index);
}

// 5
void Peak_gene_binary_looping_L_state_sampling(double **R, double **R_sample, int *RNA_Cell_Sample_vector, 
                                               double **TFA_T_R, double **TFA_T_sample, 
                                               double **B, double **L, int **L_state, 
                                               double **L_prior_mean, double L_prior_var, 
                                               double **L_prior_prob, double sigma_R_noise, 
                                               int P, int G, int M, int S) {

    // Permute indices
    int *Peak_index = (int *)malloc(P * sizeof(int));
    for (int i = 0; i < P; i++) {
        Peak_index[i] = i;
    }
    for (int i = P - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = Peak_index[i];
        Peak_index[i] = Peak_index[j];
        Peak_index[j] = temp;
    }

    int *Gene_index = (int *)malloc(G * sizeof(int));
    for (int i = 0; i < G; i++) {
        Gene_index[i] = i;
    }
    for (int i = G - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = Gene_index[i];
        Gene_index[i] = Gene_index[j];
        Gene_index[j] = temp;
    }

    // Calculate T_sample
    double **T_sample = (double **)malloc(M * sizeof(double *));
    for (int i = 0; i < M; i++) {
        T_sample[i] = (double *)malloc(S * sizeof(double));
    }

    for (int s = 0; s < S; s++) {
        int count = 0;
        double *totake_sum = (double *)calloc(M, sizeof(double));

        for (int g = 0; g < G; g++) {
            if (RNA_Cell_Sample_vector[g] == s) {
                count++;
                for (int m = 0; m < M; m++) {
                    totake_sum[m] += TFA_T_R[m][g];
                }
            }
        }

        for (int m = 0; m < M; m++) {
            if (count == 1) {
                T_sample[m][s] = TFA_T_R[m][s];
            } else if (count > 1) {
                T_sample[m][s] = totake_sum[m] / count;
            }
        }
        free(totake_sum);
    }

    // Calculate A_estimate
    double **A_estimate = (double **)malloc(P * sizeof(double *));
    for (int i = 0; i < P; i++) {
        A_estimate[i] = (double *)malloc(S * sizeof(double));
        for (int s = 0; s < S; s++) {
            A_estimate[i][s] = 0.0;
            for (int m = 0; m < M; m++) {
                A_estimate[i][s] += B[i][m] * TFA_T_sample[m][s];
            }
        }
    }

    // Main loop
    for (int i = 0; i < G; i++) {
        int g = Gene_index[i];
        double *temp = (double *)calloc(S, sizeof(double));

        for (int s = 0; s < S; s++) {
            for (int f = 0; f < P; f++) {
                temp[s] += L[f][g] * A_estimate[f][s];
            }
        }

        for (int j = 0; j < P; j++) {
            int f = Peak_index[j];
            double temp_var = 0.0;
            double mean_L = 0.0;
            double variance_L = 0.0;
            double P1 = 0.0;
            double threshold_c = 0.0;

            for (int s = 0; s < S; s++) {
                temp_var += A_estimate[f][s] * A_estimate[f][s];
            }
            temp_var = temp_var * L_prior_var / S + sigma_R_noise;

            if (L_state[f][g] > 0) {
                for (int s = 0; s < S; s++) {
                    mean_L += (R_sample[g][s] - temp[s] + L[f][g] * A_estimate[f][s]) * A_estimate[f][s];
                }
                mean_L = (mean_L * L_prior_var / S + L_prior_mean[f][g] * sigma_R_noise) / temp_var;

                variance_L = L_prior_var * sigma_R_noise / temp_var;

                double post_l1 = exp(-(L[f][g] - mean_L) * (L[f][g] - mean_L) / (2 * variance_L)) * (L_prior_prob[f][g] + 0.25) + 1e-6;
                double post_l0 = exp(-mean_L * mean_L / (2 * variance_L)) * (1 - L_prior_prob[f][g] + 0.25) + 1e-6;

                P1 = post_l1 / (post_l1 + post_l0);
                threshold_c = runif(0.0, 1.0);

                if (isnan(P1) || isnan(P1)) {
                    P1 = 0.5;
                }

                if (P1 < threshold_c) {
                    L[f][g] = 0;
                    L_state[f][g] = 0;
                } else {
                    L[f][g] = L[f][g];
                    L_state[f][g] = 1;
                }
            }

            if (L_state[f][g] == 0 && L_prior_prob[f][g] > 0) {
                for (int s = 0; s < S; s++) {
                    mean_L += (R_sample[g][s] - temp[s]) * A_estimate[f][s];
                }
                mean_L = (mean_L * L_prior_var / S + L_prior_mean[f][g] * sigma_R_noise) / temp_var;

                variance_L = L_prior_var * sigma_R_noise / temp_var;

                double ll = rnorm(0.0, 1.0);
                if (ll - 3 > 0) ll = 3;
                if (ll + 3 < 0) ll = -3;

                double L_temp = ll * sqrt(fabs(variance_L)) + mean_L;

                double post_l1 = exp(-ll * ll / 2) * (L_prior_prob[f][g] + 0.1) + 1e-6;
                double post_l0 = exp(-mean_L * mean_L / (2 * variance_L)) * (1 - L_prior_prob[f][g] + 0.1) + 1e-6;

                P1 = post_l1 / (post_l1 + post_l0);
                threshold_c = runif(0.0, 1.0);

                if (isnan(P1) || isnan(P1)) {
                    P1 = 0.5;
                }

                if (P1 < threshold_c) {
                    L[f][g] = 0;
                    L_state[f][g] = 0;
                } else {
                    L[f][g] = L_temp;
                    L_state[f][g] = 1;
                }
            }
        }
        free(temp);
    }

    // Clean up
    for (int i = 0; i < M; i++) {
        free(T_sample[i]);
    }
    free(T_sample);

    for (int i = 0; i < P; i++) {
        free(A_estimate[i]);
    }
    free(A_estimate);

    free(Peak_index);
    free(Gene_index);
}

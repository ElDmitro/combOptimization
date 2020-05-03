#include "glpk.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void readWarehouseStats(int *wSize, int *cSize,
                        double **capacity, double **openCost,
                        double **demand, double ***useCost) {
    int sstatus;
    sstatus = scanf("%d %d", wSize, cSize);

    *capacity = (double *) calloc(*wSize, sizeof(double));
    *openCost = (double *) calloc(*wSize, sizeof(double));
    for (int i = 0; i < *wSize; ++i) {
        sstatus = scanf("%lf %lf", *capacity + i, *openCost + i);
    }

    *demand = (double *) calloc(*cSize, sizeof(double));
    for (int i = 0; i < *cSize; ++i) {
        sstatus = scanf("%lf", *demand + i);
    }

    double tmp;
    *useCost = (double **) calloc(*wSize, sizeof(double *));
    for (int i = 0; i < *wSize; ++i) {
        (*useCost)[i] = (double *) calloc(*cSize, sizeof(double));

        for (int j = 0; j < *cSize; ++j) {
            sstatus = scanf("%lf", &tmp);
            (*useCost)[i][j] = tmp;
        }
    }
}


void readWarehouseStatsFile(int *wSize, int *cSize,
                        double **capacity, double **openCost,
                        double **demand, double ***useCost,
                        FILE *file) {
    int sstatus;
    sstatus = fscanf(file, "%d %d", wSize, cSize);

    *capacity = (double *) calloc(*wSize, sizeof(double));
    *openCost = (double *) calloc(*wSize, sizeof(double));
    for (int i = 0; i < *wSize; ++i) {
        sstatus = fscanf(file, "%lf %lf", *capacity + i, *openCost + i);
    }

    *demand = (double *) calloc(*cSize, sizeof(double));
    for (int i = 0; i < *cSize; ++i) {
        sstatus = fscanf(file, "%lf", *demand + i);
    }

    double tmp;
    *useCost = (double **) calloc(*wSize, sizeof(double *));
    for (int i = 0; i < *wSize; ++i) {
        (*useCost)[i] = (double *) calloc(*cSize, sizeof(double));

        for (int j = 0; j < *cSize; ++j) {
            sstatus = fscanf(file, "%lf", &tmp);
            (*useCost)[i][j] = tmp;
        }
    }
}


void buildWarehouseProb(glp_prob *problem, int wSize, int cSize,
                        double *capacity, double *openCost,
                        double *demand, double **useCost) {
    // X(w, c) addition
    long long k;
    for (int i = 0; i < wSize; ++i) {
        for (int j = 0; j < cSize; ++j) {
            k = glp_add_cols(problem, 1);

            glp_set_col_name(problem, k, "");
            glp_set_col_kind(problem, k, GLP_CV);
            glp_set_col_bnds(problem, k, GLP_LO, 0., 0.);
            glp_set_obj_coef(problem, k, useCost[i][j]);
        }
    }

    // Y(w) addition
    for (int i = 0; i < wSize; ++i) {
        k = glp_add_cols(problem, 1);

        glp_set_col_name(problem, k, "");
        glp_set_col_kind(problem, k, GLP_BV);
        glp_set_obj_coef(problem, k, openCost[i]);
    }


    // Clients constraints addition
    int *ind = (int *) calloc(1 + wSize, sizeof(int));
    double *val = (double *) calloc(1 + wSize, sizeof(double));
    for (int j = 0; j < cSize; ++j) {
        k = glp_add_rows(problem, 1);

        glp_set_row_name(problem, k, "");
        glp_set_row_bnds(problem, k, GLP_FX, 1., 1.);

        for (int i = 0; i < wSize; ++i) {
            ind[i + 1] = i * cSize + j + 1;
            val[i + 1] = 1.;
        }

        glp_set_mat_row(problem, k, wSize, ind, val);
    }

    // Warehouses constraints addition
    free(ind);
    free(val);
    ind = (int *) calloc(1 + cSize + 1, sizeof(int));
    val = (double *) calloc(1 + cSize + 1, sizeof(double));
    for (int i = 0; i < wSize; ++i) {
        k = glp_add_rows(problem, 1);

        glp_set_row_name(problem, k, "");
        glp_set_row_bnds(problem, k, GLP_UP, 0., 0.);

        for (int j = 0; j < cSize; ++j) {
            ind[j + 1] = i * cSize + j + 1;
            val[j + 1] = demand[j];
        }

        ind[cSize + 1] = wSize * cSize + i + 1;
        val[cSize + 1] = -capacity[i];

        glp_set_mat_row(problem, k, cSize + 1, ind, val);
    }

    free(ind);
    free(val);
}


void freeStoredData(int wSize, int cSize,
                    double *capacity, double *openCost,
                    double *demand, double **useCost) {
    free(capacity);
    free(openCost);
    free(demand);
    for (int i = 0; i < wSize; ++i) {
        free(useCost[i]);
    }
    free(useCost);
}

void freeWarehousePlan(int wSize, int cSize,
                       double **X, int *Y) {
    for (int i = 0; i < wSize; ++i) {
        free(X[i]);
    }
    free(X);
    free(Y);
}


int buildWarehousePlan(glp_prob *problem,
                       int wSize, int cSize,
                       double ***X, int **Y) {
    int status = glp_mip_status(problem);
    if ((status != GLP_FEAS) && (status != GLP_OPT)) {
        return 1;
    }

    *X = (double **) calloc(wSize, sizeof(double *));
    *Y = (int *) calloc(wSize, sizeof(int));
    for (int i = 0; i < wSize; ++i) {
        (*X)[i] = (double *) calloc(cSize, sizeof(double));
        for (int j = 0; j < cSize; ++j) {
            (*X)[i][j] = glp_mip_col_val(problem, 1 + i * cSize + j);
        }

        (*Y)[i] = glp_mip_col_val(problem, 1 + wSize * cSize + i);
    }

    return 0;
}


void showWarehousePlan(int wSize, int cSize, double **X, int *Y) {
    int k = 0;
    for (int i = 0; i < wSize; ++i) {
        if (Y[i] == 0) {
            continue;
        }

        ++k;
    }

    printf("%d\n", k);
    for (int i = 0; i < wSize; ++i) {
        if (Y[i] == 0) {
            continue;
        }

        printf("%d ", i + 1);
    }
    printf("\n");

    for (int i = 0; i < wSize; ++i) {
        if (Y[i] == 0) {
            continue;
        }
        for (int j = 0; j < cSize; ++j) {
            printf("%lf ", X[i][j]);
        }
        printf("\n");
    }
}


int main(int argc, char *argv[]) {
    int wSize, cSize;
    double *capacity, *openCost;
    double *demand;
    double **useCost;

    if (argc > 1) {
        FILE *input;
        input = fopen(argv[1], "r");

        readWarehouseStatsFile(&wSize, &cSize, &capacity, &openCost, &demand, &useCost, input);
    } else {
        readWarehouseStats(&wSize, &cSize, &capacity, &openCost, &demand, &useCost);
    }

    glp_prob *problem = glp_create_prob();
    glp_iocp iocp;

    buildWarehouseProb(problem, wSize, cSize, capacity, openCost, demand, useCost);
    freeStoredData(wSize, cSize, capacity, openCost, demand, useCost);

    // MIP PARAMS INIT
    glp_init_iocp(&iocp);
    iocp.msg_lev = GLP_MSG_OFF;
    iocp.br_tech = GLP_BR_MFV; /* most fractional variable */
    iocp.bt_tech = GLP_BT_BLB; /* best local bound */
    iocp.sr_heur = GLP_OFF; /* disable simple rounding heuristic */
    iocp.gmi_cuts = GLP_ON; /* enable Gomory cuts */
    iocp.presolve = GLP_ON; /* LP-relaxation pre-solution */

    glp_intopt(problem, &iocp);

    double **X;
    int *Y;
    if (buildWarehousePlan(problem, wSize, cSize, &X, &Y)) {
        return 1;
    }
    showWarehousePlan(wSize, cSize, X, Y);
    if (argc > 2) {
       printf("%lf\n", glp_mip_obj_val(problem)); 
    }
    freeWarehousePlan(wSize, cSize, X, Y);
    glp_delete_prob(problem);
}

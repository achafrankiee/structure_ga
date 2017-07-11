#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FALSE 0
#define TRUE  1
#define DATAFILE "data.txt"
#define DATALIENFILE "datalien.txt"

/// Déclarations des fonctions
void init_ga();
int nbLignes(FILE*, char []);
int batch_memory_allocation();
int allocate_memory_float(float***, int, int);
int allocate_memory_int(int***, int, int);
int matrix_fill();
void A_fill();
void reg_mat_compute(float*, int, int);
void K_assembly(float**, float**, int, int);

/// Déclarations des variables globales
int N, Nc, datacols = 4, dataliencols = 2, Acols = 3;
float** data;
int** datalien;
int** conMat;
float** A;
int lCount;
FILE* f = NULL;



int main()
{
    printf("************************** GA **************************\n");

    /// ------------------------ Initialisation ------------------------
    init_ga();


    /// ----------- Calcul de la matrice de régidité globale -----------
    float S[3] = {1,1,1};
    reg_mat_compute(S, N, lCount);


    return 0;
}

void init_ga() {
    printf("\nPress ENTER key to Continue\n");
    getchar();
    ///calcul de nb de lignes
    printf("\n1. Calcul nombre de lignes :");
    N = nbLignes(f, DATAFILE);
    Nc = nbLignes(f, DATALIENFILE);
    if (N == -1 || Nc == -1) {
        printf("\n\t[ERREUR][ECHEC_OUV_FICH]");
        exit(1);
    }
    printf("\n\t---data lignes : %d\n\t---datalien lignes : %d\n", N, Nc);
    printf("\n\t[SUCCESS]\n");

    printf("\n--------------------------------------------------------\n");

    printf("\nPress ENTER key to Continue\n");
    getchar();
    ///Allocation dynamique des matrices
    printf("\n2. allocation de memoire :");
    int allocation_res = batch_memory_allocation();
    if (allocation_res == FALSE) {
        printf("\n\t[ERREUR][MEM_INSUF]");
        exit(1);
    }
    printf("\n\t[SUCCESS]\n");

    printf("\n--------------------------------------------------------\n");

    printf("\nPress ENTER key to Continue\n");
    getchar();
    ///Remplissage des matrices à partir des fichiers
    printf("\n3. Remplissage des matrices a partir des fichiers :");
    int fill_res = matrix_fill();
    if (fill_res == FALSE) {
        printf("\n\t[ERREUR][ECHEC_OUV_FICH]");
        exit(1);
    }
    printf("\n\t[SUCCESS]\n");

    printf("\n--------------------------------------------------------\n");

    printf("\nPress ENTER key to Continue\n");
    getchar();
    ///Affichage
    printf("\n4. Affichage des matrices :");
    printf("\n\tdata\n\t\t");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < datacols; j++) {
            printf("%.2f\t", data[i][j]);
        }
        printf("\n\t\t");
    }
    printf("\n\tdatalien\n\t\t");
    for (int i = 0; i < Nc; i++) {
        for (int j = 0; j < dataliencols; j++) {
            printf("%d\t", datalien[i][j]);
        }
        printf("\n\t\t");
    }
    printf("\n\t[SUCCESS]\n");

    printf("\n--------------------------------------------------------\n");

    printf("\nPress ENTER key to Continue\n");
    getchar();
    ///Calcul de la matrice de connectivité
    printf("\n5. Calcul de la matrice de connectivite :");
    printf("\n\t5.1. Mise a zero :");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            conMat[i][j] = FALSE;
        }
    }
    printf("\n\t\t[SUCCESS]");
    printf("\n\t5.2. Remplissage :");
    int temp0, temp1;
    for (int i = 0; i < Nc; i++) {
            temp0 = datalien[i][0];
            temp1 = datalien[i][1];
            conMat[temp0][temp1] = TRUE;
    }
    printf("\n\t\t[SUCCESS]");
    printf("\n\t5.3. Aaffichage :\n\t\t");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d\t", conMat[i][j]);
        }
        printf("\n\t\t");
    }
    printf("\n\t\t[SUCCESS]");
    printf("\n\t[SUCCESS]\n");

    printf("\n--------------------------------------------------------\n");

    printf("\nPress ENTER key to Continue\n");
    getchar();
    ///Calcul du nombre de liaisons
    printf("\n6. Calcul du nombre de liaisons :");
    lCount = 0;
    for (int i = 0; i < N-1; i++) {
        for (int j = i+1; j < N; j++) {
            if (conMat[i][j] == TRUE) lCount++;
        }
    }
    printf("\n\tNombre de liasons : %d", lCount);
    printf("\n\t[SUCCESS]\n");

    printf("\nPress ENTER key to Continue\n");
    getchar();
    ///Création de la matrice A
    printf("\n7. Creation de la matrice A :");
    printf("\n\t7.1 Allocation de la memoire :");
    allocate_memory_float(&A, lCount, Acols);
    printf("\n\t\t[SUCCESS]\n");
    printf("\n\t7.2 Remplissage :");
    A_fill();
    printf("\n\t\t[SUCCESS]\n");
    printf("\n\t7.3 Affichage :");
    printf("\n\t\tA\n\t\t\t");
    for (int i = 0; i < lCount; i++) {
        for (int j = 0; j < Acols; j++) {
            printf("%.2f\t", A[i][j]);
        }
        printf("\n\t\t\t");
    }
    printf("\n\t\t[SUCCESS]\n");


}

int nbLignes(FILE* f, char fileName []) {
    f = fopen(fileName, "r");
    if (f != NULL) {
        char ch;
        int N = 1; //nb de lignes
        do {
            ch = fgetc(f);
            if(ch == '\n') N++;
        }
        while (ch != EOF);
        fclose(f);
        return(N);
    }
    return(-1);
}

int batch_memory_allocation() {
    int dt = allocate_memory_float(&data, N, datacols);
    int dtl = allocate_memory_int(&datalien, Nc, dataliencols);
    int cnm = allocate_memory_int(&conMat, N, N);
    if (dt == FALSE || dtl == FALSE || cnm == FALSE)
        return FALSE;
    return TRUE;
}

int allocate_memory_float(float*** mat, int nrows, int ncols) {
    *mat = (float**)malloc(nrows * sizeof(float*));
    if (*mat == NULL) return FALSE;
    else {
        for (int i = 0; i < nrows; i++) {
            (*mat)[i] = (float*)malloc(ncols * sizeof(float));
            if ((*mat)[i] == NULL) return FALSE;
        }
    }
    return TRUE;
}

int allocate_memory_int(int*** mat, int nrows, int ncols) {
    *mat = (int**)malloc(nrows * sizeof(int*));
    if (*mat == NULL) return FALSE;
    else {
        for (int i = 0; i < nrows; i++) {
            (*mat)[i] = (int*)malloc(ncols * sizeof(int));
            if ((*mat)[i] == NULL) return FALSE;
        }
    }
    return TRUE;
}

int matrix_fill() {
    f = fopen(DATAFILE, "r");
    if (f != NULL) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < datacols; j++) {
                float x;
                fscanf(f, "%f", &x);
                data[i][j] = x;
            }
        }
    } else return FALSE;
    fclose(f);
    f = fopen(DATALIENFILE, "r");
    if (f != NULL) {
        for (int i = 0; i < Nc; i++) {
            for (int j = 0; j < dataliencols; j++) {
                float x;
                fscanf(f, "%f", &x);
                datalien[i][j] = x;
            }
        }
    } else return FALSE;
    fclose(f);
    return TRUE;
}

void A_fill() {
    int counter = 0;
    float xi;
    float xj;
    float yi;
    float yj;
    for (int i = 0; i < N-1; i++) {
        for (int j = i+1; j < N; j++) {
            if (conMat[i][j] == TRUE) {
                A[counter][0] = i;
                A[counter][1] = j;
                xi = data[i][0];
                xj = data[j][0];
                yi = data[i][1];
                yj = data[j][1];
                A[counter][2] = sqrt(pow(xi - xj, 2) + pow(yi - yj, 2));
                counter++;
            }
        }
    }
}

void reg_mat_compute(float* S, int N, int lCount) {
    float** K;
    float** Ke_block;
    allocate_memory_float(&Ke_block, 2, 2);
    allocate_memory_float(&K, 2*N, 2*N);
    for (int i = 0; i < 2*N; i++) {
        for (int j = 0; j < 2*N; j++) {
            K[i][j] = 0;
        }
    }
    int p, i, j;
    float L, nx, ny;
    for (p = 0; p < lCount; p++) {
        //calcul de nx et ny
        i = (int)(A[p][0]);
        j = (int)(A[p][1]);
        L = A[p][2];
        nx = (data[j][0] - data[i][0]) / L;
        ny = (data[j][1] - data[i][1]) / L;
        printf("\n\nINSIDE REG_MAT_COMPUTE P = %d>>>", p);
        printf("\n\ti:%d\tj:%d\tL:%.1f\tnx:%.1f\tny:%.1f", i, j, L, nx, ny);
        Ke_block[0][0] = nx * nx;
        Ke_block[0][1] = nx * ny;
        Ke_block[1][0] = nx * ny;
        Ke_block[1][1] = ny * ny;
        printf("\n\t%.1f\t%.1f\n\t%.1f\t%.1f", Ke_block[0][0], Ke_block[0][1], Ke_block[1][0], Ke_block[1][1]);
        K_assembly(K, Ke_block, 2*i, 2*i);
        K_assembly(K, Ke_block, 2*i, 2*j);
        K_assembly(K, Ke_block, 2*j, 2*i);
        K_assembly(K, Ke_block, 2*j, 2*j);
        printf("\n\n");
        for (int i = 0; i < 2*N; i++) {
            for (int j = 0; j < 2*N; j++) {
                printf("\t%.1f", K[i][j]);
            }
            printf("\n");
        }
    }

}

void K_assembly(float** K, float** Ke_block, int row, int col) {
    int sign;
    if (row == col)
        sign = 1;
    else
        sign = -1;
    for (int i = row; i <= row + 1; i++) {
        for (int j = col; j <= col + 1; j++) {
            K[i][j] += sign * Ke_block[i%2][j%2];
        }
    }
}

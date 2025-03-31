#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void viderBuffer() {
    int c;
    while ((c = getchar()) != '\n' && c != EOF);
}

// Fonction pour allouer dynamiquement une matrice
double** allouerMatrice(int lignes, int colonnes) {
    double** matrice = malloc(lignes * sizeof(double*));
    if (!matrice) {
        printf("Erreur Allocation de la memoire pour les lignes echouee.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < lignes; i++) {
        matrice[i] = malloc(colonnes * sizeof(double));
        if (!matrice[i]) {
            printf("Erreur Allocation de la memoire pour la ligne %d echouee.\n", i);
            exit(EXIT_FAILURE);
        }
    }
    return matrice;
}

// Fonction pour allouer un vecteur
double* allouerVecteur(int taille) {
    double* vecteur = malloc(taille * sizeof(double));
    if (!vecteur) {
        printf("Erreur: Allocation de la memoire pour le vecteur echouee.\n");
        exit(EXIT_FAILURE);
    }
    return vecteur;
}

// Fonction pour liberer la memoire d'une matrice
void libererMatrice(double** matrice, int lignes) {
    for (int i = 0; i < lignes; i++) {
        free(matrice[i]);
    }
    free(matrice);
}

// Fonction pour afficher une matrice
void afficherMatrice(double** matrice, int lignes, int colonnes) {
    printf("\n");
    for (int i = 0; i < lignes; i++) {
        printf("| ");
        for (int j = 0; j < colonnes; j++) {
            double valeur = matrice[i][j];
            if (fabs(valeur) < 1e-10) {
                valeur = 0.0;
            }
            printf("%8.2f ", valeur);
        }
        printf("|\n");
    }
}

// Fonction pour afficher un vecteur
void afficherVecteur(double* vecteur, int taille) {
    printf("\n");
    for (int i = 0; i < taille; i++) {
        double valeur = vecteur[i];
        if (fabs(valeur) < 1e-10) {
            valeur = 0.0;
        }
        printf("| %8.2f |\n", valeur);
    }
}

// Fonction pour copier une matrice
void copierMatrice(double** source, double** destination, int lignes, int colonnes) {
    for (int i = 0; i < lignes; i++) {
        for (int j = 0; j < colonnes; j++) {
            destination[i][j] = source[i][j];
        }
    }
}

// Produit scalaire de deux vecteurs
double produitScalaire(double* v1, double* v2, int n) {
    double produit = 0.0;
    for (int i = 0; i < n; i++) {
        produit += v1[i] * v2[i];
    }
    return produit;
}

// Norme euclidienne d'un vecteur
double normeVecteur(double* v, int n) {
    return sqrt(produitScalaire(v, v, n));
}

// Decomposition LU (methode de Doolittle)
int decompositionLU(double** A, double** L, double** U, int n) {
    // Initialisation de L et U
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            L[i][j] = 0.0;
            U[i][j] = 0.0;
        }
        L[i][i] = 1.0; // Diagonale de L = 1
    }

    // Algorithme de decomposition LU
    for (int j = 0; j < n; j++) {
        // Calcul des elements de U
        for (int i = 0; i <= j; i++) {
            double somme = 0.0;
            for (int k = 0; k < i; k++) {
                somme += L[i][k] * U[k][j];
            }
            U[i][j] = A[i][j] - somme;
        }

        // Calcul des elements de L
        for (int i = j + 1; i < n; i++) {
            double somme = 0.0;
            for (int k = 0; k < j; k++) {
                somme += L[i][k] * U[k][j];
            }
            if (fabs(U[j][j]) < 1e-10) {
                printf("Erreur: Division par zero detectee. La matrice n'est pas decomposable par LU.\n");
                return 0;
            }
            L[i][j] = (A[i][j] - somme) / U[j][j];
        }
    }
    return 1; // Succes
}

// Resolution de Ly = b (substitution avant)
void substitutionAvant(double** L, double* b, double* y, int n) {
    for (int i = 0; i < n; i++) {
        double somme = 0.0;
        for (int j = 0; j < i; j++) {
            somme += L[i][j] * y[j];
        }
        y[i] = (b[i] - somme) / L[i][i];
    }
}

// Resolution de Ux = y (substitution arriere)
void substitutionArriere(double** U, double* y, double* x, int n) {
    for (int i = n - 1; i >= 0; i--) {
        double somme = 0.0;
        for (int j = i + 1; j < n; j++) {
            somme += U[i][j] * x[j];
        }
        if (fabs(U[i][i]) < 1e-10) {
            printf("Avertissement: Division par une valeur proche de zero.\n");
            x[i] = 0.0;
        } else {
            x[i] = (y[i] - somme) / U[i][i];
        }
    }
}

// Calculer l'erreur residuelle ||Ax - b||
double erreurResiduelle(double** A, double* x, double* b, int n) {
    double* r = allouerVecteur(n);

    // Calculer r = Ax - b
    for (int i = 0; i < n; i++) {
        r[i] = -b[i];
        for (int j = 0; j < n; j++) {
            r[i] += A[i][j] * x[j];
        }
    }

    double norme = normeVecteur(r, n);
    free(r);

    return norme;
}

void effacerEcran() {
#ifdef _WIN32
    system("cls");
#else
    system("clear");
#endif
}

void afficherLigne(int longueur) {
    for (int i = 0; i < longueur; i++) {
        printf("=");
    }
    printf("\n");
}

void afficherTitre(const char* titre) {
    int longueur = strlen(titre) + 10;
    afficherLigne(longueur);
    printf("    %s    \n", titre);
    afficherLigne(longueur);
    printf("\n");
}

int main() {
    int n, choix, choixLU;
    double **A, **L, **U;
    double *b, *x, *y;

    effacerEcran();
    afficherTitre("PROGRAMME DE DECOMPOSITION MATRICIELLE");

    // Saisie de la taille de la matrice
    printf("Entrez la taille de la matrice carree (n) : ");
    if (scanf("%d", &n) != 1 || n <= 0) {
        printf("Erreur: La taille doit etre un entier positif.\n");
        return EXIT_FAILURE;
    }
    viderBuffer();

    // Allocation des matrices
    A = allouerMatrice(n, n);
    L = allouerMatrice(n, n);
    U = allouerMatrice(n, n);

    // Saisie des elements de la matrice
    printf("\nEntrez les elements de la matrice %dx%d:\n", n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("A[%d][%d] = ", i + 1, j + 1);
            if (scanf("%lf", &A[i][j]) != 1) {
                printf("Erreur de lecture. Veuillez reessayer.\n");
                exit(EXIT_FAILURE);
            }
            viderBuffer();
        }
    }

    // Menu pour choisir l'algorithme
    do {
        effacerEcran();
        afficherTitre("DECOMPOSITION MATRICIELLE");

        printf("Matrice d'entree A (%dx%d) :\n", n, n);
        afficherMatrice(A, n, n);

        printf("\nChoisissez l'algorithme à utiliser:\n");
        printf("1. Decomposition LU\n");
        printf("2. Saisir une nouvelle matrice\n");
        printf("0. Quitter\n");
        printf("\nVotre choix : ");
        scanf("%d", &choix);
        viderBuffer();

        switch (choix) {
            case 1: {
                // Decomposition LU
                afficherTitre("RESULTAT DE LA DECOMPOSITION LU");

                int succes = decompositionLU(A, L, U, n);

                if (succes) {
                    printf("La matrice A a ete decomposee en A = L × U où:\n");

                    printf("\nMatrice L (triangulaire inferieure):");
                    afficherMatrice(L, n, n);

                    printf("\nMatrice U (triangulaire superieure):");
                    afficherMatrice(U, n, n);

                    // Verification du resultat (L × U)
                    double **produit = allouerMatrice(n, n);
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                            produit[i][j] = 0.0;
                            for (int k = 0; k < n; k++) {
                                produit[i][j] += L[i][k] * U[k][j];
                            }
                        }
                    }

                    printf("\nVerification (L × U):");
                    afficherMatrice(produit, n, n);

                    libererMatrice(produit, n);

                    // Proposer de resoudre un systeme lineaire
                    printf("\nVoulez-vous resoudre un systeme lineaire Ax = b?\n");
                    printf("1. Oui\n");
                    printf("2. Non\n");
                    printf("\nVotre choix : ");
                    scanf("%d", &choixLU);
                    viderBuffer();

                    if (choixLU == 1) {
                        // Allouer les vecteurs necessaires
                        b = allouerVecteur(n);
                        x = allouerVecteur(n);
                        y = allouerVecteur(n);

                        // Saisir le vecteur b
                        printf("\nEntrez les elements du vecteur b (%dx1):\n", n);
                        for (int i = 0; i < n; i++) {
                            printf("b[%d] = ", i + 1);
                            if (scanf("%lf", &b[i]) != 1) {
                                printf("Erreur de lecture. Veuillez reessayer.\n");
                                exit(EXIT_FAILURE);
                            }
                            viderBuffer();
                        }

                        printf("\nVecteur b:");
                        afficherVecteur(b, n);

                        // Resoudre Ly = b (substitution avant)
                        substitutionAvant(L, b, y, n);

                        printf("\nVecteur intermediaire y (de Ly = b):");
                        afficherVecteur(y, n);

                        // Resoudre Ux = y (substitution arriere)
                        substitutionArriere(U, y, x, n);

                        printf("\nSolution x du systeme Ax = b:");
                        afficherVecteur(x, n);

                        // Calculer l'erreur residuelle
                        double erreur = erreurResiduelle(A, x, b, n);
                        printf("\nErreur residuelle ||Ax - b|| = %e\n", erreur);

                        if (erreur < 1e-8) {
                            printf("La solution est tres precise!\n");
                        } else if (erreur < 1e-3) {
                            printf("La solution est acceptable.\n");
                        } else {
                            printf("Attention: La solution peut manquer de precision.\n");
                        }

                        // Liberer les vecteurs
                        free(b);
                        free(x);
                        free(y);
                    }
                }

                printf("\nAppuyez sur Entree pour continuer...");
                getchar();
                break;
            }
            case 2: {
                // Saisie d'une nouvelle matrice
                effacerEcran();
                afficherTitre("SAISIE D'UNE NOUVELLE MATRICE");

                printf("Entrez la taille de la matrice carree (n) : ");
                if (scanf("%d", &n) != 1 || n <= 0) {
                    printf("Erreur: La taille doit etre un entier positif.\n");
                    printf("\nAppuyez sur Entree pour continuer...");
                    viderBuffer();
                    getchar();
                    break;
                }
                viderBuffer();

                // Liberation et reallocation des matrices
                libererMatrice(A, n);
                libererMatrice(L, n);
                libererMatrice(U, n);

                A = allouerMatrice(n, n);
                L = allouerMatrice(n, n);
                U = allouerMatrice(n, n);

                // Saisie des elements de la nouvelle matrice
                printf("\nEntrez les elements de la matrice %dx%d:\n", n, n);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        printf("A[%d][%d] = ", i + 1, j + 1);
                        if (scanf("%lf", &A[i][j]) != 1) {
                            printf("Erreur de lecture. Veuillez reessayer.\n");
                            exit(EXIT_FAILURE);
                        }
                        viderBuffer();
                    }
                }
                break;
            }
            case 0:
                printf("\nMerci d'avoir utilise ce programme. Au revoir!\n");
                break;
            default:
                printf("\nOption invalide. Veuillez reessayer.\n");
                printf("\nAppuyez sur Entree pour continuer...");
                getchar();
        }
    } while (choix != 0);

    // Liberation de la memoire
    libererMatrice(A, n);
    libererMatrice(L, n);
    libererMatrice(U, n);

    return EXIT_SUCCESS;
}

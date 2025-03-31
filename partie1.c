#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Fonction pour vider le tampon d'entree
void viderBuffer() {
    int c;
    while ((c = getchar()) != '\n' && c != EOF);
}

// Fonction pour allouer dynamiquement une matrice
double** allouerMatrice(int lignes, int colonnes) {
    double** matrice = malloc(lignes * sizeof(double*));
    if (!matrice) {
        fprintf(stderr, "Erreur: Allocation de la memoire pour les lignes echouee.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < lignes; i++) {
        matrice[i] = malloc(colonnes * sizeof(double));
        if (!matrice[i]) {
            fprintf(stderr, "Erreur: Allocation de la memoire pour la ligne %d echouee.\n", i);
            exit(EXIT_FAILURE);
        }
    }
    return matrice;
}

// Fonction pour liberer la memoire d'une matrice
void libererMatrice(double** matrice, int lignes) {
    for (int i = 0; i < lignes; i++) {
        free(matrice[i]);
    }
    free(matrice);
}

// Fonction pour afficher une matrice avec 2 decimales et zero formate correctement
void afficherMatrice(double** matrice, int lignes, int colonnes) {
    printf("\n");
    for (int i = 0; i < lignes; i++) {
        printf("| ");
        for (int j = 0; j < colonnes; j++) {
            double valeur = matrice[i][j];
            if (fabs(valeur) < 1e-10) {  // Si la valeur est proche de zero, la forcer à zero.
                valeur = 0.0;
            }
            printf("%8.2f ", valeur);
        }
        printf("|\n");
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

// Decomposition QR (methode de Gram-Schmidt)
void decompositionQR(double** A, double** Q, double** R, int m, int n) {
    double* u = malloc(m * sizeof(double));
    if (!u) {
        fprintf(stderr, "Erreur: Allocation de memoire pour le vecteur temporaire echouee.\n");
        exit(EXIT_FAILURE);
    }
    
    // Initialiser Q et R à zero
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            Q[i][j] = 0.0;
            R[i][j] = 0.0;
        }
    }
    
    // Methode de Gram-Schmidt modifiee
    for (int j = 0; j < n; j++) {
        // Extraction de la colonne j de A dans u
        for (int i = 0; i < m; i++) {
            u[i] = A[i][j];
        }
        
        // Orthogonalisation
        for (int k = 0; k < j; k++) {
            R[k][j] = 0.0;
            for (int i = 0; i < m; i++) {
                R[k][j] += Q[i][k] * A[i][j];
            }
            for (int i = 0; i < m; i++) {
                u[i] -= R[k][j] * Q[i][k];
            }
        }
        
        // Normalisation
        R[j][j] = normeVecteur(u, m);
        if (fabs(R[j][j]) < 1e-10) {
            printf("Avertissement: Vecteur quasi-nul detecte. La decomposition peut etre imprecise.\n");
            R[j][j] = 1e-10; // Pour eviter la division par zero
        }
        
        for (int i = 0; i < m; i++) {
            Q[i][j] = u[i] / R[j][j];
        }
    }
    
    free(u);
}

// Fonction pour effacer l'ecran (multiplateforme)
void effacerEcran() {
    #ifdef _WIN32
        system("cls");
    #else
        system("clear");
    #endif
}

// Fonction pour afficher une ligne decorative
void afficherLigne(int longueur) {
    for (int i = 0; i < longueur; i++) {
        printf("=");
    }
    printf("\n");
}

// Fonction pour afficher un titre centre
void afficherTitre(const char* titre) {
    int longueur = strlen(titre) + 10;
    afficherLigne(longueur);
    printf("    %s    \n", titre);
    afficherLigne(longueur);
    printf("\n");
}

int main() {
    int n, choix;
    double **A, **L, **U, **Q, **R;
    
    effacerEcran();
    afficherTitre("PROGRAMME DE DeCOMPOSITION MATRICIELLE");
    
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
    Q = allouerMatrice(n, n);
    R = allouerMatrice(n, n);
    
    // Saisie des elements de la matrice
    printf("\nEntrez les elements de la matrice %dx%d:\n", n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("A[%d][%d] = ", i+1, j+1);
            if (scanf("%lf", &A[i][j]) != 1) {
                printf("Erreur de lecture. Veuillez reessayer.\n");
                exit(EXIT_FAILURE);
            }
            viderBuffer();
        }
    }
    
    // Menu principal
    do {
        effacerEcran();
        afficherTitre("DeCOMPOSITION MATRICIELLE");
        
        printf("Matrice d'entree A (%dx%d):\n", n, n);
        afficherMatrice(A, n, n);
        
        printf("\nChoisissez l'algorithme à utiliser :\n");
        printf("1. Decomposition LU\n");
        printf("2. Decomposition QR\n");
        printf("3. Saisir une nouvelle matrice\n");
        printf("0. Quitter\n");
        printf("\nVotre choix : ");
        
        if (scanf("%d", &choix) != 1) {
            printf("Entree invalide.\n");
            viderBuffer();
            continue;
        }
        viderBuffer();
        
        switch (choix) {
            case 1: {
                // Decomposition LU
                effacerEcran();
                afficherTitre("ReSULTAT DE LA DeCOMPOSITION LU");
                
                if (decompositionLU(A, L, U, n)) {
                    printf("La matrice A a ete decomposee en A = L × U où :\n");
                    
                    printf("\nMatrice L (triangulaire inferieure) :");
                    afficherMatrice(L, n, n);
                    
                    printf("\nMatrice U (triangulaire superieure) :");
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
                    
                    printf("\nVerification (L × U) :");
                    afficherMatrice(produit, n, n);
                    libererMatrice(produit, n);
                }
                
                printf("\nAppuyez sur Entree pour continuer...");
                getchar();
                break;
            }
            case 2: {
                // Decomposition QR
                effacerEcran();
                afficherTitre("ReSULTAT DE LA DeCOMPOSITION QR");
                
                decompositionQR(A, Q, R, n, n);
                
                printf("La matrice A a ete decomposee en A = Q × R où :\n");
                
                printf("\nMatrice Q (orthogonale) :");
                afficherMatrice(Q, n, n);
                
                printf("\nMatrice R (triangulaire superieure) :");
                afficherMatrice(R, n, n);
                
                // Verification du resultat (Q × R)
                double **produit = allouerMatrice(n, n);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        produit[i][j] = 0.0;
                        for (int k = 0; k < n; k++) {
                            produit[i][j] += Q[i][k] * R[k][j];
                        }
                    }
                }
                
                printf("\nVerification (Q × R) :");
                afficherMatrice(produit, n, n);
                
                // Verification de l'orthogonalite de Q (Q^T × Q ≈ I)
                printf("\nVerification de l'orthogonalite de Q (Q^T × Q) :\n");
                double **QTQ = allouerMatrice(n, n);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        QTQ[i][j] = 0.0;
                        for (int k = 0; k < n; k++) {
                            QTQ[i][j] += Q[k][i] * Q[k][j];
                        }
                    }
                }
                afficherMatrice(QTQ, n, n);
                
                libererMatrice(produit, n);
                libererMatrice(QTQ, n);
                
                printf("\nAppuyez sur Entree pour continuer...");
                getchar();
                break;
            }
            case 3: {
                // Saisie d'une nouvelle matrice
                effacerEcran();
                afficherTitre("SAISIE D'UNE NOUVELLE MATRICE");
                
                printf("Entrez la taille de la nouvelle matrice carree (n) : ");
                if (scanf("%d", &n) != 1 || n <= 0) {
                    printf("Erreur: La taille doit etre un entier positif.\n");
                    viderBuffer();
                    printf("\nAppuyez sur Entree pour continuer...");
                    getchar();
                    break;
                }
                viderBuffer();
                
                // Liberation et reallocation des matrices
                libererMatrice(A, n);
                libererMatrice(L, n);
                libererMatrice(U, n);
                libererMatrice(Q, n);
                libererMatrice(R, n);
                
                A = allouerMatrice(n, n);
                L = allouerMatrice(n, n);
                U = allouerMatrice(n, n);
                Q = allouerMatrice(n, n);
                R = allouerMatrice(n, n);
                
                // Saisie des elements de la nouvelle matrice
                printf("\nEntrez les elements de la matrice %dx%d:\n", n, n);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        printf("A[%d][%d] = ", i+1, j+1);
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
    libererMatrice(Q, n);
    libererMatrice(R, n);
    
    return EXIT_SUCCESS;
}

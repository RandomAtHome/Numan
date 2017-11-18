#include <stdio.h>
#include <stdlib.h>

double **get_matrix(size_t side) {
    double **matrix = calloc(side, sizeof(double *));
    matrix[0] = calloc(side * (side + 1), sizeof(double));
    size_t i;
    for (i = 1; i < side; i++) {
        matrix[i] = matrix[0] + i * (side+1);
    }
    return matrix;
}

void fill_matrix(double **matrix, size_t side) {
    size_t i, j;
    printf("Enter matrix:\n");
    for (i = 0; i < side; i++) {
        for (j = 0; j < side; j++) {
            scanf("%lf", &matrix[i][j]);
        }
    }
    printf("Enter values vector:\n");
    for (i = 0; i < side; i++) {
        scanf("%lf", &matrix[i][side]);
    }
    return;
}

void __print_matrix(double **matrix, size_t side) {
    size_t i, j;
    for (i = 0; i < side; i++) {
        for (j = 0; j < side+1; ++j) {
            printf("%lf ", matrix[i][j]);
        }
        printf("\n");
    }
    return;
}

int main() {
    size_t side;
    printf("Enter n - side of matrix:\n");
    scanf("%zu", &side);
    double **matrix = get_matrix(side);
    fill_matrix(matrix, side);
    __print_matrix(matrix, side);
    return 0;
}
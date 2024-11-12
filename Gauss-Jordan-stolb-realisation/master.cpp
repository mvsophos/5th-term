//#include <sys/sysinfo.h>
#include <sched.h>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "solve3.h"

int main(int argc, const char* argv[]) {
    int errcode = 0;
    int n, m, r, s;
    // Solving AX = B
    // V1, V2, V3 вспомогательные матрицы
    double *A  = NULL, *B  = NULL, *X  = NULL, *V1 = NULL, *V2 = NULL, *V3 = NULL;
    double r1, r2;
    double t1, t2; // t1 время вычисления решения, t2 время вычисления невязки

    if (!((argc == 5 || argc == 6) &&
        (sscanf(argv[1], "%d", &n) == 1) &&
        (sscanf(argv[2], "%d", &m) == 1) &&
        (sscanf(argv[3], "%d", &r) == 1) &&
        (sscanf(argv[4], "%d", &s) == 1) &&
        (n >= 0) && (m >  0 && m <= n) && (r >= 0 && r <= n)))
            return error(1);

    // выделяем память под матрицу и под правый вектор
    A  = new double[n * n], B  = new double[n * 1], X  = new double[n * 1];
    V1 = new double[m * m], V2 = new double[m * m], V3 = new double[m * m];
    if (A  == NULL || B  == NULL || X  == NULL || V1 == NULL || V2 == NULL || V3 == NULL) {
        delete [] A, delete [] B, delete [] X;
        delete [] V1, delete [] V2, delete [] V3;
        return error(2);
    }

    if (s == 0 && argc == 6)
        fill(A, n, m, 0, argv[5], &errcode);
    else if ((s > 0 && s < 5) && argc == 5)
        fill(A, n, m, s, NULL, NULL);
    else
        errcode = 1;

    if (errcode > 0) {
        delete [] A, delete [] B, delete [] X;
        delete [] V1, delete [] V2, delete [] V3;
        return error(errcode);
    }

    fill_right_part(A, B, n, m);

    if (r > 0) {
        print_matrix(A, n, n, m, r);
        print_matrix(B, 1, n, m, r);
    }

    t1 = clock();
    errcode = solve(n, m, A, B, X, V1, V2, V3); // после выполнения в X и в B записано решение
    t1 = (clock() - t1) / CLOCKS_PER_SEC;
    if (errcode < 0) {
        delete [] A, delete [] B, delete [] X;
        delete [] V1, delete [] V2, delete [] V3;
        return error(5);
    }

    if (s == 0) fill(A, n, m, 0, argv[5], &errcode);
    else fill(A, n, m, s, NULL, NULL);
    fill_right_part(A, B, n, m);
    
    if (s == 0) fill(A, n, m, 0, argv[5], &errcode);
    else fill(A, n, m, s, NULL, NULL);
    
    if (errcode > 0) { // это обработчик ошибок
        delete [] A, delete [] B, delete [] X;
        delete [] V1, delete [] V2, delete [] V3;
        return error(errcode);
    }
    
    t2 = clock();
    r1 = res1(A, B, X, n, m); // теперь правая часть в X (была скопирована перед решением), а в B решение
    for (int i = 0; i < n; i++){ // может переписать этот цикл
        B[i] = (i + 1) % 2;
    }
    r2 = res2(X, B, n);
    t2 = (clock() - t2) / CLOCKS_PER_SEC;

    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", argv[0], 15, r1, r2, t1, t2, s, n, m);
    

    delete [] A, delete [] B, delete [] X;
    delete [] V1, delete [] V2, delete [] V3;
    return 0;
}

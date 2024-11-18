import math

def QR_decomposition(matrix):
    n = len(matrix)
    Q = [[0] * n for _ in range(n)] 
    R = [[0] * n for _ in range(n)]  

    for i in range(n):
        column = [row[i] for row in matrix]
        
        for j in range(i):
            projection = sum(Q[k][j] * column[k] for k in range(n))
            for k in range(n):
                column[k] -= projection * Q[k][j]

        norm = math.sqrt(sum(x**2 for x in column))
        if norm < 1e-12:
            return None, None

        for k in range(n):
            Q[k][i] = column[k] / norm

        for j in range(i, n):
            R[i][j] = sum(Q[k][i] * matrix[k][j] for k in range(n))

    return Q, R

def computing_eigenvalues(matrix, max_iterations=1000, tolerance=1e-15):
    n = len(matrix)
    matrix = [row[:] for row in matrix]

    for _ in range(max_iterations):
        Q, R = QR_decomposition(matrix)
        if Q is None:
            print("Your matrix has linearly dependent rows or columns.")
            return None

    
        matrix = [[sum(R[i][k] * Q[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

        
        off_diag = sum(abs(matrix[i][j]) for i in range(n) for j in range(i))
        if off_diag < tolerance:
            break

    
    return [matrix[i][i] for i in range(n)]

if __name__ == "__main__":
    n = int(input())
    print("Enter elements row by row (space-separated values):")
    matrix = [list(map(float, input().split())) for i in range(n)]

    eigenvalues = computing_eigenvalues(matrix)

    if eigenvalues is None:
        print("Failed.")
    else:
        print("\nEigenvalues:")

        for index, eigenvalue in enumerate(eigenvalues, start=1):
          print("λ%d = %.6f" % (index, eigenvalue))


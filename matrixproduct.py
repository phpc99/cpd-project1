import time
import sys

# Displays the first 10 elements of the first row of the result matrix
def show_result_matrix(matrix, rows, cols):
    print("Result matrix:", end=" ")
    for j in range(min(10, cols)):
        print(matrix[0][j], end=" ")
    print()

# Initializes three matrices
def init_matrix(rows, cols, depth):
    pha = [[1.0 for _ in range(cols)] for _ in range(rows)]         # pha = rows * cols (filled with 1.0)
    phb = [[i + 1 for _ in range(depth)] for i in range(cols)]      # phb = cols * depth (where phb[k][j] = k + 1)
    phc = [[0.0 for _ in range(depth)] for _ in range(rows)]        # phc = rows * depth (initialized with 0.0 to store the result)
    return pha, phb, phc

# on_mult -> Algorithm 1 (rows -> depth -> cols)
def on_mult(rows, cols, depth):
    pha, phb, phc = init_matrix(rows, cols, depth)

    start_time = time.time()

    for i in range(rows):
        for j in range(depth):
            for k in range(cols):
                phc[i][j] += pha[i][k] * phb[k][j]

    elapsed_time = time.time() - start_time
    print(f"Time: {elapsed_time:.6f} seconds; Size: {rows};")

    show_result_matrix(phc, rows, depth)
    return elapsed_time

# on_mult_line -> Algorithm 2 (rows -> cols -> depth)
def on_mult_line(rows, cols, depth):
    pha, phb, phc = init_matrix(rows, cols, depth)

    start_time = time.time()

    for i in range(rows):
        for k in range(cols):
            for j in range(depth):
                phc[i][j] += pha[i][k] * phb[k][j]

    elapsed_time = time.time() - start_time
    print(f"Time: {elapsed_time:.6f} seconds; Size: {rows};")

    show_result_matrix(phc, rows, depth)
    return elapsed_time


def main():
    if len(sys.argv) < 6:
        print("Usage: python script.py <n> <m> <k> <algorithm> <output_file>")
        sys.exit(1)

    rows = int(sys.argv[1])     # n: number of rows
    cols = int(sys.argv[2])     # m: number of cols in pha (same number of rows in phb)
    depth = int(sys.argv[3])    # k: number of cols in phb
    op = int(sys.argv[4])
    output_file = sys.argv[5]
    #block_size = int(sys.argv[6]) if op == 3 and len(sys.argv) >= 7 else 0

    elapsed_time = 0.0

    if op == 1:
        elapsed_time = on_mult(rows, cols, depth)
    elif op == 2:
        elapsed_time = on_mult_line(rows, cols, depth)
    else:
        print("Invalid algorithm choice!")
        sys.exit(1)

    with open(output_file, "a") as outfile:
        outfile.write(f"{op},{rows},{elapsed_time}\n")


if __name__ == "__main__":
    main()

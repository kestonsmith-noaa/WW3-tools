from multiprocessing import Pool
def compute_value(x):
    # Example CPU-heavy or complex calculation
    return x * x

def compute_value2(x):
    # Example CPU-heavy or complex calculation
    return {
        "original": x,
        "square": x ** 2,
        "cube": x ** 3
     }


if __name__ == "__main__":
    #data_array = [1, 2, 3, 4, 5]
    data_array = range(1000)
    with Pool() as pool:
        results = pool.map(compute_value2, data_array)
    squares = [r["square"] for r in results]
    print(results)
    print(squares)


from multiprocessing import Pool, cpu_count

import unafold_python


def run_calc(i):
    print(
        unafold_python.hybrid_min(
            ["TGTGTGCGTACTGCTGCAATAT", "TGTGTGCGTACTGCTGCAATAT"],
            tmin=30,
            tmax=30,
        )
    )


if __name__ == "__main__":
    pool = Pool(cpu_count())
    try:
        test_results_ = pool.map(run_calc, list(range(10)))
    finally:
        pool.close()
        pool.join()

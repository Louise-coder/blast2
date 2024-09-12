import pstats

stats = pstats.Stats("profiling_results.prof")

stats.sort_stats("cumulative").print_stats(10)

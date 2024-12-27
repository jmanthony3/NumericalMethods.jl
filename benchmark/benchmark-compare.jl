using PkgBenchmark
using NumericalMethods

judge(NumericalMethods, "5e183b5", "7931469"; resultfile="results-compare.json")
export_markdown("results-compare.md", readresults("results-compare.json"))
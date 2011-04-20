def run(args):
  assert len(args) == 0
  from libtbx.test_utils import show_diff
  from libtbx.utils import remove_files
  from libtbx import easy_run
  import libtbx.load_env
  import os
  op = os.path
  dd = libtbx.env.dist_path(module_name="Durandal")
  def run_and_check(cmd, pdbs_file, expected_out):
    abs_paths = []
    for rel_path in open(op.join(dd, pdbs_file)).read().splitlines():
      abs_paths.append(op.normpath(op.join(dd, rel_path)))
    print >> open("list_of_pdbs", "w"), "\n".join(abs_paths)
    remove_files("out")
    print cmd
    easy_run.fully_buffered(command=cmd).raise_if_errors()
    filtered_lines = []
    for line in open("out").read().splitlines():
      sw = line.startswith
      if (sw("pole") or sw("members") or sw("min: ") or sw("max: ")):
        filtered_lines.append(line)
    assert not show_diff("\n".join(filtered_lines)+"\n", expected_out)
    print "OK"
  #
  expected_out = """\
pole position centers(5): 2 3 5 7 9
members(3): 1 3 7
members(3): 2 5 9
members(2): 0 6
members(2): 4 8
"""
  run_and_check(
    cmd="durandal.cluster_pdbs -i list_of_pdbs -d 2.0 --brute -o out -m -1",
    pdbs_file="very_few_pdbs",
    expected_out=expected_out)
  #
  run_and_check(
    cmd="durandal.cluster_pdbs -i list_of_pdbs -d 2.0 --smart -o out -m -1",
    pdbs_file="very_few_pdbs",
    expected_out=expected_out)
  #
  run_and_check(
    cmd="durandal.cluster_pdbs -i list_of_pdbs -d 2.0 -s -o out -m -1 -v",
    pdbs_file="very_few_pdbs_e",
    expected_out="""\
pole position centers(5):
min: 11 max: 17 avg: 13.6667 med: 13
min: 12 max: 18 avg: 15.3333 med: 16
members(3):
min: 12 max: 19 avg: 15.3333 med: 15
min: 10 max: 17 avg: 13.6667 med: 14
members(3):
min: 10 max: 16 avg: 13 med: 13
min: 13 max: 19 avg: 16 med: 16
members(2):
min: 14 max: 18 avg: 16 med: 16
min: 11 max: 15 avg: 13 med: 13
members(2):
""")
  #
  print "Done."

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])

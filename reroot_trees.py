import sys
from ete3 import Tree

#Author: Matthew Parks, 2016

#NOTE: change outgroup names as necessary to reflect appropriate taxa

outgroup_names = ["corethron_pennatum_L29A3","Corethron_hystrix_308"]

for line in open(sys.argv[1]):
	t = Tree(line.rstrip())
	outgroups_in_tree = list(set(t.get_leaf_names()).intersection(set(outgroup_names)))
	if len(outgroups_in_tree) > 1:
		ancestor = t.get_common_ancestor(outgroups_in_tree)
		t.set_outgroup(ancestor)
		print t.write()
	elif len(outgroups_in_tree) == 1:
		t.set_outgroup(outgroups_in_tree[0])
		print t.write()
	else:
		continue

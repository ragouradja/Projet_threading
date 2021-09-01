def pdb2str(filename):
	with open(filename, "r") as filin:
		num = -1
		for line in filin:
			if line.startswith("ATOM"):
				items = line.split()
				if num == -1:
					num = items[5]
				else:
					